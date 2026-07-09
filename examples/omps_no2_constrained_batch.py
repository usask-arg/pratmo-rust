"""Run PRATMO on OMPS weak-absorber NetCDF profiles constrained by NO2.

The script reads one or more weak-absorber retrieval files, fixes O3 to the
retrieval-assumed profile, scales PRATMO NOy until modeled NO2 matches the
retrieved NO2 at the scan local solar time, and writes a companion NetCDF with
daily mean/noon/midnight/target-LST diagnostics.

Example
-------
python examples/omps_no2_constrained_batch.py \
    /Users/djz828/dev/omps-skretrieval/local_scripts/weak_absorber/output/*.nc \
    --output-dir /private/tmp/pratmo_weak_absorber
"""

from __future__ import annotations

import argparse
import contextlib
import glob
import math
import os
from pathlib import Path
import sys
from typing import Iterable

import numpy as np
import xarray as xr

from pratmo import (
    CustomAtmosphereProfile,
    DiurnBoxSpec,
    DiurnConfig,
    No2ConstrainedDiurnConfig,
    PratmoModel,
)


NL_MAX = 41
NB_MAX = 25
NTIME_MAX = 44

DEFAULT_SPECIES = (
    "no2",
    "no",
    "bro",
    "br",
    "brono2",
    "hobr",
    "io",
    "i",
    "hoi",
    "iono2",
    "o3",
    "hno3",
    "no3",
    "n2o5",
    "clo",
    "clono2",
)

LONG_LIVED_OUTPUTS = {
    "noy": "noy",
    "bry": "brx",
    "brx": "brx",
    "iy": "iodx",
    "iodx": "iodx",
    "cly": "clx",
    "clx": "clx",
    "o3_family": "o3",
}


@contextlib.contextmanager
def suppress_process_output(enabled: bool):
    if not enabled:
        yield
        return

    sys.stdout.flush()
    sys.stderr.flush()
    old_stdout = os.dup(1)
    old_stderr = os.dup(2)
    try:
        with open(os.devnull, "w") as devnull:
            os.dup2(devnull.fileno(), 1)
            os.dup2(devnull.fileno(), 2)
            yield
            sys.stdout.flush()
            sys.stderr.flush()
    finally:
        os.dup2(old_stdout, 1)
        os.dup2(old_stderr, 2)
        os.close(old_stdout)
        os.close(old_stderr)


def parse_csv(value: str) -> list[str]:
    return [part.strip().lower() for part in value.split(",") if part.strip()]


def hours_to_hhmm(hours: float) -> int:
    if not np.isfinite(hours):
        return -1
    total_minutes = int(round((hours % 24.0) * 60.0)) % (24 * 60)
    return (total_minutes // 60) * 100 + total_minutes % 60


def hhmm_to_minutes(hhmm: int) -> int:
    return (hhmm // 100) * 60 + hhmm % 100


def cyclic_time_distance_hhmm(a: int, b: int) -> int:
    da = hhmm_to_minutes(a)
    db = hhmm_to_minutes(b)
    delta = abs(da - db)
    return min(delta, 24 * 60 - delta)


def nearest_time_index(times_hhmm: Iterable[int], target_hhmm: int) -> int:
    times = list(int(t) for t in times_hhmm)
    return min(
        range(len(times)),
        key=lambda idx: cyclic_time_distance_hhmm(times[idx], target_hhmm),
    )


def altitude_km_from_dataset(ds: xr.Dataset) -> np.ndarray:
    alt = np.asarray(ds["altitude"].values, dtype=float)
    if np.nanmax(np.abs(alt)) > 200.0:
        return alt / 1000.0
    return alt


def selected_altitude_indices(
    altitude_km: np.ndarray,
    min_alt_km: float,
    max_alt_km: float,
    max_levels: int,
) -> np.ndarray:
    idx = np.flatnonzero(
        np.isfinite(altitude_km)
        & (altitude_km > min_alt_km)
        & (altitude_km < max_alt_km)
    )
    if idx.size <= max_levels:
        return idx

    sampled = np.rint(np.linspace(0, idx.size - 1, max_levels)).astype(int)
    return idx[np.unique(sampled)]


def julian_day_for_scan(ds: xr.Dataset, scan_index: int, override: int | None) -> int:
    if override is not None:
        return override
    if "time" not in ds:
        return 120

    value = np.asarray(ds["time"].values)[scan_index]
    if np.issubdtype(np.asarray(value).dtype, np.datetime64):
        day = np.datetime64(value, "D")
        year = str(day)[:4]
        first_day = np.datetime64(f"{year}-01-01", "D")
        return int((day - first_day).astype("timedelta64[D]").astype(int)) + 1

    return 120


def output_path_for(input_path: Path, output_dir: Path, suffix: str) -> Path:
    if input_path.name.endswith(".nc"):
        name = input_path.name[:-3] + suffix
    else:
        name = input_path.name + suffix
    return output_dir / name


def init_output_dataset(
    source: xr.Dataset,
    altitude_km: np.ndarray,
    species: list[str],
    long_lived: list[str],
    write_time_series: bool,
) -> xr.Dataset:
    nscan = source.sizes["scanno"]
    nalt = source.sizes["altitude"]

    coords = {
        "scanno": source["scanno"].values if "scanno" in source else np.arange(nscan),
        "altitude": source["altitude"].values,
        "diurn_time_step": np.arange(NTIME_MAX),
    }
    out = xr.Dataset(coords=coords)
    out["altitude_km"] = ("altitude", altitude_km)

    two_d = (nscan, nalt)
    out["pratmo_model_level_selected"] = (
        ("scanno", "altitude"),
        np.zeros(two_d, dtype=bool),
    )
    out["pratmo_observed_no2_target_lst_cm3"] = (
        ("scanno", "altitude"),
        np.full(two_d, np.nan),
    )
    out["pratmo_modeled_no2_target_lst_cm3"] = (
        ("scanno", "altitude"),
        np.full(two_d, np.nan),
    )
    out["pratmo_noy_scale"] = (("scanno", "altitude"), np.full(two_d, np.nan))
    out["pratmo_no2_noy_ratio_target_lst"] = (
        ("scanno", "altitude"),
        np.full(two_d, np.nan),
    )
    out["pratmo_no2_modeled_observed_ratio_target_lst"] = (
        ("scanno", "altitude"),
        np.full(two_d, np.nan),
    )
    out["pratmo_target_hhmm"] = ("scanno", np.full(nscan, -1, dtype=np.int32))
    out["pratmo_time_hhmm"] = (
        ("scanno", "diurn_time_step"),
        np.full((nscan, NTIME_MAX), -1, dtype=np.int32),
    )

    for name in species:
        out[f"pratmo_{name}_daily_cm3"] = (
            ("scanno", "altitude"),
            np.full(two_d, np.nan),
        )
        out[f"pratmo_{name}_noon_cm3"] = (
            ("scanno", "altitude"),
            np.full(two_d, np.nan),
        )
        out[f"pratmo_{name}_midnight_cm3"] = (
            ("scanno", "altitude"),
            np.full(two_d, np.nan),
        )
        out[f"pratmo_{name}_target_lst_cm3"] = (
            ("scanno", "altitude"),
            np.full(two_d, np.nan),
        )
        if write_time_series:
            out[f"pratmo_{name}_time_cm3"] = (
                ("scanno", "altitude", "diurn_time_step"),
                np.full((nscan, nalt, NTIME_MAX), np.nan),
            )

    for name in long_lived:
        out[f"pratmo_{name}_daily_vmr"] = (
            ("scanno", "altitude"),
            np.full(two_d, np.nan),
        )
        out[f"pratmo_{name}_daily_cm3"] = (
            ("scanno", "altitude"),
            np.full(two_d, np.nan),
        )

    return out


def finite_positive_profile(values: np.ndarray) -> np.ndarray:
    return np.isfinite(values) & (values > 0.0)


def process_scan(
    model: PratmoModel,
    ds: xr.Dataset,
    out: xr.Dataset,
    scan_index: int,
    altitude_indices: np.ndarray,
    species: list[str],
    long_lived: list[str],
    args: argparse.Namespace,
) -> None:
    pressure_mb = np.asarray(ds["pressure_pa"].isel(scanno=scan_index).values, dtype=float) / 100.0
    temperature_k = np.asarray(ds["temperature_k"].isel(scanno=scan_index).values, dtype=float)
    o3_vmr = np.asarray(ds["assumed_o3_vmr"].isel(scanno=scan_index).values, dtype=float)
    no2_vmr = np.asarray(ds["no2_vmr"].isel(scanno=scan_index).values, dtype=float)
    air_m3 = np.asarray(ds["background_numberdensity"].isel(scanno=scan_index).values, dtype=float)
    altitude_km = np.asarray(out["altitude_km"].values, dtype=float)

    atmosphere_ok = (
        finite_positive_profile(pressure_mb)
        & finite_positive_profile(temperature_k)
        & np.isfinite(o3_vmr)
        & (o3_vmr >= 0.0)
    )
    atm_indices = np.array([idx for idx in altitude_indices if atmosphere_ok[idx]], dtype=int)
    if atm_indices.size == 0:
        return

    # PRATMO expects pressure to decrease with altitude and altitude to increase.
    order = np.argsort(altitude_km[atm_indices])
    atm_indices = atm_indices[order]

    rel_lookup = {int(original_idx): rel for rel, original_idx in enumerate(atm_indices)}
    no2_cm3 = no2_vmr * air_m3 / 1.0e6
    obs_ok = np.isfinite(no2_cm3) & (no2_cm3 > 0.0)
    box_original_indices = np.array(
        [idx for idx in atm_indices if obs_ok[idx]],
        dtype=int,
    )
    if box_original_indices.size == 0:
        return

    latitude = float(ds["latitude"].isel(scanno=scan_index).values)
    local_solar_time = float(ds["local_solar_time_hours"].isel(scanno=scan_index).values)
    target_hhmm = hours_to_hhmm(local_solar_time)
    if not np.isfinite(latitude) or target_hhmm < 0:
        return

    out["pratmo_target_hhmm"].values[scan_index] = target_hhmm

    atmosphere = CustomAtmosphereProfile(
        pressure_mb=pressure_mb[atm_indices].tolist(),
        temperature_k=temperature_k[atm_indices].tolist(),
        o3=o3_vmr[atm_indices].tolist(),
        o3_kind="mixing_ratio",
        altitude_km=altitude_km[atm_indices].tolist(),
    )

    julian_day = julian_day_for_scan(ds, scan_index, args.julian_day)

    for start in range(0, box_original_indices.size, args.max_boxes_per_run):
        chunk_original = box_original_indices[start : start + args.max_boxes_per_run]
        relative_levels = [rel_lookup[int(idx)] + 1 for idx in chunk_original]
        boxes = [
            DiurnBoxSpec(altitude_level=int(level), albedo=args.albedo)
            for level in relative_levels
        ]

        diurn = DiurnConfig(
            latitude_deg=latitude,
            julian_day=julian_day,
            integration_days=args.integration_days,
            boxes=boxes,
            bromine=True,
            iodine=not args.disable_iodine,
            parallel_boxes=True,
            atmosphere=atmosphere,
        )
        constrained = No2ConstrainedDiurnConfig(
            diurn=diurn,
            observed_no2_cm3=no2_cm3[chunk_original].tolist(),
            target_hhmm=target_hhmm,
            iterations=args.iterations,
        )
        with suppress_process_output(not args.verbose_pratmo):
            result = model.run_diurn_no2_constrained(constrained)
        output = result.output
        if not output.time_series:
            continue

        times_hhmm = np.array(
            [step.time_hhmm for step in output.time_series[0].steps],
            dtype=np.int32,
        )
        ntimes = min(times_hhmm.size, NTIME_MAX)
        out["pratmo_time_hhmm"].values[scan_index, :ntimes] = times_hhmm[:ntimes]

        noon_idx = nearest_time_index(times_hhmm, 1200)
        midnight_idx = nearest_time_index(times_hhmm, 0)
        target_idx = nearest_time_index(times_hhmm, target_hhmm)

        species_grids = {name: np.asarray(output.species_grid(name)) for name in species}

        for box_index, original_idx in enumerate(chunk_original):
            out["pratmo_model_level_selected"].values[scan_index, original_idx] = True
            out["pratmo_observed_no2_target_lst_cm3"].values[scan_index, original_idx] = no2_cm3[
                original_idx
            ]
            out["pratmo_modeled_no2_target_lst_cm3"].values[scan_index, original_idx] = (
                result.modeled_no2_cm3[box_index]
            )
            if no2_cm3[original_idx] > 0.0:
                out["pratmo_no2_modeled_observed_ratio_target_lst"].values[
                    scan_index, original_idx
                ] = result.modeled_no2_cm3[box_index] / no2_cm3[original_idx]
            out["pratmo_noy_scale"].values[scan_index, original_idx] = result.noy_scale[box_index]

            snap = output.boxes[box_index]
            implicit = snap.implicit.to_dict()
            long_lived_dict = snap.long_lived.to_dict()
            air_cm3 = snap.air_density_cm3

            for name in species:
                grid = species_grids[name]
                out[f"pratmo_{name}_daily_cm3"].values[scan_index, original_idx] = implicit[name]
                out[f"pratmo_{name}_noon_cm3"].values[scan_index, original_idx] = grid[
                    box_index, noon_idx
                ]
                out[f"pratmo_{name}_midnight_cm3"].values[scan_index, original_idx] = grid[
                    box_index, midnight_idx
                ]
                out[f"pratmo_{name}_target_lst_cm3"].values[scan_index, original_idx] = grid[
                    box_index, target_idx
                ]
                if args.write_time_series:
                    out[f"pratmo_{name}_time_cm3"].values[
                        scan_index, original_idx, :ntimes
                    ] = grid[box_index, :ntimes]

            for output_name in long_lived:
                field = LONG_LIVED_OUTPUTS[output_name]
                vmr = long_lived_dict[field]
                out[f"pratmo_{output_name}_daily_vmr"].values[scan_index, original_idx] = vmr
                out[f"pratmo_{output_name}_daily_cm3"].values[scan_index, original_idx] = (
                    vmr * air_cm3
                )

            noy_cm3 = long_lived_dict["noy"] * air_cm3
            if noy_cm3 > 0.0 and math.isfinite(noy_cm3):
                out["pratmo_no2_noy_ratio_target_lst"].values[scan_index, original_idx] = (
                    result.modeled_no2_cm3[box_index] / noy_cm3
                )


def process_file(input_path: Path, output_path: Path, args: argparse.Namespace) -> None:
    open_kwargs = {}
    if args.engine:
        open_kwargs["engine"] = args.engine
    ds = xr.open_dataset(input_path, decode_times=True, **open_kwargs)
    try:
        altitude_km = altitude_km_from_dataset(ds)
        altitude_indices = selected_altitude_indices(
            altitude_km,
            args.min_alt_km,
            args.max_alt_km,
            args.max_atmosphere_levels,
        )
        species = parse_csv(args.species)
        long_lived = parse_csv(args.long_lived)
        unknown_long_lived = sorted(set(long_lived) - set(LONG_LIVED_OUTPUTS))
        if unknown_long_lived:
            raise ValueError(f"Unknown long-lived outputs: {', '.join(unknown_long_lived)}")

        out = init_output_dataset(
            ds,
            altitude_km,
            species,
            long_lived,
            write_time_series=args.write_time_series,
        )
        out.attrs.update(
            {
                "title": "PRATMO NO2-constrained diurnal diagnostics",
                "source_file": str(input_path),
                "pratmo_min_alt_km": args.min_alt_km,
                "pratmo_max_alt_km": args.max_alt_km,
                "pratmo_iterations": args.iterations,
                "pratmo_integration_days": args.integration_days,
                "pratmo_note": (
                    "NO2 VMR was converted to cm-3 using "
                    "background_numberdensity in m-3."
                ),
            }
        )

        model = PratmoModel.with_defaults()
        nscan = ds.sizes["scanno"]
        if args.scan_limit is not None:
            nscan = min(nscan, args.scan_limit)
        for scan_index in range(nscan):
            print(f"  scan {scan_index + 1}/{nscan}")
            process_scan(
                model,
                ds,
                out,
                scan_index,
                altitude_indices,
                species,
                long_lived,
                args,
            )

        output_path.parent.mkdir(parents=True, exist_ok=True)
        out.to_netcdf(output_path)
    finally:
        ds.close()


def expand_inputs(patterns: list[str]) -> list[Path]:
    paths: list[Path] = []
    for pattern in patterns:
        expanded = sorted(glob.glob(pattern)) if any(ch in pattern for ch in "*?[") else []
        if expanded:
            paths.extend(Path(p) for p in expanded)
        else:
            paths.append(Path(pattern))
    return paths


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", help="Input weak-absorber NetCDF files or globs.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for companion PRATMO NetCDF files. Defaults to each input directory.",
    )
    parser.add_argument(
        "--suffix",
        default="_pratmo_no2constrained.nc",
        help="Suffix appended to each input stem.",
    )
    parser.add_argument("--engine", default=None, help="Optional xarray NetCDF engine.")
    parser.add_argument(
        "--min-alt-km",
        type=float,
        default=8.0,
        help="Lower altitude bound in km; selected levels must be greater than this.",
    )
    parser.add_argument(
        "--max-alt-km",
        type=float,
        default=40.0,
        help="Upper altitude bound in km; selected levels must be less than this.",
    )
    parser.add_argument("--max-atmosphere-levels", type=int, default=NL_MAX)
    parser.add_argument(
        "--max-boxes-per-run",
        type=int,
        default=1,
        help=(
            "Number of constrained altitude boxes solved together. The safe default is 1; "
            "larger chunks are faster but can perturb the NO2-constrained NOy iteration."
        ),
    )
    parser.add_argument("--iterations", type=int, default=4)
    parser.add_argument("--integration-days", type=int, default=3)
    parser.add_argument("--julian-day", type=int, default=None)
    parser.add_argument("--albedo", type=float, default=0.0)
    parser.add_argument("--disable-iodine", action="store_true")
    parser.add_argument(
        "--verbose-pratmo",
        action="store_true",
        help="Do not suppress PRATMO's internal stdout/stderr messages.",
    )
    parser.add_argument(
        "--species",
        default=",".join(DEFAULT_SPECIES),
        help="Comma-separated implicit species to write.",
    )
    parser.add_argument(
        "--long-lived",
        default="noy,bry,iy,cly",
        help="Comma-separated long-lived families to write.",
    )
    parser.add_argument(
        "--no-time-series",
        dest="write_time_series",
        action="store_false",
        help="Skip full diurnal time-series variables.",
    )
    parser.set_defaults(write_time_series=True)
    parser.add_argument("--scan-limit", type=int, default=None)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.max_atmosphere_levels = max(1, min(args.max_atmosphere_levels, NL_MAX))
    args.max_boxes_per_run = max(1, min(args.max_boxes_per_run, NB_MAX))

    for input_path in expand_inputs(args.inputs):
        if not input_path.exists():
            raise FileNotFoundError(input_path)
        output_dir = args.output_dir if args.output_dir is not None else input_path.parent
        output_path = output_path_for(input_path, output_dir, args.suffix)
        print(f"Processing {input_path}")
        process_file(input_path, output_path, args)
        print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
