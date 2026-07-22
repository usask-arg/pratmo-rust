"""Compare Rust PRATMO with a legacy C++ PRATMO/OSIRIS diurnal profile.

This reproduces the C++ wrapper configuration: an 81-level radiative column,
one chemistry box, monthly PRATMO N2O/NOy/ASA climatologies, measured O3 and
ECMWF state, fixed-mu time steps, endpoint convergence, iodine off, sulfate
heterogeneous chemistry on, file albedo, and no rainout.
"""

from __future__ import annotations

import argparse
from datetime import date
from pathlib import Path

import matplotlib
import numpy as np
import xarray as xr

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pratmo import (
    CustomAtmosphereProfile,
    DiurnBoxSpec,
    DiurnConfig,
    PratmoClimatology,
    PratmoModel,
)


TARGET_ALTITUDE_KM = 26.5
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def vertical_interp(
    source_altitude: np.ndarray, values: np.ndarray, target_altitude: np.ndarray
) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim == 1:
        return np.interp(target_altitude, source_altitude, array)
    if array.ndim == 2 and array.shape[0] == len(source_altitude):
        return np.column_stack(
            [
                np.interp(target_altitude, source_altitude, array[:, column])
                for column in range(array.shape[1])
            ]
        )
    raise ValueError("values must have altitude as their first dimension")


def linear_extrapolate(
    source_altitude: np.ndarray, values: np.ndarray, target_altitude: np.ndarray
) -> np.ndarray:
    """Interpolate inside the measured range and linearly extend both ends."""
    source = np.asarray(source_altitude, dtype=float)
    values = np.asarray(values, dtype=float)
    target = np.asarray(target_altitude, dtype=float)
    result = np.interp(target, source, values)
    below = target < source[0]
    above = target > source[-1]
    result[below] = values[0] + (target[below] - source[0]) * (
        values[1] - values[0]
    ) / (source[1] - source[0])
    result[above] = values[-1] + (target[above] - source[-1]) * (
        values[-1] - values[-2]
    ) / (source[-1] - source[-2])
    return result


def normalized_rmse(model: np.ndarray, reference: np.ndarray) -> float:
    return 100.0 * np.sqrt(np.mean((model - reference) ** 2)) / np.sqrt(
        np.mean(reference**2)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input", type=Path, help="archived C++ PRATMO/OSIRIS NetCDF file"
    )
    parser.add_argument(
        "--output", type=Path, default=Path("pratmo_osiris_cpp_reconciled.png")
    )
    parser.add_argument(
        "--climatology-dir", type=Path, default=REPOSITORY_ROOT / "fortran"
    )
    parser.add_argument("--integration-days", type=int, default=150)
    parser.add_argument(
        "--no-radiative-aerosol",
        action="store_true",
        help="leave aerosol out of the photolysis optical depth (ablation run)",
    )
    args = parser.parse_args()

    source = xr.open_dataset(args.input).load()
    source_altitude = np.asarray(source.Altitude, dtype=float) / 1000.0

    # Preserve the exact C++ radiative shells. The chemistry box is off-grid
    # and obtains a wavelength-by-wavelength interpolation of the 26/27-km
    # shell fluxes.
    atmosphere_altitude = np.linspace(0.0, 80.0, 81)
    target_level = int(np.searchsorted(atmosphere_altitude, TARGET_ALTITUDE_KM) - 1)
    box_altitude = np.array([TARGET_ALTITUDE_KM])
    target_box = 0

    # The archived file contains the ECMWF column only over the retrieval
    # range. Extend it to the C++ radiative shells; pressure is log-linear,
    # temperature is linear, and the user-defined O3 profile clamps at its ends.
    pressure_mb = np.exp(
        linear_extrapolate(
            source_altitude,
            np.log(np.asarray(source.pressure) * 0.01),
            atmosphere_altitude,
        )
    )
    temperature_k = linear_extrapolate(
        source_altitude, source.temperature, atmosphere_altitude
    )
    o3_vmr = vertical_interp(source_altitude, source.o3, atmosphere_altitude)

    climatology = PratmoClimatology(args.climatology_dir)
    radiative_climatology = climatology.sample(
        float(source.latitude), date(2008, 7, 2), atmosphere_altitude
    )
    atmosphere = CustomAtmosphereProfile(
        pressure_mb=pressure_mb,
        temperature_k=temperature_k,
        altitude_km=atmosphere_altitude,
        o3=o3_vmr,
        o3_kind="mixing_ratio",
        aerosol_surface_area_um2_cm3=radiative_climatology.aerosol_surface_area_um2_cm3,
    )
    climatology_profile = climatology.sample(
        float(source.latitude), date(2008, 7, 2), box_altitude
    )
    observed_box_o3 = vertical_interp(source_altitude, source.o3, box_altitude)
    boxes = [
        DiurnBoxSpec(
            altitude_level=int(target_level + 1),
            altitude_km=TARGET_ALTITUDE_KM,
            aerosol_surface_area_um2_cm3=float(aerosol_area),
            sea_salt_surface_area_um2_cm3=0.0,
        )
        for aerosol_area in climatology_profile.aerosol_surface_area_um2_cm3
    ]
    config = DiurnConfig(
        latitude_deg=float(source.latitude),
        julian_day=184,
        integration_days=args.integration_days,
        atmosphere=atmosphere,
        boxes=boxes,
        bromine=True,
        iodine=False,
        parallel_boxes=False,
        cpp_compatibility=True,
        # Integrate at the exact 49 half-hour timestamps stored in the C++
        # archive, rather than interpolating that sharp twilight transition
        # onto PRATMO's generated fixed-mu grid.
        elapsed_time_hours=(np.asarray(source.LST, dtype=float) + 12.0).tolist(),
        surface_albedo=float(source.albedo),
        heterogeneous_chemistry=True,
        radiative_aerosol=not args.no_radiative_aerosol,
        initial_mixing_ratios=climatology_profile.initial_mixing_ratios(
            o3=observed_box_o3
        ),
    )

    model = PratmoModel.with_defaults()
    corrected = model.run_diurn(config)

    elapsed_hours = np.asarray(corrected.elapsed_seconds) / 3600.0
    reference_hours = np.asarray(source.LST, dtype=float) + 12.0

    reference_no2_native = vertical_interp(
        source_altitude, source.NO2, np.array([TARGET_ALTITUDE_KM])
    )[0]
    reference_no_native = vertical_interp(
        source_altitude, source.NO, np.array([TARGET_ALTITUDE_KM])
    )[0]
    reference_no2 = np.interp(elapsed_hours, reference_hours, reference_no2_native)
    reference_no = np.interp(elapsed_hours, reference_hours, reference_no_native)

    corrected_no2 = np.asarray(corrected.species_grid("no2"))[target_box]
    corrected_no = np.asarray(corrected.species_grid("no"))[target_box]

    series = [
        ("NO$_2$", reference_no2, corrected_no2, "molecules cm$^{-3}$"),
        ("NO", reference_no, corrected_no, "molecules cm$^{-3}$"),
        (
            "NO/NO$_2$",
            reference_no / reference_no2,
            corrected_no / corrected_no2,
            "ratio",
        ),
    ]

    fig, axes = plt.subplots(3, 1, figsize=(7.8, 8.8), sharex=True)
    for axis, (name, reference, model_values, unit) in zip(axes, series):
        axis.plot(
            elapsed_hours,
            reference,
            color="black",
            linestyle="--",
            label="Archived C++ PRATMO",
        )
        axis.plot(elapsed_hours, model_values, label="Rust, matched C++ grids")
        axis.set_ylabel(f"{name}\n{unit}")
        axis.grid(alpha=0.25)
        amplitude_ratio = np.dot(model_values, reference) / np.dot(reference, reference)
        print(
            f"{name}: NRMSE={normalized_rmse(model_values, reference):.2f}% "
            f"least-squares amplitude ratio={amplitude_ratio:.4f}"
        )

    event_hour = float(source.Event_LST) + 12.0
    event_no2 = float(
        vertical_interp(
            source_altitude, source.OSIRIS_NO2, np.array([TARGET_ALTITUDE_KM])
        )[0]
    )
    axes[0].scatter(
        [event_hour],
        [event_no2],
        marker="o",
        s=28,
        color="tab:red",
        label="OSIRIS event NO$_2$",
        zorder=4,
    )

    axes[0].legend()
    axes[0].set_title(
        f"2008-07-02, {TARGET_ALTITUDE_KM:.1f} km, {float(source.latitude):.2f}°\n"
        "C++ 0–80 km radiative shells · flux interpolated to 26.5 km\n"
        "exact archived 0.5-hour integration grid\n"
        f"albedo={float(source.albedo):.3f} · hetero on · aerosol RT "
        f"{'off' if args.no_radiative_aerosol else 'on'} · iodine/rainout off",
        fontsize=11,
    )
    axes[-1].set_xlabel("Elapsed hours from local noon")
    axes[-1].set_xticks(
        [0, 6, 12, 18, 24], ["12:00", "18:00", "00:00", "06:00", "12:00"]
    )
    fig.tight_layout()
    fig.savefig(args.output, dpi=300)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
