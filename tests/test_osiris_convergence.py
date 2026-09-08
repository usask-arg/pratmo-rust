"""Real OSIRIS baseline inputs; no retrieval packages or external data needed."""

from dataclasses import replace
import json
from pathlib import Path

import numpy as np
import pytest

import pratmo


FIXTURES = Path(__file__).parent / "fixtures" / "osiris"
SLOW_SPECIES = ("hbr", "ch3o2h", "n2o5", "hno3", "hno4", "h2o2", "hcl")
FAMILIES = {
    "noy": {"no": 1, "no2": 1, "no3": 1, "n2o5": 2, "hno3": 1,
            "hno2": 1, "hno4": 1, "clono2": 1, "brono2": 1},
    "clx": {"hcl": 1, "cl": 1, "cl2": 2, "clo": 1, "clono2": 1,
            "hocl": 1, "oclo": 1, "cl2o2": 2, "brcl": 1},
    "brx": {"bro": 1, "br": 1, "hbr": 1, "brono2": 1, "hobr": 1, "brcl": 1},
}


def scan_arguments(scan, *, all_boxes=False, parallel=True, days=20):
    saved = json.loads((FIXTURES / f"o{scan}-pratmo-inputs.json").read_text())
    indices = range(23) if all_boxes else [5]
    return dict(
        latitude=saved["latitude"],
        day=saved["day"],
        atmosphere=pratmo.Atmosphere(**saved["atmosphere"]),
        boxes=[pratmo.Box(**saved["boxes"][i]) for i in indices],
        initial_mixing_ratios=[
            pratmo.LongLivedMixingRatios(**saved["initial_mixing_ratios"][i])
            for i in indices
        ],
        chemistry=pratmo.ChemistryOptions(**saved["chemistry"]),
        photolysis=pratmo.PhotolysisOptions(**saved["photolysis"]),
        options=replace(pratmo.DiurnalOptions(**saved["options"]),
                        parallel_boxes=parallel, integration_days=days),
    )


def assert_converged_chemistry(out):
    diag = out.diagnostics
    assert diag.newraf_nonconvergence_count == 0
    assert diag.rafday_nonconvergence_count == 0
    assert diag.rafday_warnings == []
    assert np.isfinite(diag.rafday_max_final_relative_correction)
    assert diag.rafday_max_final_relative_correction < 3.0e-5
    assert 1 <= diag.rafday_max_correction_iterations <= 9

    # Validate the returned orbit independently of the solver's stopping flag.
    # For these slow species, the fractional noon-to-noon change is the
    # integrated daily P-L divided by the starting density.
    for name in SLOW_SPECIES:
        grid = out.species_grid(name)
        assert np.all(grid > 0.0), name
        np.testing.assert_allclose(grid[:, -1] / grid[:, 0], 1.0, rtol=3.0e-5)
    for name in pratmo.IMPLICIT_SPECIES_NAMES:
        grid = out.species_grid(name)
        assert np.all(np.isfinite(grid)), name
        assert np.all(grid >= 0.0), name
        snapshot = out.species_profile(name)
        assert np.all(np.isfinite(snapshot)) and np.all(snapshot >= 0.0), name
    for family, members in FAMILIES.items():
        target = out.long_lived_profile(family) * out.air_density_cm3
        snapshot = sum(weight * out.species_profile(name) for name, weight in members.items())
        cycle = sum(weight * out.species_grid(name) for name, weight in members.items())
        np.testing.assert_allclose(snapshot / target, 1.0, rtol=1.0e-8)
        np.testing.assert_allclose(cycle / target[:, None], 1.0, rtol=1.0e-8)


@pytest.mark.parametrize("scan", [37375021, 37375022])
@pytest.mark.parametrize("all_boxes", [False, True], ids=["17.5km", "full-column"])
@pytest.mark.parametrize("parallel", [False, True], ids=["serial", "parallel"])
def test_osiris_scan_converges(scan, all_boxes, parallel):
    out = pratmo.Model().diurnal(**scan_arguments(scan, all_boxes=all_boxes, parallel=parallel))
    assert len(out.boxes) == (23 if all_boxes else 1)
    assert_converged_chemistry(out)


@pytest.mark.parametrize("days", [3, 40, 100])
def test_failing_box_does_not_require_more_integration_days(days):
    out = pratmo.Model().diurnal(**scan_arguments(37375021, days=days))
    assert_converged_chemistry(out)


def test_neighboring_scan_retains_released_solution():
    expected = json.loads((FIXTURES / "o37375022-v0.3.0-species.json").read_text())
    out = pratmo.Model().diurnal(**scan_arguments(37375022, all_boxes=True))
    for name, values in expected.items():
        np.testing.assert_allclose(out.species_profile(name), values, rtol=5.0e-5, atol=1.0e-20,
                                   err_msg=name)


@pytest.mark.parametrize("failure", ["iteration-limit", "invalid-jacobian"])
@pytest.mark.parametrize("parallel", [False, True], ids=["serial", "parallel"])
def test_solver_failure_warns_and_retains_usable_column(tmp_path, failure, parallel):
    # Exercise the actual public recovery path using the legacy numerical
    # controls, without changing the captured atmosphere or composition.
    import shutil

    source = Path(__file__).resolve().parents[1] / "pratmo-core" / "data"
    for name in ("fort01.x", "fort02.x", "fort10_cam06.x", "fort11_jpl09.x",
                 "fort13.x", "fort14.x", "J_H2O_SZA0.dat"):
        shutil.copyfile(source / name, tmp_path / name)
    config = tmp_path / "fort01.x"
    lines = config.read_text().splitlines()
    for i, line in enumerate(lines):
        if failure == "iteration-limit" and line.startswith("XRF/RLX/LB"):
            lines[i] = line[:10] + f"{0:5d}" + line[15:]
        elif failure == "invalid-jacobian" and line.startswith("RALERR"):
            values = line[10:].split()
            values[5] = "0.0"
            lines[i] = line[:10] + " ".join(values)
    config.write_text("\n".join(lines) + "\n")

    with pytest.warns(RuntimeWarning, match="Returning the last valid diurnal cycle") as caught:
        out = pratmo.Model(tmp_path).diurnal(
            **scan_arguments(37375021, all_boxes=True, parallel=parallel)
        )
    assert len(caught) == 1  # Aggregate boxes into one warning per public run.
    assert len(out.boxes) == 23
    assert out.diagnostics.rafday_nonconvergence_count == 23
    assert len(out.diagnostics.rafday_warnings) == 23
    assert "box 6 at 17.500 km" in out.diagnostics.rafday_warnings[5]
    reason = "iteration limit" if failure == "iteration-limit" else "invalid Jacobian perturbation"
    assert reason in out.diagnostics.rafday_warnings[5]
    assert out.diagnostics.newraf_nonconvergence_count == 0
    for name in pratmo.IMPLICIT_SPECIES_NAMES:
        grid = out.species_grid(name)
        assert np.all(np.isfinite(grid)) and np.all(grid >= 0.0), name
    for family, members in FAMILIES.items():
        target = out.long_lived_profile(family) * out.air_density_cm3
        snapshot = sum(weight * out.species_profile(name) for name, weight in members.items())
        cycle = sum(weight * out.species_grid(name) for name, weight in members.items())
        np.testing.assert_allclose(snapshot / target, 1.0, rtol=1.0e-8)
        np.testing.assert_allclose(cycle / target[:, None], 1.0, rtol=1.0e-8)
    # A restored output must be labelled honestly, not reported as equilibrium.
    assert any(np.max(np.abs(out.species_grid(name)[:, -1] /
                             out.species_grid(name)[:, 0] - 1.0)) > 3.0e-5
               for name in SLOW_SPECIES)
