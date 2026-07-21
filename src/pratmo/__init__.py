"""
pratmo — Python bindings for the PRATMO stratospheric photochemical box model.

Typical usage::

    from pratmo import PratmoModel, DiurnConfig, DiurnBoxSpec

    model = PratmoModel.with_defaults()
    cfg = DiurnConfig(
        latitude_deg=0.0,
        julian_day=120,
        integration_days=20,
        boxes=[DiurnBoxSpec(altitude_level=25)],
    )
    out = model.run_diurn(cfg)
    o3 = out.species_grid("o3")   # numpy array, shape (n_boxes, n_timesteps)
    print(out.boxes[0].implicit.o3)
"""

from pratmo._pratmo import (
    IMPLICIT_SPECIES_NAMES,
    JVALUE_NAMES,
    LONG_LIVED_NAMES,
    BoxSnapshot,
    CtmBoxSpec,
    CtmConfig,
    CtmOutput,
    CustomAtmosphereProfile,
    Diagnostics,
    DiurnBoxSpec,
    DiurnBoxTimeSeries,
    DiurnConfig,
    DiurnOutput,
    DiurnTimeStep,
    ImplicitSpecies,
    JValues,
    LongLivedMixingRatios,
    No2ConstrainedDiurnConfig,
    No2ConstrainedDiurnOutput,
    PratmoModel,
)
from pratmo.climatology import PratmoClimatology, PratmoClimatologyProfile

__all__ = [
    "IMPLICIT_SPECIES_NAMES",
    "LONG_LIVED_NAMES",
    "JVALUE_NAMES",
    "PratmoModel",
    "DiurnConfig",
    "CtmConfig",
    "DiurnBoxSpec",
    "CtmBoxSpec",
    "CustomAtmosphereProfile",
    "No2ConstrainedDiurnConfig",
    "No2ConstrainedDiurnOutput",
    "DiurnOutput",
    "CtmOutput",
    "BoxSnapshot",
    "DiurnBoxTimeSeries",
    "DiurnTimeStep",
    "ImplicitSpecies",
    "LongLivedMixingRatios",
    "JValues",
    "Diagnostics",
    "PratmoClimatology",
    "PratmoClimatologyProfile",
]
