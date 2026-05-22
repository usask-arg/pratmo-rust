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
    BoxSnapshot,
    CtmBoxSpec,
    CtmConfig,
    CtmOutput,
    Diagnostics,
    DiurnBoxSpec,
    DiurnBoxTimeSeries,
    DiurnConfig,
    DiurnOutput,
    DiurnTimeStep,
    ImplicitSpecies,
    JValues,
    LongLivedMixingRatios,
    PratmoModel,
)

__all__ = [
    "PratmoModel",
    "DiurnConfig",
    "CtmConfig",
    "DiurnBoxSpec",
    "CtmBoxSpec",
    "DiurnOutput",
    "CtmOutput",
    "BoxSnapshot",
    "DiurnBoxTimeSeries",
    "DiurnTimeStep",
    "ImplicitSpecies",
    "LongLivedMixingRatios",
    "JValues",
    "Diagnostics",
]
