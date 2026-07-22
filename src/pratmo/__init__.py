"""Python interface for the PRATMO stratospheric photochemical box model.

Typical usage::

    from pratmo import Model

    model = Model()
    out = model.diurnal(latitude=0.0, day="2026-04-30")
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
from pratmo.interface import (
    Atmosphere,
    Box,
    ChemistryOptions,
    CtmOptions,
    DiurnalOptions,
    ExperimentalFeatureWarning,
    Model,
    PhotolysisOptions,
    PratmoWarning,
    background_mixing_ratios,
)
from pratmo.quantities import (
    altitude,
    mixing_ratio,
    mixing_ratio_as,
    number_density,
    number_density_as,
    ppbv,
    ppmv,
    pptv,
    pressure,
    surface_area_density,
    temperature,
)

__all__ = [
    "IMPLICIT_SPECIES_NAMES",
    "LONG_LIVED_NAMES",
    "JVALUE_NAMES",
    "Model",
    "Atmosphere",
    "Box",
    "ChemistryOptions",
    "PhotolysisOptions",
    "DiurnalOptions",
    "CtmOptions",
    "PratmoWarning",
    "ExperimentalFeatureWarning",
    "background_mixing_ratios",
    "mixing_ratio",
    "mixing_ratio_as",
    "ppmv",
    "ppbv",
    "pptv",
    "pressure",
    "temperature",
    "altitude",
    "number_density",
    "number_density_as",
    "surface_area_density",
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
