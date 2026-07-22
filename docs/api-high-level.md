# High-level interface

## Model and atmospheric inputs

```{eval-rst}
.. autoclass:: pratmo.Model
   :members:

.. autoclass:: pratmo.Atmosphere
   :members:

.. autoclass:: pratmo.Box
   :members:
```

## Options

```{eval-rst}
.. autoclass:: pratmo.ChemistryOptions
   :members:

.. autoclass:: pratmo.PhotolysisOptions
   :members:

.. autoclass:: pratmo.DiurnalOptions
   :members:

.. autoclass:: pratmo.CtmOptions
   :members:
```

## Initialization and warnings

```{eval-rst}
.. autofunction:: pratmo.background_mixing_ratios

.. autoclass:: pratmo.PratmoWarning

.. autoclass:: pratmo.ExperimentalFeatureWarning
```

## Quantity helpers

```{eval-rst}
.. autofunction:: pratmo.mixing_ratio
.. autofunction:: pratmo.ppmv
.. autofunction:: pratmo.ppbv
.. autofunction:: pratmo.pptv
.. autofunction:: pratmo.pressure
.. autofunction:: pratmo.temperature
.. autofunction:: pratmo.altitude
.. autofunction:: pratmo.number_density
.. autofunction:: pratmo.surface_area_density
.. autofunction:: pratmo.mixing_ratio_as
.. autofunction:: pratmo.number_density_as
```
