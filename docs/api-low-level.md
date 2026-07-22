# Lower-level compatibility API

These classes expose the native Rust-extension configuration directly. Values
use canonical units and fewer convenience checks. Prefer {class}`pratmo.Model`
for new code.

## Model and configuration

```{eval-rst}
.. autoclass:: pratmo.PratmoModel
   :members:

.. autoclass:: pratmo.DiurnConfig
   :members:

.. autoclass:: pratmo.CtmConfig
   :members:

.. autoclass:: pratmo.DiurnBoxSpec
   :members:

.. autoclass:: pratmo.CtmBoxSpec
   :members:

.. autoclass:: pratmo.CustomAtmosphereProfile
   :members:

.. autoclass:: pratmo.LongLivedMixingRatios
   :members:

.. autoclass:: pratmo.No2ConstrainedDiurnConfig
   :members:
```

## Detailed records

```{eval-rst}
.. autoclass:: pratmo.BoxSnapshot
   :members:

.. autoclass:: pratmo.DiurnBoxTimeSeries
   :members:

.. autoclass:: pratmo.DiurnTimeStep
   :members:

.. autoclass:: pratmo.ImplicitSpecies
   :members:

.. autoclass:: pratmo.JValues
   :members:
```

## Discoverable field names

```{eval-rst}
.. autodata:: pratmo.IMPLICIT_SPECIES_NAMES
.. autodata:: pratmo.LONG_LIVED_NAMES
.. autodata:: pratmo.JVALUE_NAMES
```
