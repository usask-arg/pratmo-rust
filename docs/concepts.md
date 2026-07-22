# Model concepts

PRATMO is a photochemical box model. It follows chemistry in one or more
vertical boxes, but it does not move air or chemical tracers between them.

## A box

Each chemistry box has:

- pressure, temperature, altitude, and air number density;
- prescribed long-lived gases and chemical families;
- 40 short-lived or implicit species solved by the chemical mechanism;
- photolysis rates calculated from the surrounding radiative atmosphere; and
- optional sulfate-aerosol and sea-salt surface area for heterogeneous
  reactions.

Boxes in the same call share the latitude, date, and radiation options. They
remain chemically independent and may be evaluated in parallel.

## CTM and DIURN

| Workflow | What it returns | Atmosphere | Use it for |
|---|---|---|---|
| CTM | One converged snapshot per box | Embedded climatology and standard levels | Vertical profiles, climatological sensitivity tests |
| DIURN | A noon-to-noon time series plus daily summaries | Embedded or custom atmosphere | Local-time behavior, instrument grids, observation matching |

CTM is not a three-dimensional chemistry-transport model despite its legacy
name. Both workflows are collections of independent boxes.

## Three vertical selections

1. **Radiative levels** define pressure, temperature, ozone, and optical depth.
2. **Chemistry boxes** select the locations where chemistry is integrated.
3. **Exact altitudes** are DIURN chemistry boxes interpolated between explicit
   custom radiative levels.

On the embedded atmosphere, `Box.at_level(15)` selects standard level 15. On a
custom atmosphere, it selects the fifteenth supplied profile row. Exact
altitudes require a custom atmosphere with an explicit altitude coordinate.

## Inputs and outputs

Atmospheric ozone is used by radiation and can be supplied as a mixing ratio or
number density. The modeled implicit species named `o3` is an output number
density. The long-lived field also named `o3` is an input/diagnostic mixing
ratio. See {doc}`species-and-units` for the distinction.

## Equilibrium and diagnostics

An integration ending without a numerical error does not prove photochemical
equilibrium. Integration length, solver corrections, initialization, and the
stability of quantities important to the study must all be checked. See
{doc}`results` and {doc}`science-scope` before interpreting a run.
