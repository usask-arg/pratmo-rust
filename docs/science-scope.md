# Scientific scope and assumptions

PRATMO is intended for stratospheric photochemistry experiments in independent
boxes. It is not a complete atmosphere model.

## Included

- A full solar-day photochemical calculation.
- Gas-phase chemistry for the documented implicit and long-lived fields.
- Bromine and chlorine chemistry.
- Sulfate-aerosol and sea-salt heterogeneous reactions when enabled.
- Radiative transfer and photolysis based on an embedded or caller-provided
  vertical atmosphere.
- An opt-in experimental inorganic iodine extension.

## Not represented as prognostic processes

- Horizontal or vertical transport between boxes.
- Meteorological dynamics or feedback on pressure and temperature.
- Interactive emissions, deposition, or aerosol microphysics.
- A year-specific reanalysis or emissions inventory.
- A scientifically validated tropospheric or operational forecast system.

Long-lived inputs, atmosphere profiles, and aerosol surface areas are boundary
conditions. A multi-box result is therefore a profile of independent
calculations, not a coupled atmospheric column.

## Time and location

Latitude and day of year control solar geometry and climatological sampling.
The calendar year in an ISO date is not otherwise used. Local time in DIURN is
local solar time, not a civil time zone.

## Atmospheric domain

The embedded data and documented defaults focus on the stratosphere. The
standard grid extends outside it for legacy compatibility, but availability of
a level does not establish scientific validity there. For an instrument grid,
provide measured pressure, temperature, altitude, and ozone and document the
source of every substituted field.

## What warnings mean

`PratmoWarning` identifies likely unit mistakes, unusually short integrations,
or suspicious but technically possible settings. It is not a scientific
quality-control system. A run without warnings can still be inappropriate for
the intended study.
