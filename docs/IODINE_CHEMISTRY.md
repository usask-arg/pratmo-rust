# Iodine chemistry

```{warning}
This mechanism is experimental and has not been scientifically validated. The
tests described below check implementation properties such as stoichiometry,
conservation, and Jacobian consistency; they do not establish atmospheric
accuracy or fitness for research use.
```

## Scope and sources

PRATMO carries ten gas-phase inorganic iodine species: `I`, `IO`, `HOI`,
`IONO2`, `HI`, `OIO`, `I2`, `I2O2`, `I2O3`, and `I2O4`.

The reaction mechanism follows Tables 1–3 of Saiz-Lopez et al. (2014), with
higher-oxide absorption spectra from Lewis et al. (2020). The mechanism is a
stratospheric box-model subset: species that require expanding the state space
(`INO`, `INO2`, `ICl`, `IBr`, `IO3`, `HIO3`, and particulate iodate) are not yet
carried explicitly.

Primary references:

- Saiz-Lopez et al. (2014), https://doi.org/10.5194/acp-14-13119-2014
- Lewis et al. (2020), https://doi.org/10.5194/acp-20-10865-2020
- Cuevas et al. (2022), https://doi.org/10.1073/pnas.2110864119

## Numerical implementation

Iodine reactions have one stoichiometric definition in `chemistry.rs`.
Production/loss assembly and the analytic Newton Jacobian both consume that
definition. Tests check iodine-atom conservation for every listed reaction and
finite-difference the IONO2 photolysis Jacobian column.

The pressure-dependent IO self-reaction, IO + OIO, OIO + OIO, and thermal
I2O2/I2O4 decomposition expressions are evaluated directly from the
Saiz-Lopez parameterizations using pressure in hPa.

The prescribed gas-phase family is

```
Iy = I + IO + HOI + IONO2 + HI + OIO
   + 2 (I2 + I2O2 + I2O3 + I2O4)
```

Because particulate iodine is not represented, irreversible IxOy aerosol loss
is not applied. This avoids removing gas iodine in CHEMPL and then silently
restoring it through family normalization.

## Photolysis

Eight iodine J-values are required whenever iodine species are active. The
reader rejects partial iodine spectral extensions instead of silently leaving
missing tables at zero.

- OIO and I2Oy spectra are labelled with their Lewis (2020) provenance.
- The current IO, HOI, and IONO2 spectra remain legacy estimates and are marked
  `LEGACY-EST` in `fort10_cam06.x`.
- The calculation assumes unit photolysis quantum yield. This is a documented
  uncertainty, especially for higher-oxide product branching.

The implemented higher-oxide products follow the Saiz-Lopez mechanism:

- `I2O2 + hv -> I + OIO`
- `I2O3 + hv -> IO + OIO`
- `I2O4 + hv -> 2 OIO`

Lewis et al. (2020) shows that IO3-containing channels may also occur. Adding
those channels requires adding IO3/HIO3 to the model state.

## Heterogeneous recycling

Generic aerosol and sea-salt surface areas are separate inputs:

- `aerosol_surface_area_um2_cm3` drives legacy sulfate/general aerosol
  reactions.
- `sea_salt_surface_area_um2_cm3` drives iodine sea-salt recycling only.

Saiz-Lopez represents HOI and IONO2 uptake as equal ICl and IBr production.
Until ICl and IBr are explicit species, PRATMO uses their rapid-photolysis
limit: iodine returns as I and the corresponding Cl/Br radicals are returned.
The ICl branch of IO + ClO is similarly combined with the atomic branch as a
documented proxy.

## Iodine-off controls and convergence

For standard-grid and CTM sensitivity runs, `iodine = false` retains the
40-row Newton structure and sets total iodine to a `1e-20` mixing-ratio
numerical floor. This is eight orders of magnitude below the standard 1 ppt
case and avoids changing nonlinear-system dimension between paired runs.
Custom-atmosphere runs retain their established 30-species iodine-off path so
box-slot invariance is preserved.

`Diagnostics` exposes `newraf_nonconvergence_count` and
`rafday_nonconvergence_count`. Scientific sensitivity runs should require both
to be zero.

Normal-mode RAFDAY applies its multi-day relaxation once, before the
daily-mean Newton corrections. Repeating that relaxation after every correction
integrated HNO3 forward again and erased its update; `FIXRAT` then propagated
the error into BrONO2 through the coupled NOy/Bry constraint. The one-time
relaxation keeps the fast diurnal orbit initialized without resetting the slow
family solve.

## Known boundaries

- IO, HOI, and IONO2 cross sections still need fully traced evaluated tables
  and explicit quantum yields.
- ICl/IBr reservoirs and their finite photolysis lifetimes are approximated.
- IO3/HIO3 and gas-particle iodate cycling are not represented.
- The original Fortran model has no iodine mechanism. Mechanism-level tests and
  literature comparisons are available, but the iodine extension has not been
  independently or scientifically validated.
