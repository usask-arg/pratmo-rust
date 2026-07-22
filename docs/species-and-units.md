# Species, families, and units

PRATMO exposes three categories with different units.

| Category | Accessor | Unit | Meaning |
|---|---|---|---|
| Implicit species | `species_profile`, `species_grid` | molecules cm⁻³ | Short-lived chemical state |
| Long-lived species/families | `long_lived_profile` | Dimensionless volume mixing ratio | Prescribed or slowly varying family state |
| J-values | `jvalue_profile`, `jvalue_grid` | s⁻¹ | Photolysis frequencies |

Use `mixing_ratio_as(values, "ppbv")` or
`number_density_as(values, "m-3")` for display conversions. Arrays returned by
the model remain in canonical units.

## Implicit species

| Group | Python field names |
|---|---|
| Ox and hydrogen | `h`, `oh`, `ho2`, `h2o2`, `o`, `o3` |
| Reactive nitrogen | `no`, `no2`, `no3`, `n2o5`, `hno2`, `hno3`, `hno4` |
| Chlorine | `cl`, `cl2`, `clo`, `hcl`, `hocl`, `clono2`, `oclo`, `cl2o2` |
| Bromine | `br`, `bro`, `hbr`, `hobr`, `brono2`, `brcl` |
| Carbon/peroxy | `h2co`, `ch3o2`, `ch3o2h` |
| Experimental iodine | `i`, `io`, `hoi`, `iono2`, `hi`, `oio`, `i2`, `i2o2`, `i2o3`, `i2o4` |

## Long-lived fields and families

| Field | Chemical quantity |
|---|---|
| `o3` | Ozone mixing ratio used for initialization/reference |
| `n2o` | Nitrous oxide |
| `noy` | Total reactive nitrogen family |
| `ch4` | Methane |
| `co` | Carbon monoxide |
| `clx` | Total inorganic/active chlorine family used by PRATMO |
| `cf2cl2` | CFC-12 |
| `cfcl3` | CFC-11 |
| `ccl4` | Carbon tetrachloride |
| `ch3cl` | Methyl chloride |
| `ch3ccl3` | Methyl chloroform |
| `h2` | Molecular hydrogen |
| `h2o` | Water vapour |
| `nh3` | Ammonia |
| `c5h8` | Isoprene |
| `brx` | Total bromine family used by PRATMO |
| `ch3br` | Methyl bromide |
| `ocs` | Carbonyl sulfide |
| `iodx` | Total iodine family used by the experimental extension |

The name `o3` therefore appears in both categories. Use
`species_profile("o3")` for modeled O3 number density and
`long_lived_profile("o3")` for the O3 mixing-ratio field.

## Discover names programmatically

```python
from pratmo import IMPLICIT_SPECIES_NAMES, JVALUE_NAMES, LONG_LIVED_NAMES

print(IMPLICIT_SPECIES_NAMES)
print(LONG_LIVED_NAMES)
print(JVALUE_NAMES)
```

Names are case-insensitive when passed to output accessors. The constants are
the authoritative complete lists for the installed version.
