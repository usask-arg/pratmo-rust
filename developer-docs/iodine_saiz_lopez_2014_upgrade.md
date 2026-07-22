# Iodine chemistry upgrade: Saiz-Lopez et al. (2014) — historical plan

> The implemented mechanism and its remaining limitations are documented in
> [IODINE_CHEMISTRY.md](IODINE_CHEMISTRY.md). This file is retained as planning
> history and contains superseded placeholders and proposed values.

**Reference:** Saiz-Lopez, A. et al., "Iodine chemistry in the troposphere and its effect on ozone,"
*Atmos. Chem. Phys.*, 14, 13119–13143, 2014. https://doi.org/10.5194/acp-14-13119-2014

The paper adds IO self-reactions, OIO chemistry, and higher iodine oxides (I₂O₂, I₂O₃, I₂O₄)
with fast photolysis. Its key result is that IxOy photolysis (using cross-sections from Gomez Martín
et al. 2013) rapidly recycles iodine back to I/IO rather than letting it accumulate in reservoir
species, amplifying catalytic ozone destruction.

The current pratmo scheme has 5 iodine species (I, IO, HOI, IONO₂, HI) and 12 reactions
(r[221]–r[232]). This upgrade adds 5 species, 5 photolysis channels, and 13 reactions.

---

## Data still needed from the papers

Before starting implementation, source the following. Everything else in this document is
derived from IUPAC 2008/2013 or JPL 2019 and can be implemented now.

### REQUIRED — Gomez Martín et al. 2013 (PCCP 15, 15612–15622, doi:10.1039/c3cp51217g)

The White Rose eprint is at https://eprints.whiterose.ac.uk/87171/ (accepted manuscript PDF).

Need from this paper:
- UV-vis absorption cross-sections σ(λ) for **I₂O₂**, **I₂O₃**, and **I₂O₄** as a function of
  wavelength, at one or two temperatures (250 K and/or 298 K).
- Confirm the assumed photolysis product channels:
  - I₂O₂ + hν → 2 IO  (or I + OIO?)
  - I₂O₃ + hν → IO + OIO
  - I₂O₄ + hν → 2 OIO  (or I₂O₂ + O₂?)

These go into `pratmo-core/data/fort10_cam06.x` as new cross-section tables, in the same 77-bin
format as the existing J(IO)/J(HOI)/J(IONO₂) entries (lines 1438–1550 of that file).

### OPTIONAL — Saiz-Lopez et al. 2014 Table 1 (reaction mechanism)

Verify the exact rate constants and branching ratios used for:
- IO + IO → OIO + I  (branching fraction α₁)
- IO + IO → I₂ + O₂  (branching fraction α₂)
- IO + IO + M → I₂O₂  (branching fraction α₃, and Troe parameters k₀, k∞)
- IO + OIO → I₂O₃  (k)
- OIO + OIO → I₂O₄  (k)
- OIO + OH → products  (k and product channel)

The values used below are from IUPAC 2008 (IO+IO total rate) with branching from the Saiz-Lopez
paper. If the paper's Table 1 has different values, use those instead.

---

## 1. Constants (`pratmo-core/src/constants.rs`)

| Constant | Current | New | Notes |
|----------|---------|-----|-------|
| `NDEN`   | 35      | 40  | +5 iodine species |
| `NXS`    | 43      | 48  | +5 iodine photolysis channels |

`NPMEAN = 460 + NDEN` updates automatically. `NRCOLUM` is a documentation constant not used
in the Rust struct; no change needed.

```rust
pub const NXS:  usize = 48;   // cross-section species (40 + 8 iodine J-values)
pub const NDEN: usize = 40;   // implicit species (30 base + 10 iodine)
```

---

## 2. New species

Slots n[35]–n[39] (0-based), mapped to density arrays in `state.rs` and `init.rs`.

| n index | Species | Role |
|---------|---------|------|
| n[35]   | OIO     | produced by IO+IO self-reaction; recycled by photolysis and OIO+NO |
| n[36]   | I₂      | oceanic/boundary source; photolysed to 2I in < 1 s in sunlight |
| n[37]   | I₂O₂    | termolecular IO+IO+M product; photolysis → 2IO |
| n[38]   | I₂O₃    | from IO+OIO; photolysis → IO+OIO |
| n[39]   | I₂O₄    | from OIO+OIO; photolysis → 2OIO |

### `state.rs` additions

After the existing `di_`, `dio`, `dhoi`, `diono2`, `dhi` density arrays, add:

```rust
pub doio:   [f64; NB],
pub di2:    [f64; NB],
pub di2o2:  [f64; NB],
pub di2o3:  [f64; NB],
pub di2o4:  [f64; NB],
```

After the existing `vio`, `vhoi`, `viono2` J-value arrays, add:

```rust
pub voio:   [f64; NL],   // J(OIO)
pub vi2:    [f64; NL],   // J(I₂)
pub vi2o2:  [f64; NL],   // J(I₂O₂)
pub vi2o3:  [f64; NL],   // J(I₂O₃)
pub vi2o4:  [f64; NL],   // J(I₂O₄)
```

The proposed standalone `fi2` source field was not retained. I2 is part of the
gas-phase `fiodx` family and is initialized/repartitioned by chemistry.

### `init.rs` changes

Assign the new species slots:

```rust
// OIO, I2, I2O2, I2O3, I2O4 — all initialise to zero; built up by chemistry
n[35] = pointer to doio array;
n[36] = pointer to di2 array;
n[37] = pointer to di2o2 array;
n[38] = pointer to di2o3 array;
n[39] = pointer to di2o4 array;
```

Set species names in `tname`/`tnamet`:

```rust
tname[35] = "OIO".to_string();
tname[36] = "I2".to_string();
tname[37] = "I2O2".to_string();
tname[38] = "I2O3".to_string();
tname[39] = "I2O4".to_string();
```

The `fiodx` family conservation (currently I + IO + HOI + IONO₂ + HI) does **not** need to
include the new species at initialisation — they start at zero and grow from chemistry. However,
the long-run total Iy diagnostic should eventually sum them all.

---

## 3. New photolysis cross-sections (`pratmo-core/data/fort10_cam06.x`)

Five new tables appended after the existing J(IONO₂) block (currently ending around line 1550).
Format matches the existing iodine entries: title line, then two rows of 77 wavelength-bin values
at 250 K and 298 K.

### J(OIO) — source: Ingham et al. (2000) + Spietz et al. (2005)

OIO has two absorption regions:
- UV: broad continuum 300–370 nm, σ_peak ≈ 4×10⁻¹⁸ cm²
- Visible: structured bands 480–700 nm, σ_peak ≈ 2×10⁻¹⁷ cm²

These cross-sections are tabulated in Ingham et al. (2000), *J. Phys. Chem. A* 104, 2561, and
Spietz et al. (2005), *J. Photochem. Photobiol. A* 176, 50. Both papers are open-access or
available through standard library access.

### J(I₂) — source: JPL 2019 Table 4-1 (freely available)

I₂ absorbs across 400–700 nm, peaking near 500 nm (σ_max ≈ 6×10⁻¹⁸ cm²). The cross-sections
are well-established and tabulated in the JPL Data Evaluation (Table 4-1, p. 4-82 of JPL 19-5).

### J(I₂O₂), J(I₂O₃), J(I₂O₄) — **source: Gomez Martín et al. (2013) — DATA NOT YET SOURCED**

These three tables are the only blockers. Once the cross-section values are available from the
paper, format them as:

```
J(I2O2)  250K  298K
<77 values at 250 K, space-separated, 10 per line>
<77 values at 298 K, space-separated, 10 per line>
```

If only one temperature is reported, duplicate the row (the difference is expected to be small
for these species). If the paper reports σ only at specific wavelengths rather than across all
77 bins, interpolate linearly between reported points and set σ = 0 outside the measured range.

**Placeholder approach (to unblock code development):** Use a uniform σ = 1×10⁻¹⁷ cm² across
290–400 nm and zero elsewhere as a temporary substitute. Label the block clearly:

```
J(I2O2) [PLACEHOLDER - replace with Gomez Martin 2013 Table values]
```

---

## 4. New rate constants (`pratmo-core/data/fort11_jpl09.x`)

Add after the current last iodine entry (`RATE/ 230`, IO+BrO). File uses 1-based indices, so
code r[N] corresponds to `RATE/ N+1`.

Note on IO+IO branching: the total rate is k_total = 5.4×10⁻¹¹·exp(+180/T) cm³/s (IUPAC 2008).
The three entries below share this total with branching fractions α₁ = α₂ = 0.35, α₃ = 0.30.
**Verify these branching fractions against Saiz-Lopez 2014 Table 1.**

```
RATE/  235  1.89E-11     -180.     IO+IO=OIO+I          --IUPAC08,alpha=0.35--
RATE/  236  1.89E-11     -180.     IO+IO=I2+O2          --IUPAC08,alpha=0.35--
RATE/2 237  1.62E-31       5.0     IO+IO+M=I2O2         --IUPAC08,alpha=0.30--(k0)
       237  1.00E-10       0.0     I2O2                 (kinf, approximate)
       237    -.5108        0.     I2O2 EXP=LN(0.6)
RATE/  238  6.70E-12        0.     OIO+NO=IO+NO2        --GomezM09--
RATE/  239  1.00E-10        0.     OIO+OH=HOI+O2        --estimate,uncertain--
RATE/  240  5.00E-11        0.     IO+OIO=I2O3          --SaizL14--
RATE/  241  2.00E-11        0.     OIO+OIO=I2O4         --SaizL14--
RATE/  242  2.10E-10        0.     I2+OH=HOI+I          --IUPAC13--
```

Reaction index mapping (0-based code ↔ 1-based file):

| Code r[N] | File RATE/ | Reaction |
|-----------|------------|----------|
| r[234]    | 235        | IO + IO → OIO + I |
| r[235]    | 236        | IO + IO → I₂ + O₂ |
| r[236]    | 237        | IO + IO + M → I₂O₂ (termolecular) |
| r[237]    | 238        | OIO + NO → IO + NO₂ |
| r[238]    | 239        | OIO + OH → HOI + O₂ |
| r[239]    | 240        | IO + OIO → I₂O₃ |
| r[240]    | 241        | OIO + OIO → I₂O₄ |
| r[241]    | 242        | I₂ + OH → HOI + I |

Photolysis reactions (r[243]–r[247]) are computed directly from J-values in `chems()` and
do not need fort11 entries.

---

## 5. `chemistry.rs` changes

### `setupr()` — rate constant setup

Inside the `if s.liod {` block, after the existing `s.ratek[229]` line:

```rust
s.ratek[234] = fk(tt, 234);                     // IO + IO → OIO + I
s.ratek[235] = fk(tt, 235);                     // IO + IO → I₂ + O₂
s.ratek[236] = fknasa(tt300, zdnum, 236, s);    // IO + IO + M → I₂O₂
s.ratek[237] = fk(tt, 237);                     // OIO + NO → IO + NO₂
s.ratek[238] = fk(tt, 238);                     // OIO + OH → HOI + O₂
s.ratek[239] = fk(tt, 239);                     // IO + OIO → I₂O₃
s.ratek[240] = fk(tt, 240);                     // OIO + OIO → I₂O₄
s.ratek[241] = fk(tt, 241);                     // I₂ + OH → HOI + I
```

### `chems()` — species aliases and J-values

Add after the existing iodine species aliases (around line 361):

```rust
let xoio  = if s.liod && n[35] > 0 { xr(s, n[35]) } else { 0.0 };
let xi2   = if s.liod && n[36] > 0 { xr(s, n[36]) } else { 0.0 };
let xi2o2 = if s.liod && n[37] > 0 { xr(s, n[37]) } else { 0.0 };
let xi2o3 = if s.liod && n[38] > 0 { xr(s, n[38]) } else { 0.0 };
let xi2o4 = if s.liod && n[39] > 0 { xr(s, n[39]) } else { 0.0 };
```

Add after the existing J-value aliases (around line 422):

```rust
let jvoio  = s.voio[iv];
let jvi2   = s.vi2[iv];
let jvi2o2 = s.vi2o2[iv];
let jvi2o3 = s.vi2o3[iv];
let jvi2o4 = s.vi2o4[iv];
```

### `chems()` — reaction rates

Add after `r[233]` inside `if s.liod {`:

```rust
// New kinetic reactions
r[234] = s.ratek[234] * xio   * xio;    // IO + IO → OIO + I
r[235] = s.ratek[235] * xio   * xio;    // IO + IO → I₂ + O₂
r[236] = s.ratek[236] * xio   * xio;    // IO + IO + M → I₂O₂
r[237] = s.ratek[237] * xoio  * xno;    // OIO + NO → IO + NO₂
r[238] = s.ratek[238] * xoio  * xoh;    // OIO + OH → HOI + O₂
r[239] = s.ratek[239] * xio   * xoio;   // IO + OIO → I₂O₃
r[240] = s.ratek[240] * xoio  * xoio;   // OIO + OIO → I₂O₄
r[241] = s.ratek[241] * xi2   * xoh;    // I₂ + OH → HOI + I

// New photolysis reactions
r[243] = jvoio  * xoio;                 // OIO + hν → I + O₂
r[244] = jvi2   * xi2;                  // I₂ + hν → 2I
r[245] = jvi2o2 * xi2o2;               // I₂O₂ + hν → 2IO
r[246] = jvi2o3 * xi2o3;               // I₂O₃ + hν → IO + OIO
r[247] = jvi2o4 * xi2o4;               // I₂O₄ + hν → 2OIO
```

Note: r[234] and r[235] are both `xio * xio` (second-order in IO). Each reaction event
consumes **two** IO molecules, which is handled in the loss terms below.

### `chempl()` — production/loss updates

**Patch existing species** (add new terms to existing P/L sums):

```rust
// I — additional production
rp[idx(n[30])] += r[234]          // IO + IO → OIO + I
                + r[241]          // I₂ + OH → HOI + I  (yields 1 I)
                + 2.0*r[244]      // I₂ + hν → 2I
                + r[243];         // OIO + hν → I + O₂

// IO — additional loss (note ×2: each IO+IO event consumes 2 IO)
rl[idx(n[31])] += 2.0*r[234]
                + 2.0*r[235]
                + 2.0*r[236]
                + r[239];         // IO + OIO → I₂O₃

// IO — additional production
rp[idx(n[31])] += r[237]          // OIO + NO → IO + NO₂
                + 2.0*r[245]      // I₂O₂ + hν → 2IO
                + r[246];         // I₂O₃ + hν → IO + OIO

// HOI — additional production
rp[idx(n[32])] += r[238]          // OIO + OH → HOI + O₂
                + r[241];         // I₂ + OH → HOI + I

// OH — additional loss
rl[idx(n[6])]  += r[238] + r[241];

// NO — additional loss
rl[idx(n[0])]  += r[237];

// NO₂ — additional production
rp[idx(n[1])]  += r[237];
```

**New species P/L blocks** (add after the existing HI block):

```rust
// OIO
if n[35] > 0 {
    rp[idx(n[35])] = r[234] + r[246] + 2.0*r[247];
    rl[idx(n[35])] = r[237] + r[238] + r[239] + 2.0*r[240] + r[243];
}

// I₂
if n[36] > 0 {
    rp[idx(n[36])] = r[235];
    rl[idx(n[36])] = r[241] + r[244];
}

// I₂O₂
if n[37] > 0 {
    rp[idx(n[37])] = r[236];
    rl[idx(n[37])] = r[245];
}

// I₂O₃
if n[38] > 0 {
    rp[idx(n[38])] = r[239];
    rl[idx(n[38])] = r[246];
}

// I₂O₄
if n[39] > 0 {
    rp[idx(n[39])] = r[240];
    rl[idx(n[39])] = r[247];
}
```

---

## 6. `reader.rs` — J-value label lookup

The fort10 reader recognises iodine J-values by title string (around line 369). Extend the
key lookup to map the five new labels to slots 43–47:

```rust
"J(OIO)"  => NTAB + 3,   // slot 43
"J(I2)"   => NTAB + 4,   // slot 44
"J(I2O2)" => NTAB + 5,   // slot 45
"J(I2O3)" => NTAB + 6,   // slot 46
"J(I2O4)" => NTAB + 7,   // slot 47
```

---

## 7. `tests/iodine_chemistry.rs` — test updates

- **`test_iodine_conservation`**: include OIO, I₂, I₂O₂, I₂O₃, I₂O₄ in the Iy sum
  (the 30% tolerance may need widening during spin-up; check after integration)
- **Add `test_oio_daytime`**: daytime OIO should be non-zero but smaller than IO
  (OIO is a transient intermediate, not a large reservoir)
- **Add `test_i2_photolysed`**: daytime I₂ should be near zero (photolysed in < 1 s)
- **Add `test_ixoy_do_not_accumulate`**: I₂O₂, I₂O₃, I₂O₄ should each be ≪ total Iy
  (the whole point of the fast-photolysis result in Saiz-Lopez 2014)
- **Update non-regression**: O₃ depletion rate should be equal or greater than before
  (new scheme adds ozone loss pathways; a regression would be a bug)

---

## 8. Implementation order

1. `constants.rs` — bump NDEN and NXS first; everything else follows from this
2. `state.rs` — add arrays (compiler will flag all missing initialisations)
3. `init.rs` — initialise new species slots
4. `fort11_jpl09.x` — add kinetic rate entries (no code change needed to read them)
5. `reader.rs` — extend J-value key lookup
6. `fort10_cam06.x` — add J(OIO) and J(I₂) now; add IxOy blocks as stubs until
   Gomez Martín data is sourced
7. `chemistry.rs` — setupr, chems, chempl changes
8. Tests

The model will compile and run with placeholder IxOy cross-sections; the key science result
(fast IxOy recycling) will only be correct once the real Gomez Martín 2013 values are entered.
