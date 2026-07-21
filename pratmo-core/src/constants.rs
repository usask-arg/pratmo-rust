// Grid dimensions (PARAMETER statements from bcomm.h)
pub const NB: usize = 25; // max boxes
/// Storage capacity for radiative-atmosphere levels. The Fortran model uses
/// 41 levels, while the later C++ box model defaults to 81 one-kilometre
/// shells from 0 to 80 km.
pub const NL: usize = 81;
/// Atmospheric-state capacity includes the radiative shells plus one optional
/// off-grid chemistry level per box. Only the first `ModelState::nc` entries
/// participate in radiative transfer.
pub const NATM: usize = NL + NB;
pub const LEGACY_NL: usize = 41;
pub const NTAU: usize = NL + 1;
pub const NQ: usize = 77; // wavelength bins
pub const NXS: usize = 48; // cross-section species (40 + 8 iodine J-values)
pub const NSPEC: usize = 59; // total species in file (40 implicit + 19 long-lived)
pub const NSLOWM: usize = 11;
pub const NNDXPQ: usize = 25;
pub const NTAB: usize = 40; // tabulation points
pub const NJH2O: usize = 89; // H2O photolysis table points

// Derived sizes
pub const NJVAL: usize = NXS + 4; // total J-values (52)
pub const NDEN: usize = 40; // implicit species count (30 base + 10 iodine)
pub const NR: usize = 264; // reaction rate array size (250 legacy + 14 iodine extensions)
pub const NRCOLUM: usize = 430; // legacy Fortran diagnostic column layout
pub const NPMEAN: usize = 460 + NDEN; // daily mean storage (P-L block starts at 460)
pub const PPMEAN_FAMILY_OFFSET: usize = NDEN;
pub const PPMEAN_RATE_OFFSET: usize = NDEN + 20;
pub const NPPMEAN: usize = PPMEAN_RATE_OFFSET + NNDXPQ;
// The legacy Fortran grid uses at most 44 points. The default C++ port uses a
// 49-point fixed-cos(SZA) grid, so retain enough room for either integration.
pub const NXNOFT: usize = 64;
