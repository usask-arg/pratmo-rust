// Grid dimensions (PARAMETER statements from bcomm.h)
pub const NB: usize = 25;       // max boxes
pub const NL: usize = 41;       // altitude levels
pub const NQ: usize = 77;       // wavelength bins
pub const NXS: usize = 48;      // cross-section species (40 + 8 iodine J-values)
pub const NSPEC: usize = 59;    // total species in file (40 implicit + 19 long-lived)
pub const NSLOWM: usize = 11;
pub const NNDXPQ: usize = 25;
pub const NTAB: usize = 40;     // tabulation points
pub const NJH2O: usize = 89;    // H2O photolysis table points

// Derived sizes
pub const NJVAL: usize = NXS + 4;   // total J-values (52)
pub const NDEN: usize = 40;         // implicit species count (30 base + 10 iodine)
pub const NR: usize = 250;          // reaction rate array size
pub const NRCOLUM: usize = 430;     // P-L column storage (XR + RP + RL + RPF + RLF + RQF)
pub const NPMEAN: usize = 460 + NDEN; // daily mean storage (P-L block starts at 460)
pub const NXNOFT: usize = 44;       // time steps per day (max)
