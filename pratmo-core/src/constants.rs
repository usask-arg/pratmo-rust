// Grid dimensions (PARAMETER statements from bcomm.h)
pub const NB: usize = 25;       // max boxes
pub const NL: usize = 41;       // altitude levels
pub const NQ: usize = 77;       // wavelength bins
pub const NXS: usize = 43;      // cross-section species (40 + 3 iodine: J_IO, J_HOI, J_IONO2)
pub const NSPEC: usize = 54;    // total species in file (35 implicit + 19 long-lived)
pub const NSLOWM: usize = 11;
pub const NNDXPQ: usize = 25;
pub const NTAB: usize = 40;     // tabulation points
pub const NJH2O: usize = 89;    // H2O photolysis table points

// Derived sizes
pub const NJVAL: usize = NXS + 4;   // total J-values (47)
pub const NDEN: usize = 35;         // implicit species count (30 base + 5 iodine: I, IO, HOI, IONO2, HI)
pub const NR: usize = 250;          // reaction rate array size
pub const NRCOLUM: usize = 430;     // P-L column storage (XR + RP + RL + RPF + RLF + RQF)
pub const NPMEAN: usize = 460 + NDEN; // daily mean storage (P-L block starts at 460)
pub const NXNOFT: usize = 44;       // time steps per day (max)
