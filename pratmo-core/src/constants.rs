// Grid dimensions (PARAMETER statements from bcomm.h)
pub const NB: usize = 25;       // max boxes
pub const NL: usize = 41;       // altitude levels
pub const NQ: usize = 77;       // wavelength bins
pub const NXS: usize = 40;      // cross-section species
pub const NSPEC: usize = 49;    // total species in file
pub const NSLOWM: usize = 11;
pub const NNDXPQ: usize = 25;
pub const NTAB: usize = 40;     // tabulation points
pub const NJH2O: usize = 89;    // H2O photolysis table points

// Derived sizes
pub const NJVAL: usize = NXS + 4;   // total J-values (44)
pub const NDEN: usize = 30;         // implicit species count
pub const NR: usize = 250;          // reaction rate array size
pub const NRCOLUM: usize = 430;     // P-L column storage (XR + RP + RL + RPF + RLF + RQF)
pub const NPMEAN: usize = 490;      // daily mean storage
pub const NXNOFT: usize = 44;       // time steps per day (max)
