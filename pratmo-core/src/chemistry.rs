// bchem.f → chemistry module
// SETUPR, CHEMS, CHEMPL, RHSLHS, REACT, FKNASA

use nalgebra::{DMatrix, DVector};

use crate::state::ModelState;

// ── Index helpers (Fortran 1-based N1..N30 stored as s.n[0]..s.n[29]) ────────

/// Read XR for species Ni (where ni = s.n[i-1], 1-based Fortran value).
#[inline(always)]
fn xr(s: &ModelState, ni: usize) -> f64 {
    s.xr[ni - 1]
}

/// 0-based index from Fortran 1-based species index.
#[inline(always)]
fn idx(ni: usize) -> usize {
    ni - 1
}

// ── SETUPR — temperature-dependent rate constants ─────────────────────────────

/// Compute all RATEK(1..221) for current T, M, altitude.
/// Fortran: SUBROUTINE SETUPR
pub fn setupr(s: &mut ModelState) {
    let ialt  = s.ialt;   // 0-based
    let ibox  = s.ibox;   // 0-based
    let zalt  = s.z[ialt];

    // Temperature with optional box offset
    let tt = s.t[ialt] + s.boxtt[ibox];
    let tt300 = (tt / 300.0_f64).ln();
    s.ztemp = tt;
    s.zdnum = s.dm[ialt];
    s.zalt  = zalt;
    s.izalt = ialt;

    let zdnum = s.zdnum;
    let zdh2o = s.fh2o[ibox] * zdnum;

    // Aerosol surface area: use climatology ASA if available, else BOXAA
    if s.asa[ialt] >= 0.0 {
        s.boxaa[ibox] = s.asa[ialt];
    }
    let zaersl = 675.0 * tt.sqrt() * s.boxaa[ibox] * 1.0e-8;
    let zraino = s.boxrn[ibox] / 86400.0;

    // Arrhenius helpers
    let fk  = |temp: f64, i: usize| s.rk[[0, i]] * (-s.rk[[1, i]] / temp).exp();
    let fkx = |tl:   f64, i: usize| s.rk[[0, i]] * (-tl * s.rk[[1, i]]).exp();

    // ── Rate constants (0-based arrays; Fortran RATEK(J) → s.ratek[j-1]) ────

    // R4: O(1D) + M (two channels, N2 and O2)
    let ndx4 = s.ndxrat[3] as usize; // NDXRAT(4) 0-based
    s.ratek[3] = (s.pn2 * fk(tt, 3) + s.po2 * s.rkadd[[0, ndx4]] * (-s.rkadd[[1, ndx4]] / tt).exp())
                 * zdnum;
    s.ratek[4]  = fk(tt, 4);   // O(1D)+H2O
    s.ratek[5]  = fk(tt, 5);   // O(1D)+H2
    s.ratek[6]  = fk(tt, 6);   // O(1D)+CH4 → OH+CH3
    s.ratek[7]  = fk(tt, 7);   // O(1D)+CH4 → H2+CH2O
    s.ratek[8]  = fk(tt, 8);   // O(1D)+N2O → NO+NO
    s.ratek[9]  = fk(tt, 9);   // O(1D)+N2O → N2+O2
    s.ratek[10] = fkx(tt300, 10) * zdnum; // O+O2+M
    s.ratek[11] = fk(tt, 11);  // O+O3
    s.ratek[12] = fk(tt, 12) * zdnum; // O+O+M
    s.ratek[13] = fk(tt, 13);  // O+H2
    s.ratek[15] = fk(tt, 15);  // O+OH
    s.ratek[16] = fk(tt, 16);  // O+HO2
    s.ratek[17] = fk(tt, 17);  // O+H2O2
    s.ratek[18] = fk(tt, 18);  // O+NO2
    s.ratek[20] = fk(tt, 20);  // O3+NO
    s.ratek[21] = fk(tt, 21);  // O(1D)+O3
    let ndx23 = s.ndxrat[22] as usize;
    s.ratek[22] = fkx(tt300, 22) * zdnum; // O+NO+M (3-body)
    let _ = ndx23;
    s.ratek[23] = fkx(tt300, 23) * zdnum; // O+NO2+M
    s.ratek[24] = fk(tt, 24);  // H+O3
    s.ratek[25] = fk(tt, 25);  // OH+O3
    s.ratek[26] = fk(tt, 26);  // O3+HO2
    s.ratek[27] = fk(tt, 27);  // O3+NO2
    // OCS reactions
    s.ratek[28] = fk(tt, 28);  // O+OCS
    s.ratek[29] = fk(tt, 29);  // OH+OCS
    s.ratek[31] = fkx(tt300, 31) * zdnum; // H+O2+M
    s.ratek[32] = fk(tt, 32);  // H+HO2 → OH+OH
    s.ratek[33] = fk(tt, 33);  // H+HO2 → H2+O2
    s.ratek[34] = fkx(tt300, 34) * s.pn2 * zdnum * zdnum; // H+HO2+N2 (3rd order)
    s.ratek[35] = fk(tt, 35);  // O+NO3
    // Source terms (ppt/day) → #/cm³/s
    s.ratek[37] = s.rk[[0, 37]] * zdnum / 86400.0e12; // H2O2 source
    s.ratek[38] = fk(tt, 38);  // OH+OH
    let ndx40 = s.ndxrat[39] as usize;
    s.ratek[39] = fk(tt, 39) + zdnum * s.rkadd[[0, ndx40]] * (-s.rkadd[[1, ndx40]] / tt).exp();
    s.ratek[40] = fk(tt, 40);  // OH+H2O2
    s.ratek[41] = fk(tt, 41);  // HO2+NO
    s.ratek[42] = fk(tt, 42) * zaersl; // HO2 aerosol loss
    let ndx44 = s.ndxrat[43] as usize;
    s.ratek[43] = (fk(tt, 43) + zdnum * s.rkadd[[0, ndx44]] * (-s.rkadd[[1, ndx44]] / tt).exp())
                  * (1.0 + zdh2o * s.rkadd[[0, ndx44 + 1]]); // HO2+HO2
    s.ratek[44] = zraino / s.rk[[0, 44]]; // H2O2 rainout
    s.ratek[45] = fk(tt, 45) * zaersl;   // H2O2 aerosol loss
    s.ratek[46] = fk(tt, 46);  // OH+H2
    let ndx48 = s.ndxrat[47] as usize;
    s.ratek[47] = fk(tt, 47) * (1.0 + s.rkadd[[0, ndx48]] * zdnum); // OH+CO
    s.ratek[48] = fk(tt, 48);  // OH+CH4
    s.ratek[50] = fk(tt, 50);  // CH3O2+NO
    s.ratek[51] = fk(tt, 51);  // CH3O2+HO2
    s.ratek[52] = fk(tt, 52);  // CH3O2+CH3O2
    s.ratek[53] = s.rk[[0, 53]]; // ROOH photolysis (constant)
    s.ratek[54] = fk(tt, 54);  // OH+ROOH
    s.ratek[55] = zraino / s.rk[[0, 55]]; // ROOH rainout
    s.ratek[56] = fk(tt, 56) * zaersl; // ROOH aerosol loss
    s.ratek[57] = fk(tt, 57);  // OH+H2CO
    s.ratek[60] = zraino / s.rk[[0, 60]]; // H2CO rainout
    s.ratek[61] = fknasa(tt300, zdnum, 61, s); // OH+NO2 → HNO3
    let ndx64 = s.ndxrat[63] as usize;
    let ek2 = s.rkadd[[0, ndx64]] * (-s.rkadd[[1, ndx64]] / tt).exp();
    let ek3 = s.rkadd[[0, ndx64 + 1]] * zdnum * (-s.rkadd[[1, ndx64 + 1]] / tt).exp();
    s.ratek[63] = fk(tt, 63) + ek2 * ek3 / (ek2 + ek3); // OH+HNO3
    s.ratek[64] = zraino / s.rk[[0, 64]]; // HNO3 rainout
    s.ratek[66] = fknasa(tt300, zdnum, 66, s); // OH+NO → HNO2
    s.ratek[67] = fk(tt, 67);  // HO2+NO2
    s.ratek[69] = fk(tt, 69);  // OH+HNO2
    s.ratek[70] = zraino / s.rk[[0, 70]]; // HNO2 rainout
    s.ratek[71] = fknasa(tt300, zdnum, 71, s); // HO2+NO2 → HNO4
    s.ratek[72] = s.ratek[71] * fk(tt, 72); // HNO4 thermal decomp
    s.ratek[74] = zraino / s.rk[[0, 74]]; // HNO4 rainout
    s.ratek[75] = fk(tt, 75);  // OH+HNO4
    s.ratek[78] = fk(tt, 78);  // NO+NO3
    s.ratek[79] = fk(tt, 79);  // NO2+NO3
    s.ratek[80] = s.rk[[0, 80]] * zdnum / 86400.0e12; // H2CO source
    s.ratek[81] = fknasa(tt300, zdnum, 81, s); // NO2+NO3 → N2O5
    s.ratek[82] = s.ratek[81] * fk(tt, 82); // N2O5 thermal decomp
    s.ratek[84] = s.rk[[0, 84]] * zaersl; // NOx aerosol loss
    s.ratek[86] = fk(tt, 86);  // N(4S)+O2
    s.ratek[87] = fk(tt, 87);  // N(4S)+O3
    s.ratek[88] = fk(tt, 88);  // N(4S)+NO
    s.ratek[89] = fk(tt, 89);  // N(4S)+NO2
    s.ratek[93] = s.rk[[0, 93]] * zdnum / 86400.0e12; // XXXX source
    s.ratek[94] = fk(tt, 94);  // OH+XXXX
    s.ratek[95] = fk(tt, 95);  // XXXX aerosol loss
    s.ratek[96] = fk(tt, 96);  // Cl+XXXX
    s.ratek[99] = fk(tt, 99);  // OH+HCl
    s.ratek[100] = fk(tt, 100); // O+HCl
    s.ratek[101] = fk(tt, 101); // HCl+ClONO2 (gas)
    s.ratek[102] = fk(tt, 102); // Cl+HO2
    s.ratek[103] = fk(tt, 103); // Cl+CH4
    s.ratek[104] = fk(tt, 104); // Cl+HNO4
    s.ratek[105] = fk(tt, 105); // Cl+H2CO
    s.ratek[106] = fk(tt, 106); // Cl+H2
    s.ratek[107] = fk(tt, 107); // Cl+HO2
    s.ratek[108] = fk(tt, 108); // Cl+H2O2
    s.ratek[109] = fk(tt, 109); // O(1D)+CCl4
    s.ratek[110] = fk(tt, 110); // O(1D)+CFCl3
    s.ratek[111] = fk(tt, 111); // O(1D)+CF2Cl2
    s.ratek[112] = fk(tt, 112); // O(1D)+HCl
    s.ratek[113] = fk(tt, 113); // O3+Cl
    s.ratek[114] = fk(tt, 114); // O+ClO
    s.ratek[115] = fk(tt, 115); // ClO+NO
    s.ratek[117] = fknasa(tt300, zdnum, 117, s); // ClO+ClO+M → Cl2O2
    s.ratek[118] = s.ratek[117] * fk(tt, 118); // Cl2O2 thermal decomp
    s.ratek[119] = fk(tt, 119); // OH+ClO
    s.ratek[120] = fk(tt, 120); // HO2+ClO
    s.ratek[121] = fknasa(tt300, zdnum, 121, s); // ClO+NO2+M → ClONO2
    s.ratek[122] = s.ratek[121] * fk(tt, 122); // ClONO2 thermal decomp
    s.ratek[123] = s.rk[[0, 123]]; // ClONO2+hv (photolysis rate stored as const)
    s.ratek[124] = s.rk[[0, 124]]; // ClONO2+hv (second channel)
    s.ratek[125] = fk(tt, 125); // ClONO2+O
    s.ratek[126] = fk(tt, 126); // ClONO2+OH
    s.ratek[127] = fk(tt, 127); // HO2+ClO (HOCl channel)
    s.ratek[129] = fk(tt, 129); // OH+HOCl
    s.ratek[142] = fk(tt, 142); // OH+CH3Cl
    s.ratek[144] = fk(tt, 144); // OH+CH3CCl3
    s.ratek[147] = fk(tt, 147); // OH+CH3Br
    s.ratek[148] = fk(tt, 148); // BrONO2+O
    s.ratek[149] = fk(tt, 149); // BrO+OH
    s.ratek[156] = fk(tt, 156); // OH+ROOH (second)
    s.ratek[158] = fk(tt, 158); // Cl+Cl+M
    s.ratek[160] = fk(tt, 160); // ClO+ClO (OClO channel)
    s.ratek[162] = fk(tt, 162); // OClO+NO
    s.ratek[163] = fk(tt, 163); // OClO+OH
    s.ratek[164] = fk(tt, 164); // OClO+Cl
    s.ratek[165] = fk(tt, 165); // Cl+ClONO2 (gas)
    s.ratek[166] = fk(tt, 166); // Cl+NO3
    s.ratek[167] = s.rk[[0, 167]] * zdnum / 86400.0e9; // O3 source (ppb/day)
    // Heterogeneous rates — computed via hetprob (called from chems/setupr)
    // apply_het_rates() sets ratek[169..176]
    s.ratek[177] = s.rk[[0, 177]] * zdnum / 86400.0e12; // NOx rain source
    s.ratek[178] = s.rk[[0, 178]] * zdnum / 86400.0e12; // NOx rain source
    s.ratek[180] = fk(tt, 180);
    s.ratek[182] = fk(tt, 182);
    s.ratek[184] = fk(tt, 184);
    s.ratek[189] = fk(tt, 189);
    s.ratek[192] = fk(tt, 192);
    s.ratek[193] = fk(tt, 193);
    s.ratek[194] = fk(tt, 194);
    s.ratek[195] = fk(tt, 195);
    s.ratek[196] = fk(tt, 196);
    s.ratek[197] = fk(tt, 197);
    s.ratek[198] = fk(tt, 198);
    s.ratek[199] = fk(tt, 199);

    if s.nrates >= 201 {
        s.ratek[200] = fk(tt, 200); // HBr+OH
        s.ratek[201] = fk(tt, 201); // HBr+O
        s.ratek[203] = fk(tt, 203); // Br+HO2
        s.ratek[204] = fk(tt, 204); // Br+O3
        s.ratek[205] = fk(tt, 205); // BrO+O
        s.ratek[206] = fk(tt, 206); // BrO+NO
        s.ratek[207] = fk(tt, 207); // BrO+O3
        let ndx209 = s.ndxrat[208] as usize;
        s.ratek[208] = fk(tt, 208)
            + s.rkadd[[0, ndx209]] * (-s.rkadd[[1, ndx209]] / tt).exp(); // BrO+BrO
        s.ratek[210] = fk(tt, 210); // BrO+HO2
        s.ratek[212] = fk(tt, 212); // OH+HOBr
        s.ratek[213] = fknasa(tt300, zdnum, 213, s); // BrO+NO2+M → BrONO2
        s.ratek[214] = s.rk[[0, 214]]; // BrONO2 photolysis k
        s.ratek[215] = s.rk[[0, 215]]; // BrONO2 photolysis k2
        s.ratek[216] = fk(tt, 216);    // ClO+BrO
        s.ratek[217] = fk(tt, 217);    // Br+H2CO
        s.ratek[218] = fk(tt, 218);    // ClO+BrO (second channel)
        s.ratek[219] = s.rk[[0, 219]]; // VBrCl/VCl2 scale
        s.ratek[220] = fk(tt, 220);    // ClO+BrO (third channel)
    }

    // Heterogeneous reactions 170–177 via hetprob
    apply_het_rates_setupr(s, zaersl);
}

/// Compute heterogeneous rate constants using hetprob (called from setupr).
fn apply_het_rates_setupr(s: &mut ModelState, zaersl: f64) {
    use crate::heterogeneous::hetprob;
    let ialt  = s.ialt;
    let ibox  = s.ibox;
    let tt    = s.ztemp;
    let zdnum = s.zdnum;
    let ph2o    = s.fh2o[ibox] * s.pstd[ialt];
    let phcl    = if zdnum > 0.0 { s.dhcl[ibox] / zdnum * s.pstd[ialt] / 1013.0 } else { 0.0 };
    let pclono2 = if zdnum > 0.0 { s.dclno3[ibox] / zdnum * s.pstd[ialt] / 1013.0 } else { 0.0 };
    let phbr    = if zdnum > 0.0 { s.dhbr[ibox] / zdnum * s.pstd[ialt] / 1013.0 } else { 0.0 };

    let g = hetprob(s.pstd[ialt], tt, 1.0e-5, ph2o, phcl, pclono2, phbr, false);
    // g[0]=γ_ClONO2+HCl(172), g[1]=γ_ClONO2+H2O(171), g[2]=γ_HOCl+HCl(174),
    // g[3]=γ_N2O5+HCl(173),   g[4]=γ_N2O5+H2O(170),   g[5]=γ_BrONO2+H2O(177),
    // g[6]=γ_HOBr+HCl(175),   g[7]=γ_HOBr+HBr(176)
    let g1 = g[0].min(1.0); // ClONO2+HCl
    let g2 = g[1];          // ClONO2+H2O
    let g3 = g[2];          // HOCl+HCl
    let g4 = g[3];          // N2O5+HCl (→ 0 per code)
    let g5 = g[4];          // N2O5+H2O
    let g6 = g[5];          // BrONO2+H2O
    let g7 = g[6];          // HOBr+HCl
    let _g8 = 0.0_f64;      // HOBr+HBr (set to 0 in Fortran)

    s.ratek[169] = g5 * zaersl * 0.5;  // N2O5+H2O (r170)
    s.ratek[170] = g2 * zaersl * 0.5;  // ClONO2+H2O (r171)
    s.ratek[171] = g1 * zaersl * 0.5;  // ClONO2+HCl (r172)
    s.ratek[172] = g4 * zaersl * 0.5;  // N2O5+HCl (r173)
    s.ratek[173] = g3 * zaersl * 2.0;  // HOCl+HCl (r174)
    s.ratek[174] = g7 * zaersl * 1.5;  // HOBr+HCl (r175)
    s.ratek[175] = _g8 * zaersl * 1.5; // HOBr+HBr (r176, =0)
    s.ratek[176] = g6 * zaersl;        // BrONO2+H2O (r177)
}

// ── FKNASA — 3-body (fall-off) rate constant ─────────────────────────────────

/// NASA/JPL Troe formula for three-body reactions.
/// Fortran: FUNCTION FKNASA(T300LG, DENSTY, II)
fn fknasa(t300lg: f64, densty: f64, ii: usize, s: &ModelState) -> f64 {
    let r3 = s.rk[[0, ii]] * densty * (-t300lg * s.rk[[1, ii]]).exp();
    let ndx = s.ndxrat[ii] as usize;
    let r2  = s.rkadd[[0, ndx]] * (-t300lg * s.rkadd[[1, ndx]]).exp();
    let ratio = r3 / r2;
    (r3 / (1.0 + ratio))
        * (s.rkadd[[0, ndx + 1]] / (1.0 + ratio.log10().powi(2))).exp()
}

// ── CHEMS — compute all reaction rates ───────────────────────────────────────

/// Compute all 221 reaction rates (s.r[0..220]) from current species densities.
/// Fortran: SUBROUTINE CHEMS
pub fn chems(s: &mut ModelState) {
    let iv = s.ibox; // 0-based; J-values loaded per box
    let ib = s.ibox;

    // Recompute RATEK if altitude changed
    if s.ialt != s.izalt {
        setupr(s);
    } else {
        setupr(s); // always recompute (matches Fortran with GOTO 2 commented out)
    }

    // Local aliases for species densities (all use 1-based N-indices)
    let n = s.n; // copy so we can borrow s again
    let xno    = xr(s, n[0]);
    let xno2   = xr(s, n[1]);
    let xno3   = xr(s, n[2]);
    let xn2o5  = xr(s, n[3]);
    let xhno3  = xr(s, n[4]);
    let xhno2  = xr(s, n[14]);
    let xhno4  = xr(s, n[20]);
    let xh     = xr(s, n[5]);
    let xoh    = xr(s, n[6]);
    let xho2   = xr(s, n[7]);
    let xh2o2  = xr(s, n[8]);
    let xo     = xr(s, n[9]);
    let xo3    = xr(s, n[10]);
    let xnum   = s.zdnum;
    let xo2    = s.po2 * xnum;
    let xh2    = s.fh2[ib] * xnum;
    let xh2o   = s.fh2o[ib] * xnum;
    let xch4   = s.fch4[ib] * xnum;
    let _xnoy  = s.fnoy[ib] * xnum;
    let xco    = s.fco[ib] * xnum;
    let xco2   = s.pco2 * xnum;
    let _xn2o  = s.fn2o[ib] * xnum;
    let xocs   = s.focs[ib] * xnum;
    let xtcl   = s.fclx[ib] * xnum;
    let _xcf2cl  = s.fcf2cl[ib] * xnum;
    let _xcfcl3  = s.fcfcl3[ib] * xnum;
    let _xccl4   = s.fccl4[ib] * xnum;
    let _xch3cl  = s.fch3cl[ib] * xnum;
    let _xmecl   = s.fmecl[ib] * xnum;
    let xhcl   = xr(s, n[15]);
    let xcl    = xr(s, n[16]);
    let xcl2   = xr(s, n[17]);
    let xclo   = xr(s, n[18]);
    let xhocl  = xr(s, n[21]);
    let xclno3 = xr(s, n[19]);
    let xh2co  = xr(s, n[24]);
    let xroo   = xr(s, n[25]);
    let xrooh  = xr(s, n[26]);
    let xoclo  = xr(s, n[27]);
    let xcl2o2 = xr(s, n[28]);
    let xbro   = xr(s, n[11]);
    let xbr    = xr(s, n[12]);
    let xhbr   = xr(s, n[13]);
    let xbrno3 = xr(s, n[22]);
    let xhobr  = xr(s, n[23]);
    let xbrcl  = xr(s, n[29]);
    let _xtbr  = s.fbrx[ib] * xnum;
    let _xch3br = s.fch3br[ib] * xnum;
    let xxxx   = s.fxxx[ib] * xnum;
    let mut xo1d = 1.0_f64;
    let mut xn4s = 1.0_f64;
    let _ = (xtcl, xco2);

    // Zero all rates
    s.r.iter_mut().for_each(|v| *v = 0.0);
    s.rp.iter_mut().for_each(|v| *v = 0.0);
    s.rl.iter_mut().for_each(|v| *v = 0.0);
    s.rpf.iter_mut().for_each(|v| *v = 0.0);
    s.rlf.iter_mut().for_each(|v| *v = 0.0);
    s.rqf.iter_mut().for_each(|v| *v = 1.0);

    // Macro-like helper: get J-value by Fortran 1-based field index from v-arrays
    let jno    = s.vno[iv];
    let jvo2   = s.vo2[iv];
    let jvo3   = s.vo3[iv];
    let jvo3d  = s.vo3d[iv];
    let jvh2coa = s.vh2coa[iv];
    let jvh2cob = s.vh2cob[iv];
    let jvh2o2  = s.vh2o2[iv];
    let jvrooh  = s.vrooh[iv];
    let jvno2   = s.vno2[iv];
    let jvno3x  = s.vno3x[iv];
    let jvno3l  = s.vno3l[iv];
    let jvn2o5  = s.vn2o5[iv];
    let jvhno2  = s.vhno2[iv];
    let jvhno3  = s.vhno3[iv];
    let jvhno4  = s.vhno4[iv];
    let jvclno3 = s.vclno3[iv];
    let jvcl2   = s.vcl2[iv];
    let jvhocl  = s.vhocl[iv];
    let jvoclo  = s.voclo[iv];
    let jvcl2o2 = s.vcl2o2[iv];
    let jvclo   = s.vclo[iv];
    let jvbro   = s.vbro[iv];
    let jvbrno3 = s.vbrno3[iv];
    let jvhobr  = s.vhobr[iv];
    let jvn2o   = s.vn2o[iv];
    let jvcfcl3 = s.vcfcl3[iv];
    let jvf2cl2 = s.vf2cl2[iv];
    let jvf113  = s.vf113[iv];
    let jvf114  = s.vf114[iv];
    let jvf115  = s.vf115[iv];
    let jvccl4  = s.vccl4[iv];
    let jvch3cl = s.vch3cl[iv];
    let jvmecf  = s.vmecf[iv];
    let jvch3br = s.vch3br[iv];
    let jvh1211 = s.vh1211[iv];
    let jvh1301 = s.vh1301[iv];
    let jvh2402 = s.vh2402[iv];
    let jvh22   = s.vh22[iv];
    let jvh123  = s.vh123[iv];
    let jvh141b = s.vh141b[iv];
    let jvchbr3 = s.vchbr3[iv];
    let jvch3i  = s.vch3i[iv];
    let jvcf3i  = s.vcf3i[iv];
    let jvocs   = s.vocs[iv];
    let _ = (jvchbr3, jvch3i, jvcf3i); // not used in CHEMS directly

    // ── Reaction rates (0-based: r[0]=R(1), r[1]=R(2), ...) ─────────────────
    let r = &mut s.r;

    r[0]   = jvo2 * xo2;
    r[1]   = jvo3d * xo3;
    r[2]   = jvo3 * xo3;
    r[3]   = s.ratek[3] * xo1d;
    r[4]   = s.ratek[4] * xh2o * xo1d;
    r[5]   = s.ratek[5] * xo1d * xh2;
    r[6]   = s.ratek[6] * xo1d * xch4;
    r[7]   = s.ratek[7] * xo1d * xch4;
    r[8]   = s.ratek[8] * xo1d * s.fn2o[ib] * xnum;
    r[9]   = s.ratek[9] * xo1d * s.fn2o[ib] * xnum;
    r[10]  = s.ratek[10] * xo * xo2;
    r[11]  = s.ratek[11] * xo * xo3;
    r[12]  = s.ratek[12] * xo * xo;
    r[13]  = s.ratek[13] * xo * xh2;
    // R15: H2O+hv approximation
    r[14]  = s.xjdo;
    r[15]  = s.ratek[15] * xo * xoh;
    r[16]  = s.ratek[16] * xo * xho2;
    r[17]  = s.ratek[17] * xo * xh2o2;
    r[18]  = s.ratek[18] * xo * xno2;
    r[19]  = jvno2 * xno2;
    r[20]  = s.ratek[20] * xo3 * xno;
    r[21]  = s.ratek[21] * xo1d * xo3;
    r[22]  = s.ratek[22] * xo * xno;
    r[23]  = s.ratek[23] * xo * xno2;
    r[24]  = s.ratek[24] * xh * xo3;
    r[25]  = s.ratek[25] * xoh * xo3;
    r[26]  = s.ratek[26] * xo3 * xho2;
    r[27]  = s.ratek[27] * xo3 * xno2;
    r[28]  = s.ratek[28] * xo * xocs;  // OCS
    r[29]  = s.ratek[29] * xoh * xocs; // OCS
    r[30]  = jvh2o2 * xh2o2;
    r[31]  = s.ratek[31] * xh * xo2;
    r[32]  = s.ratek[32] * xh * xho2;
    r[33]  = s.ratek[33] * xh * xho2;
    r[34]  = s.ratek[34] * xo1d;
    r[35]  = s.ratek[35] * xo * xno3;
    // r[36] = R(37) = O(1D) density — set after steady-state below
    r[37]  = s.ratek[37];
    r[38]  = s.ratek[38] * xoh * xoh;
    r[39]  = s.ratek[39] * xoh * xho2;
    r[40]  = s.ratek[40] * xoh * xh2o2;
    r[41]  = s.ratek[41] * xho2 * xno;
    r[42]  = s.ratek[42] * xho2;
    r[43]  = s.ratek[43] * xho2 * xho2;
    r[44]  = s.ratek[44] * xh2o2;
    r[45]  = s.ratek[45] * xh2o2;
    r[46]  = s.ratek[46] * xoh * xh2;
    r[47]  = s.ratek[47] * xoh * xco;
    r[48]  = s.ratek[48] * xoh * xch4;
    // r[49] = R(50) — unused
    r[50]  = s.ratek[50] * xroo * xno;
    r[51]  = s.ratek[51] * xroo * xho2;
    r[52]  = s.ratek[52] * xroo * xroo;
    r[53]  = jvrooh * xrooh;
    r[54]  = s.ratek[54] * xoh * xrooh;
    r[55]  = s.ratek[55] * xrooh;
    r[56]  = s.ratek[56] * xrooh;
    r[57]  = s.ratek[57] * xoh * xh2co;
    r[58]  = jvh2coa * xh2co;
    r[59]  = jvh2cob * xh2co;
    r[60]  = s.ratek[60] * xh2co;
    r[61]  = s.ratek[61] * xoh * xno2 * 1.2;
    r[62]  = jvhno3 * xhno3;
    r[63]  = s.ratek[63] * xoh * xhno3;
    r[64]  = s.ratek[64] * xhno3;
    // r[65], r[66] unused
    r[66]  = s.ratek[66] * xoh * xno;
    r[67]  = s.ratek[67] * xho2 * xno2;
    r[68]  = jvhno2 * xhno2;
    r[69]  = s.ratek[69] * xoh * xhno2;
    r[70]  = s.ratek[70] * xhno2;
    r[71]  = s.ratek[71] * xho2 * xno2;
    r[72]  = s.ratek[72] * xhno4;
    r[73]  = jvhno4 * xhno4;
    r[74]  = s.ratek[74] * xhno4;
    r[75]  = s.ratek[75] * xoh * xhno4;
    r[76]  = jvno3x * xno3;
    r[77]  = jvno3l * xno3;
    r[78]  = s.ratek[78] * xno * xno3;
    r[79]  = s.ratek[79] * xno2 * xno3;
    r[80]  = s.ratek[80];
    r[81]  = s.ratek[81] * xno2 * xno3;
    r[82]  = s.ratek[82] * xn2o5;
    r[83]  = jvn2o5 * xn2o5;
    r[84]  = s.ratek[84] * xno2;
    r[85]  = jno * xno;
    r[86]  = s.ratek[86] * xn4s * xo2;
    r[87]  = s.ratek[87] * xn4s * xo3;
    r[88]  = s.ratek[88] * xn4s * xno;
    r[89]  = s.ratek[89] * xn4s * xno2;
    r[93]  = s.ratek[93];
    r[94]  = s.ratek[94] * xoh * xxxx;
    r[95]  = s.ratek[95] * xxxx;
    r[96]  = s.ratek[96] * xcl * xxxx;
    r[99]  = s.ratek[99] * xoh * xhcl;
    r[100] = s.ratek[100] * xo * xhcl;
    r[101] = s.ratek[101] * xhcl * xclno3;
    r[102] = s.ratek[102] * xcl * xho2;
    r[103] = s.ratek[103] * xcl * xch4;
    r[104] = s.ratek[104] * xcl * xhno4;
    r[105] = s.ratek[105] * xcl * xh2co;
    r[106] = s.ratek[106] * xcl * xh2;
    r[107] = s.ratek[107] * xcl * xho2;
    r[108] = s.ratek[108] * xcl * xh2o2;
    r[109] = s.ratek[109] * xo1d * s.fccl4[ib] * xnum;
    r[110] = s.ratek[110] * xo1d * s.fcfcl3[ib] * xnum;
    r[111] = s.ratek[111] * xo1d * s.fcf2cl[ib] * xnum;
    r[112] = s.ratek[112] * xo1d * xhcl;
    r[113] = s.ratek[113] * xo3 * xcl;
    r[114] = s.ratek[114] * xo * xclo;
    r[115] = s.ratek[115] * xclo * xno;
    r[116] = jvclo * xclo;
    r[117] = s.ratek[117] * xclo * xclo;
    r[118] = s.ratek[118] * xcl2o2;
    r[119] = s.ratek[119] * xoh * xclo;
    r[120] = s.ratek[120] * xho2 * xclo;
    r[121] = s.ratek[121] * xclo * xno2;
    r[122] = s.ratek[122] * xclno3;
    r[123] = s.ratek[123] * jvclno3 * xclno3;
    r[124] = s.ratek[124] * jvclno3 * xclno3;
    r[125] = s.ratek[125] * xclno3 * xo;
    r[126] = s.ratek[126] * xclno3 * xoh;
    r[127] = s.ratek[127] * xho2 * xclo;
    r[128] = jvhocl * xhocl;
    r[129] = s.ratek[129] * xoh * xhocl;
    r[136] = jvn2o * s.fn2o[ib] * xnum;
    r[138] = jvf2cl2 * s.fcf2cl[ib] * xnum;
    r[139] = jvcfcl3 * s.fcfcl3[ib] * xnum;
    r[140] = jvccl4 * s.fccl4[ib] * xnum;
    r[141] = jvch3cl * s.fch3cl[ib] * xnum;
    r[142] = s.ratek[142] * xoh * s.fch3cl[ib] * xnum;
    r[143] = jvmecf * s.fmecl[ib] * xnum;
    r[144] = s.ratek[144] * xoh * s.fmecl[ib] * xnum;
    r[145] = jvocs * xocs;  // OCS photolysis
    r[146] = jvch3br * s.fch3br[ib] * xnum;
    r[147] = s.ratek[147] * s.fch3br[ib] * xnum * xoh;
    r[148] = s.ratek[148] * xbrno3 * xo;
    r[149] = s.ratek[149] * xbro * xoh;
    r[156] = s.ratek[156] * xoh * xrooh;
    r[157] = jvcl2o2 * xcl2o2;
    r[158] = s.ratek[158] * xcl * xcl;
    r[159] = jvcl2 * xcl2;
    r[160] = s.ratek[160] * xclo * xclo;
    r[161] = jvoclo * xoclo;
    r[162] = s.ratek[162] * xoclo * xno;
    r[163] = s.ratek[163] * xoclo * xoh;
    r[164] = s.ratek[164] * xoclo * xcl;
    r[165] = s.ratek[165] * xcl * xclno3;
    r[166] = s.ratek[166] * xcl * xno3;
    r[167] = s.ratek[167];
    // Heterogeneous (170–177, 0-based 169–176)
    r[169] = s.ratek[169] * xn2o5;
    r[170] = s.ratek[170] * xclno3;
    r[171] = s.ratek[171] * xclno3;
    r[172] = s.ratek[172] * xn2o5;
    r[173] = s.ratek[173] * xhocl;
    r[174] = s.ratek[174] * xhobr;
    r[175] = s.ratek[175] * xhobr;
    r[176] = s.ratek[176] * xbrno3;
    r[177] = s.ratek[177];
    r[178] = s.ratek[178];
    // F-113/114/115 photolysis (loss freqs only — not density-weighted)
    r[179] = jvf113;
    r[180] = s.ratek[180] * xo1d;
    r[181] = jvf114;
    r[182] = s.ratek[182] * xo1d;
    r[183] = jvf115;
    r[184] = s.ratek[184] * xo1d;
    r[185] = jvh1211;
    r[186] = jvh1301;
    r[187] = jvh2402;
    r[188] = jvh22;
    r[189] = s.ratek[189] * xo1d;
    r[190] = jvh123;
    r[191] = jvh141b;
    r[192] = s.ratek[192] * xoh;
    r[193] = s.ratek[193] * xoh;
    r[194] = s.ratek[194] * xoh;
    r[195] = s.ratek[195] * xo1d;
    r[196] = s.ratek[196] * xo1d;
    r[197] = s.ratek[197] * xo1d;
    r[198] = s.ratek[198] * xclo * xclo;
    r[199] = s.ratek[199] * xclo * xclo;

    // Bromine reactions (200–220)
    if s.lbrom {
        r[200] = s.ratek[200] * xhbr * xoh;
        r[201] = s.ratek[201] * xhbr * xo;
        r[203] = s.ratek[203] * xbr * xho2;
        r[204] = s.ratek[204] * xbr * xo3;
        r[205] = s.ratek[205] * xbro * xo;
        r[206] = s.ratek[206] * xbro * xno;
        r[207] = s.ratek[207] * xbro * xo3;
        r[208] = s.ratek[208] * xbro * xbro;
        r[209] = jvbro * xbro;
        r[210] = s.ratek[210] * xbro * xho2;
        r[211] = jvhobr * xhobr;
        r[212] = s.ratek[212] * xoh * xhobr;
        r[213] = s.ratek[213] * xbro * xno2;
        r[214] = s.ratek[214] * jvbrno3 * xbrno3;
        r[215] = s.ratek[215] * jvbrno3 * xbrno3;
        r[216] = s.ratek[216] * xclo * xbro;
        r[217] = s.ratek[217] * xbr * xh2co;
        r[218] = s.ratek[218] * xclo * xbro;
        r[219] = s.ratek[219] * jvcl2 * xbrcl; // VBrCl scaled to VCl2
        r[220] = s.ratek[220] * xclo * xbro;
    }

    // ── O(1D) instant steady-state ───────────────────────────────────────────
    let tempp = r[1]; // R(2) = jO3→O1D * O3
    let templ = r[3] + r[4] + r[5] + r[6] + r[7] + r[8] + r[9] + r[21] + r[112];
    xo1d = if templ > 0.0 { tempp / templ } else { 0.0 };
    r[36] = xo1d; // R(37)
    // Rescale O(1D)-dependent rates
    for i in [3usize, 4, 5, 6, 7, 8, 9, 21, 112] { r[i] *= xo1d; }
    for i in [34usize, 109, 110, 111, 180, 182, 184, 189, 195, 196, 197] { r[i] *= xo1d; }

    // ── N(4S) instant steady-state ────────────────────────────────────────────
    let tempp = r[85]; // R(86) = J(NO)*NO
    let templ = r[86] + r[87] + r[88] + r[89];
    xn4s = if templ > 0.0 { tempp / templ } else { 0.0 };
    r[90] = xn4s; // R(91)
    for i in [86usize, 87, 88, 89] { r[i] *= xn4s; }

    chempl(s);
}

// ── CHEMPL — production/loss for each implicit species ───────────────────────

/// Compute s.rp, s.rl, s.rpf, s.rlf from s.r.
/// Fortran: SUBROUTINE CHEMPL
pub fn chempl(s: &mut ModelState) {
    let n  = s.n;
    let r  = s.r;  // copy (250 × f64 = 2 kB) to avoid borrow conflicts
    let ib = s.ibox;

    // Intermediates
    let rpch3o = r[50] + r[52] + r[52] + r[53];
    let rpcho  = r[57] + r[58] + r[105] + r[217];

    let rp = &mut s.rp;
    let rl = &mut s.rl;

    // Indexed by Fortran 1-based N-index (converted with idx())

    // O
    rp[idx(n[9])]  = r[0]+r[0]+r[19]+r[38]+r[76]+r[116]+r[209]+r[161]
                    +r[2]-r[4]-r[5]-r[6]-r[7]-r[8]-r[9]-r[21]-r[112];
    rl[idx(n[9])]  = r[12]+r[12]+r[15]+r[16]+r[18]+r[10]+r[11]+r[35]
                    +r[22]+r[23]+r[17]+r[13]+r[100]+r[114]+r[125]+r[201]+r[205]+r[148];

    // CH3O2 = ROO
    rp[idx(n[25])] = r[6]+r[48]+r[54]+r[103]+r[141]+r[143]+r[146];
    rl[idx(n[25])] = r[50]+r[51]+r[52]+r[52];

    // CH3O2H = ROOH
    rp[idx(n[26])] = r[51] + r[93];
    rl[idx(n[26])] = r[55]+r[54]+r[56]+r[53]+r[156];

    // H2CO
    rp[idx(n[24])] = r[7]+r[156] + rpch3o + r[80];
    rl[idx(n[24])] = r[57]+r[58]+r[59]+r[60]+r[105]+r[217];

    s.r[130] = rpcho + rpch3o + r[6] + r[58] + r[53] - r[48] - r[54] - r[51] - r[57];

    // H2O2
    rp[idx(n[8])]  = r[43] + r[37];
    rl[idx(n[8])]  = r[30]+r[17]+r[40]+r[44]+r[45]+r[108];

    // HNO2
    rp[idx(n[14])] = r[66]+r[67];
    rl[idx(n[14])] = r[68]+r[69] + r[70];

    // HNO3
    rp[idx(n[4])]  = r[61]+r[101] +r[169]+r[169]+r[170]+r[171]+r[172]+r[176];
    rl[idx(n[4])]  = r[62]+r[63] + r[64];

    // HNO4
    rp[idx(n[20])] = r[71];
    rl[idx(n[20])] = r[73]+r[75]+r[104]+r[72] + r[74];

    // N2O5
    rp[idx(n[3])]  = r[81];
    rl[idx(n[3])]  = r[83]+r[82] +r[169]+r[172];

    // ClONO2
    rp[idx(n[19])] = r[121];
    rl[idx(n[19])] = r[122]+r[123]+r[124]+r[125]+r[126]+r[101]+r[165]+r[170]+r[171];

    // HOCl
    rp[idx(n[21])] = r[127]+r[126]+r[163] + r[170];
    rl[idx(n[21])] = r[128]+r[129] +r[173];

    // HCl
    rp[idx(n[15])] = r[106]+r[108]+r[107]+r[103]+r[105]+r[104]
                    +r[96]+r[119]+r[120];
    rl[idx(n[15])] = r[99]+r[101]+r[100]+r[112]+r[171]+r[172]+r[173] +r[174];

    // ClO
    rp[idx(n[18])] = r[113]+r[118]+r[118]+r[122]+r[123]+r[125]+r[129]
                    +r[161]+r[162]+r[164]+r[164]+r[166]+r[102];
    rl[idx(n[18])] = r[114]+r[115]+r[116]+r[117]+r[117]+r[119]+r[120]+r[121]
                    +r[127]+r[216]+r[218]+r[160]+r[160]+r[220]
                    +2.0*(r[198]+r[199]);

    // Cl
    rp[idx(n[16])] = r[99]+r[100]+r[115]+r[116]+r[114]+r[124]+r[220]+r[172]
                    +r[219]+r[112]+r[128]+r[157]+r[157]+r[159]+r[159]+r[160]
                    +2.0*r[198];
    rl[idx(n[16])] = r[96]+r[103]+r[104]+r[105]+r[106]+r[107]+r[108]+r[113]
                    +r[158]+r[158]+r[164]+r[165]+r[166]+r[102];

    // OClO
    rp[idx(n[27])] = r[218]+r[160];
    rl[idx(n[27])] = r[161]+r[162]+r[163]+r[164];

    // Cl2O2
    rp[idx(n[28])] = r[117];
    rl[idx(n[28])] = r[118]+r[157];

    // Cl2
    rp[idx(n[17])] = r[158]+r[101]+r[165] +r[171]+r[173] + r[199];
    rl[idx(n[17])] = r[159];

    if s.lbrom {
        // HBr
        rp[idx(n[13])] = r[203]+r[217];
        rl[idx(n[13])] = r[200]+r[201] +r[175];

        // BrONO2
        rp[idx(n[22])] = r[213];
        rl[idx(n[22])] = r[215]+r[214] +r[176]+r[148];

        // HOBr
        rp[idx(n[23])] = r[210] +r[176];
        rl[idx(n[23])] = r[211]+r[212] +r[174]+r[175];

        // BrO
        rp[idx(n[11])] = r[204]+r[212]+r[214] +r[148];
        rl[idx(n[11])] = r[205]+r[206]+r[213]+2.0*r[208]+r[207]+r[210]+r[216]
                        +r[209]+r[218]+r[220] +r[149];

        // Br
        rp[idx(n[12])] = r[205]+r[206]+2.0*r[208]+r[207]+r[200]+r[201]
                        +r[211]+r[209]+r[215]+r[218]+r[219]+r[220]
                        +2.0*r[175] +r[149];
        rl[idx(n[12])] = r[204]+r[203]+r[217];

        // BrCl
        rp[idx(n[29])] = r[216] +r[174];
        rl[idx(n[29])] = r[219];
    }

    // NO
    rp[idx(n[0])]  = r[79]+r[77]+r[19]+r[18]+r[68]     +r[177];
    rl[idx(n[0])]  = r[22]+r[78]+r[41]+r[20]+r[66]+r[50]+r[206]+r[115]+r[162];

    // NO2
    rp[idx(n[1])]  = r[62]+r[22]+r[78]+r[78]+r[76]+r[83]+r[82]+r[122]
                    +r[35]+r[41]+r[20]+r[69]+r[50]+r[75]+r[104]+r[206]
                    +r[72]+r[214]+r[115]+r[123]+r[162]+r[166]+r[172]+r[178];
    rl[idx(n[1])]  = r[27]+r[23]+r[81]+r[19]+r[18]+r[67]+r[71]+r[121]+r[213]
                    +r[61];

    // NO3
    rp[idx(n[2])]  = r[27]+r[23]+r[83]+r[82]+r[73]+r[215]+r[125]+r[126]+r[124]
                    +r[165] + r[63] +r[148];
    rl[idx(n[2])]  = r[35]+r[78]+r[79]+r[77]+r[76]+r[81]+r[166];

    // H
    rp[idx(n[5])]  = r[5]+r[15]+r[46]+r[47]+r[58]+r[13]+r[106] +r[14];
    rl[idx(n[5])]  = r[24]+r[33]+r[32]+r[31];

    // OH
    rp[idx(n[6])]  = r[5]+r[24]+2.0*(r[32]+r[30]+r[4])
                    +r[17]+r[16]+r[62]+r[41]+r[26]
                    +r[6]+r[13]+r[68]+r[53]+r[73]
                    +r[100]+r[201]+r[211]+r[128]+r[112] + r[42]+r[102]+r[14];
    rl[idx(n[6])]  = r[15]+r[46]+r[47]+r[25]+r[38]+r[38]+r[40]
                    +r[39]+r[61]+r[63]+r[48]+r[54]+r[57]+r[66]+r[69]
                    +r[94]+r[119]+r[75]+r[200]+r[212]+r[147]
                    +r[99]+r[142]+r[144]+r[126]+r[129]+r[163] +r[149];

    // HO2
    rp[idx(n[7])]  = r[31]+r[17]+r[25]+r[40]+r[72]+r[108] +r[149]
                    +rpch3o+rpcho;
    rl[idx(n[7])]  = r[33]+r[32]+r[16]+r[41]+r[39]+r[43]+r[43]+r[26]
                    +r[67]+r[51]+r[71]+r[42]+r[210]+r[203]+r[127]+r[107]+r[120]
                    +r[102];

    // O3
    rp[idx(n[10])] = r[10]+r[167]+r[120];
    rl[idx(n[10])] = r[24]+r[25]+r[2]+r[11]+r[20]+r[26]+r[27]+r[113]+r[204]
                    +r[207]+r[21];

    // ── Diagnostic rates ─────────────────────────────────────────────────────

    // OddH P=R(132), L=R(133)
    s.r[131] = 2.0*(r[5]+r[4]+r[13])+r[62]+r[6]+r[58]+r[72]
              +rpch3o+rpcho+r[68]+r[53]+r[73]+r[106]+r[100]
              +r[201]+r[211]+r[128] +2.0*r[14];
    s.r[132] = 2.0*(r[33]+r[38]+r[40]+r[39]+r[45]+r[44])
              +r[61]+r[63]+r[48]+r[57]+r[66]+r[69]+r[67]+r[51]+r[71]
              +r[54]+r[75]+      r[210]+r[129]+r[212]
              +r[99]+r[200]+r[203]+r[147]+r[127]+r[108]+r[94]
              +r[107]+r[120]+r[142]+r[144]+r[126]+r[163]+r[119];

    // OddO P=R(135), L=R(136)
    let (r135, r136) = {
        let r = &s.r;
        let mut r135 = r[38]+r[0]+r[0]+r[76]+r[116]+r[209]+r[19]+r[161]+r[167]+r[120];
        let mut r136 = r[5]+r[15]+r[24]+r[4]+r[17]+r[16]+r[21]+r[21]+r[25]+r[18]
                      +r[35]+r[11]+r[11]+r[12]+r[12]+r[26]+r[7]+r[6]+r[23]+r[8]
                      +r[9]+r[204]+r[205]+r[207]+r[100]+r[113]+r[201]+r[114]
                      +r[112]+r[125]+r[13]+r[20]+r[27]+r[22];
        let trupox = r[0]+r[0]+r[41]+r[50]+r[73]+r[63] + r[167];
        r136 = r136 - r135 + trupox;
        r135 = trupox;
        (r135, r136)
    };
    s.r[134] = r135;
    s.r[135] = r136;

    // P-L O3 (mixing ratio/day) = R(169)
    let ialt = s.ialt;
    s.r[168] = (r135 - r136) * 86400.0 / s.dm[ialt];

    // Cl/Br catalytic cycles
    s.r[133] = 2.0*(r[114]+r[157]+r[198]+r[199]) +r[100]+r[119]+r[128]
              +r[163]+r[216]+r[220];
    s.r[137] = 2.0*(r[205]+r[207]+r[208]) +r[211]+r[216]+r[220];

    // ── Long-lived (explicit integration) P-L ────────────────────────────────
    let rpf = &mut s.rpf;
    let rlf = &mut s.rlf;

    // 1 = O3
    rpf[0] = r135;
    rlf[0] = r136;
    // 2 = N2O
    rpf[1] = r[89]+r[34];
    rlf[1] = r[8]+r[9]+r[136];
    // 3 = NOy
    rpf[2] = r[8]+r[8] +r[177]+r[178];
    rlf[2] = r[64]+r[70]+r[74]+r[84]+r[88]+r[88]+r[89]+r[89];
    // 4 = CH4
    rlf[3] = r[7]+r[6]+r[48]+r[103];
    // 5 = CO
    rpf[4] = r[57]+r[58]+r[59];
    rlf[4] = r[47];
    // 6 = CLX
    rpf[5] = 2.*r[138]+3.*r[139]+4.*r[140]+r[141]+r[142]+3.*r[143]+3.*r[144];
    // 7 = CF2Cl2
    rlf[6] = r[138]+r[111];
    // 8 = CFCl3
    rlf[7] = r[139]+r[110];
    // 9 = CCl4
    rlf[8] = r[140]+r[109];
    // 10 = CH3Cl
    rlf[9] = r[142]+r[141];
    // 11 = MeCl = CH3CCl3
    rlf[10] = r[144]+r[143];
    // 12 = H2
    rpf[11] = r[7]+r[33]+r[59];
    rlf[11] = r[5]+r[13]+r[46]+r[106];
    // 13 = H2O
    rpf[12] = r[5]+r[46]+r[13]+r[106]+r[6]+r[48]+r[103]
             +r[57]+r[58] + 1.5*(rlf[9]+rlf[10]);
    rlf[12] = r[33];
    // 16 = total inorganic Br
    rpf[15] = r[147]+r[146];
    // 17 = CH3Br
    rlf[16] = r[147]+r[146];
    // 18 = XXX (C2H6)
    rlf[17] = r[94]+r[95]+r[96];
    // 19 = H2O w/o CH4
    rpf[18] = r[5]+r[46]+r[13]+r[106]+r[57]+r[58] + 1.5*(rlf[9]+rlf[10]);
    rlf[18] = r[33];
    // 20 = H2O literal
    rpf[19] = r[38]+r[39]+r[40]+r[46]+r[48]+r[54]+r[57]+r[63]+r[69]
             +r[75]+r[99]+r[129]+r[200]+r[212]+r[173]+r[174]+r[175];
    rlf[19] = r[4]+r[169]+r[176]+r[170] +r[14];
    // 21 = OCS
    rpf[20] = 0.0;
    rlf[20] = r[28] + r[29] + r[145];
}

// ── RHSLHS — RHS vector and Jacobian ─────────────────────────────────────────

/// Build either:
///   ido=0: RHS vector fxo[i] = RL[i] - RP[i] + (XR[i] - XNOLD[i]) * DELTT
///   ido=1: Jacobian A (NTOT×NTOT) via REACT calls
/// Fortran: SUBROUTINE RHSLHS(FXO, IDO)
pub fn rhslhs(s: &mut ModelState, ido: i32) -> DVector<f64> {
    let ntot = s.ntot;
    let n    = s.n;

    // Closure species row indices (Fortran 1-based; used as rows in A)
    let innoy = n[1];  // N2
    let inclx = n[18]; // N19 (ClO)
    let inbrx = n[11]; // N12 (BrO)

    if ido == 0 {
        // Steady-state closure when DELTT < 1e-20
        if s.deltt < 1.0e-20 {
            let zdnum = s.zdnum;
            let ib    = s.ibox;
            s.rp[idx(innoy)] = s.fnoy[ib] * zdnum;
            s.rl[idx(innoy)] = (0..ntot)
                .filter(|&i| [n[0],n[1],n[2],2*n[3],n[4],n[14],n[20],n[19],n[22]]
                    .contains(&(i+1)))
                .map(|i| s.xr[i])
                .sum::<f64>();
            // Simplified: sum the listed species
            s.rl[idx(innoy)] = s.xr[idx(n[0])] + s.xr[idx(n[1])] + s.xr[idx(n[2])]
                + 2.0*s.xr[idx(n[3])] + s.xr[idx(n[4])] + s.xr[idx(n[14])]
                + s.xr[idx(n[20])] + s.xr[idx(n[19])] + s.xr[idx(n[22])];
            s.rp[idx(inclx)] = s.fclx[ib] * zdnum;
            s.rl[idx(inclx)] = s.xr[idx(n[15])] + s.xr[idx(n[16])]
                + s.xr[idx(n[18])] + s.xr[idx(n[19])] + s.xr[idx(n[21])]
                + s.xr[idx(n[27])]
                + 2.0*(s.xr[idx(n[17])] + s.xr[idx(n[28])]);
            s.rp[idx(inbrx)] = s.fbrx[ib] * zdnum;
            s.rl[idx(inbrx)] = s.xr[idx(n[11])] + s.xr[idx(n[12])]
                + s.xr[idx(n[13])] + s.xr[idx(n[22])] + s.xr[idx(n[23])];
        }

        let mut fxo = DVector::zeros(ntot);
        let deltt = s.deltt;
        for i in 0..ntot {
            fxo[i] = s.rl[i] - s.rp[i] + (s.xr[i] - s.xnold[i]) * deltt;
        }
        fxo
    } else {
        // Build Jacobian A(NTOT×NTOT)
        let a = &mut s.a_mat;
        for j in 0..ntot {
            for i in 0..ntot {
                a[[i, j]] = 0.0;
            }
            a[[j, j]] = -s.deltt * s.xr[j];
        }

        let r = s.r;
        let n = s.n;

        // Each REACT call: react(s, rrate, nl1, nl2, nl3, np1, np2) with 1-based indices
        macro_rules! rct {
            ($rt:expr, $l1:expr, $l2:expr, $l3:expr, $p1:expr, $p2:expr) => {
                react(&mut s.a_mat, $rt, $l1, $l2, $l3, $p1, $p2)
            };
        }

        let r5to10 = r[4]+r[5]+r[6]+r[7]+r[8]+r[9]+r[21]+r[112];
        rct!(r[2]-r5to10, n[10] as i32, 0, 0, n[9] as i32, 0);
        rct!(r[4], n[10] as i32, 0, 0, n[6] as i32, n[6] as i32);
        rct!(r[5], n[10] as i32, 0, 0, n[6] as i32, n[5] as i32);
        rct!(r[6], n[10] as i32, 0, 0, n[6] as i32, n[25] as i32);
        rct!(r[7], n[10] as i32, 0, 0, n[24] as i32, 0);
        rct!(r[8], n[10] as i32, 0, 0, 0, 0);
        rct!(r[9], n[10] as i32, 0, 0, 0, 0);
        rct!(r[10], n[9] as i32, 0, 0, n[10] as i32, 0);
        rct!(r[11], n[9] as i32, n[10] as i32, 0, 0, 0);
        rct!(r[12], n[9] as i32, n[9] as i32, 0, 0, 0);
        rct!(r[13], n[9] as i32, 0, 0, n[5] as i32, n[6] as i32);
        rct!(r[15], n[9] as i32, n[6] as i32, 0, n[5] as i32, 0);
        rct!(r[16], n[9] as i32, n[7] as i32, 0, n[6] as i32, 0);
        rct!(r[17], n[9] as i32, n[8] as i32, 0, n[6] as i32, n[7] as i32);
        rct!(r[18], n[9] as i32, n[1] as i32, 0, n[0] as i32, 0);
        rct!(r[19], n[1] as i32, 0, 0, n[0] as i32, n[9] as i32);
        rct!(r[20], n[10] as i32, n[0] as i32, 0, n[1] as i32, 0);
        rct!(r[21], n[10] as i32, n[10] as i32, 0, 0, 0);
        rct!(r[22], n[9] as i32, n[0] as i32, 0, n[1] as i32, 0);
        rct!(r[23], n[9] as i32, n[1] as i32, 0, n[2] as i32, 0);
        rct!(r[24], n[10] as i32, n[5] as i32, 0, n[6] as i32, 0);
        rct!(r[25], n[10] as i32, n[6] as i32, 0, n[7] as i32, 0);
        rct!(r[26], n[10] as i32, n[7] as i32, 0, n[6] as i32, 0);
        rct!(r[27], n[10] as i32, n[1] as i32, 0, n[2] as i32, 0);
        rct!(r[28], n[9] as i32, 0, 0, 0, 0); // OCS
        rct!(r[29], n[6] as i32, 0, 0, 0, 0); // OCS
        rct!(r[30], n[8] as i32, 0, 0, n[6] as i32, n[6] as i32);
        rct!(r[31], n[5] as i32, 0, 0, n[7] as i32, 0);
        rct!(r[32], n[5] as i32, n[7] as i32, 0, n[6] as i32, n[6] as i32);
        rct!(r[33], n[5] as i32, n[7] as i32, 0, 0, 0);
        rct!(r[35], n[9] as i32, n[2] as i32, 0, n[1] as i32, 0);
        rct!(r[38], n[6] as i32, n[6] as i32, 0, n[9] as i32, 0);
        rct!(r[39], n[6] as i32, n[7] as i32, 0, 0, 0);
        rct!(r[40], n[6] as i32, n[8] as i32, 0, n[7] as i32, 0);
        rct!(r[41], n[7] as i32, n[0] as i32, 0, n[6] as i32, n[1] as i32);
        rct!(r[42], n[7] as i32, 0, 0, n[6] as i32, 0);
        rct!(r[43], n[7] as i32, n[7] as i32, 0, n[8] as i32, 0);
        rct!(r[44], n[8] as i32, 0, 0, 0, 0);
        rct!(r[45], n[8] as i32, 0, 0, 0, 0);
        rct!(r[46], n[6] as i32, 0, 0, n[5] as i32, 0);
        rct!(r[47], n[6] as i32, 0, 0, n[5] as i32, 0);
        rct!(r[48], n[6] as i32, 0, 0, n[25] as i32, 0);
        rct!(r[50], n[25] as i32, n[0] as i32, 0, n[7] as i32, n[24] as i32);
        rct!(r[50], -(n[25] as i32), -(n[0] as i32), 0, n[1] as i32, 0);
        rct!(r[51], n[25] as i32, n[7] as i32, 0, n[26] as i32, 0);
        rct!(r[52], n[25] as i32, n[25] as i32, 0, n[24] as i32, n[24] as i32);
        rct!(r[52], -(n[25] as i32), -(n[25] as i32), 0, n[7] as i32, n[7] as i32);
        rct!(r[53], n[26] as i32, 0, 0, n[6] as i32, 0);
        rct!(r[53], -(n[26] as i32), 0, 0, n[24] as i32, n[7] as i32);
        rct!(r[54], n[26] as i32, n[6] as i32, 0, n[25] as i32, 0);
        rct!(r[55], n[26] as i32, 0, 0, 0, 0);
        rct!(r[56], n[26] as i32, 0, 0, 0, 0);
        rct!(r[57], n[24] as i32, n[6] as i32, 0, n[7] as i32, 0);
        rct!(r[58], n[24] as i32, 0, 0, n[5] as i32, n[7] as i32);
        rct!(r[59], n[24] as i32, 0, 0, 0, 0);
        rct!(r[60], n[24] as i32, 0, 0, 0, 0);
        rct!(r[61], n[1] as i32, n[6] as i32, 0, n[4] as i32, 0);
        rct!(r[62], n[4] as i32, 0, 0, n[6] as i32, n[1] as i32);
        rct!(r[63], n[6] as i32, n[4] as i32, 0, n[2] as i32, 0);
        rct!(r[64], n[4] as i32, 0, 0, 0, 0);
        rct!(r[66], n[0] as i32, n[6] as i32, 0, n[14] as i32, 0);
        rct!(r[67], n[7] as i32, n[1] as i32, 0, n[14] as i32, 0);
        rct!(r[68], n[14] as i32, 0, 0, n[6] as i32, n[0] as i32);
        rct!(r[69], n[6] as i32, n[14] as i32, 0, n[1] as i32, 0);
        rct!(r[70], n[14] as i32, 0, 0, 0, 0);
        rct!(r[71], n[7] as i32, n[1] as i32, 0, n[20] as i32, 0);
        rct!(r[72], n[20] as i32, 0, 0, n[7] as i32, n[1] as i32);
        rct!(r[73], n[20] as i32, 0, 0, n[6] as i32, n[2] as i32);
        rct!(r[74], n[20] as i32, 0, 0, 0, 0);
        rct!(r[75], n[6] as i32, n[20] as i32, 0, n[1] as i32, 0);
        rct!(r[76], n[2] as i32, 0, 0, n[1] as i32, n[9] as i32);
        rct!(r[77], n[2] as i32, 0, 0, n[0] as i32, 0);
        rct!(r[78], n[2] as i32, n[0] as i32, 0, n[1] as i32, n[1] as i32);
        rct!(r[79], n[1] as i32, n[2] as i32, 0, n[1] as i32, n[0] as i32);
        rct!(r[81], n[1] as i32, n[2] as i32, 0, n[3] as i32, 0);
        rct!(r[82], n[3] as i32, 0, 0, n[1] as i32, n[2] as i32);
        rct!(r[83], n[3] as i32, 0, 0, n[1] as i32, n[2] as i32);
        rct!(r[94], n[6] as i32, 0, 0, 0, 0);
        rct!(r[96], n[16] as i32, 0, 0, n[15] as i32, 0);
        rct!(r[99], n[15] as i32, n[6] as i32, 0, n[16] as i32, 0);
        rct!(r[100], n[15] as i32, n[9] as i32, 0, n[16] as i32, n[6] as i32);
        rct!(r[101], n[15] as i32, n[19] as i32, 0, n[17] as i32, n[4] as i32);
        rct!(r[102], n[16] as i32, n[7] as i32, 0, n[18] as i32, n[6] as i32);
        rct!(r[103], n[16] as i32, 0, 0, n[15] as i32, n[25] as i32);
        rct!(r[104], n[16] as i32, n[20] as i32, 0, n[15] as i32, n[1] as i32);
        rct!(r[105], n[24] as i32, n[16] as i32, 0, n[15] as i32, n[7] as i32);
        rct!(r[106], n[16] as i32, 0, 0, n[15] as i32, n[5] as i32);
        rct!(r[107], n[16] as i32, n[7] as i32, 0, n[15] as i32, 0);
        rct!(r[108], n[16] as i32, n[8] as i32, 0, n[15] as i32, n[7] as i32);
        rct!(r[112], n[10] as i32, n[15] as i32, 0, n[6] as i32, n[16] as i32);
        rct!(r[113], n[16] as i32, n[10] as i32, 0, n[18] as i32, 0);
        rct!(r[114], n[18] as i32, n[9] as i32, 0, n[16] as i32, 0);
        rct!(r[115], n[18] as i32, n[0] as i32, 0, n[16] as i32, n[1] as i32);
        rct!(r[116], n[18] as i32, 0, 0, n[16] as i32, n[9] as i32);
        rct!(r[117], n[18] as i32, n[18] as i32, 0, n[28] as i32, 0);
        rct!(r[118], n[28] as i32, 0, 0, n[18] as i32, n[18] as i32);
        rct!(r[119], n[18] as i32, n[6] as i32, 0, n[15] as i32, 0);
        rct!(r[120], n[18] as i32, n[7] as i32, 0, n[15] as i32, n[10] as i32);
        rct!(r[121], n[18] as i32, n[1] as i32, 0, n[19] as i32, 0);
        rct!(r[122], n[19] as i32, 0, 0, n[18] as i32, n[1] as i32);
        rct!(r[123], n[19] as i32, 0, 0, n[18] as i32, n[1] as i32);
        rct!(r[124], n[19] as i32, 0, 0, n[16] as i32, n[2] as i32);
        rct!(r[125], n[19] as i32, n[9] as i32, 0, n[18] as i32, n[2] as i32);
        rct!(r[126], n[19] as i32, n[6] as i32, 0, n[21] as i32, n[2] as i32);
        rct!(r[127], n[18] as i32, n[7] as i32, 0, n[21] as i32, 0);
        rct!(r[128], n[21] as i32, 0, 0, n[6] as i32, n[16] as i32);
        rct!(r[129], n[21] as i32, n[6] as i32, 0, n[18] as i32, 0);
        rct!(r[142], n[6] as i32, 0, 0, 0, 0);
        rct!(r[144], n[6] as i32, 0, 0, 0, 0);
        rct!(r[147], n[6] as i32, 0, 0, 0, 0);
        rct!(r[148], n[22] as i32, n[9] as i32, 0, n[11] as i32, n[2] as i32);
        rct!(r[149], n[11] as i32, n[6] as i32, 0, n[12] as i32, n[7] as i32);
        rct!(r[156], n[26] as i32, n[6] as i32, 0, n[24] as i32, n[6] as i32);
        rct!(r[157], n[28] as i32, 0, 0, n[16] as i32, n[16] as i32);
        rct!(r[158], n[16] as i32, n[16] as i32, 0, n[17] as i32, 0);
        rct!(r[159], n[17] as i32, 0, 0, n[16] as i32, n[16] as i32);
        rct!(r[160], n[18] as i32, n[18] as i32, 0, n[27] as i32, 0);
        rct!(r[161], n[27] as i32, 0, 0, n[18] as i32, n[9] as i32);
        rct!(r[162], n[27] as i32, n[0] as i32, 0, n[18] as i32, n[1] as i32);
        rct!(r[163], n[27] as i32, n[6] as i32, 0, n[21] as i32, 0);
        rct!(r[164], n[27] as i32, n[16] as i32, 0, n[18] as i32, n[18] as i32);
        rct!(r[165], n[16] as i32, n[19] as i32, 0, n[17] as i32, n[2] as i32);
        rct!(r[166], n[16] as i32, n[2] as i32, 0, n[18] as i32, n[1] as i32);
        rct!(r[169], n[3] as i32, 0, 0, n[4] as i32, n[4] as i32);
        rct!(r[170], n[19] as i32, 0, 0, n[21] as i32, n[4] as i32);
        rct!(r[171], n[19] as i32, n[15] as i32, 0, n[17] as i32, n[4] as i32);
        rct!(r[172], n[3] as i32, n[15] as i32, 0, n[16] as i32, n[4] as i32);
        rct!(r[173], n[21] as i32, n[15] as i32, 0, n[17] as i32, 0);
        rct!(r[174], n[23] as i32, n[15] as i32, 0, n[29] as i32, 0);
        rct!(r[175], n[23] as i32, n[13] as i32, 0, n[12] as i32, n[12] as i32);
        rct!(r[176], n[22] as i32, 0, 0, n[23] as i32, n[4] as i32);
        rct!(r[198], n[18] as i32, n[18] as i32, 0, n[16] as i32, n[16] as i32);
        rct!(r[199], n[18] as i32, n[18] as i32, 0, n[17] as i32, 0);

        if s.lbrom {
            rct!(r[200], n[13] as i32, n[6] as i32, 0, n[12] as i32, 0);
            rct!(r[201], n[13] as i32, n[9] as i32, 0, n[12] as i32, n[6] as i32);
            rct!(r[203], n[12] as i32, n[7] as i32, 0, n[13] as i32, 0);
            rct!(r[204], n[12] as i32, n[10] as i32, 0, n[11] as i32, 0);
            rct!(r[205], n[11] as i32, n[9] as i32, 0, n[12] as i32, 0);
            rct!(r[206], n[11] as i32, n[0] as i32, 0, n[12] as i32, n[1] as i32);
            rct!(r[207], n[11] as i32, n[10] as i32, 0, n[12] as i32, 0);
            rct!(r[208], n[11] as i32, n[11] as i32, 0, n[12] as i32, n[12] as i32);
            rct!(r[209], n[11] as i32, 0, 0, n[12] as i32, n[9] as i32);
            rct!(r[210], n[11] as i32, n[7] as i32, 0, n[23] as i32, 0);
            rct!(r[211], n[23] as i32, 0, 0, n[12] as i32, n[6] as i32);
            rct!(r[212], n[23] as i32, n[6] as i32, 0, n[11] as i32, 0);
            rct!(r[213], n[11] as i32, n[1] as i32, 0, n[22] as i32, 0);
            rct!(r[214], n[22] as i32, 0, 0, n[11] as i32, n[1] as i32);
            rct!(r[215], n[22] as i32, 0, 0, n[12] as i32, n[2] as i32);
            rct!(r[216], n[11] as i32, n[18] as i32, 0, n[29] as i32, 0);
            rct!(r[217], n[24] as i32, n[12] as i32, 0, n[13] as i32, n[7] as i32);
            rct!(r[218], n[11] as i32, n[18] as i32, 0, n[12] as i32, n[27] as i32);
            rct!(r[219], n[29] as i32, 0, 0, n[12] as i32, n[16] as i32);
            rct!(r[220], n[11] as i32, n[18] as i32, 0, n[12] as i32, n[16] as i32);
        }

        // Divide by species densities to get d/dX out of loss
        let ntot = s.ntot;
        let a = &mut s.a_mat;
        for icol in 0..ntot {
            let xc = s.xr[icol];
            if xc != 0.0 {
                for irow in 0..ntot {
                    if a[[irow, icol]] != 0.0 {
                        a[[irow, icol]] /= xc;
                    }
                }
            }
        }

        // Closure Jacobian rows (steady-state)
        if s.deltt < 1.0e-20 {
            let ntot = s.ntot;
            for icol in 0..ntot {
                s.a_mat[[idx(innoy), icol]] = 0.0;
            }
            s.a_mat[[idx(innoy), idx(n[0])]]  = -1.0;
            s.a_mat[[idx(innoy), idx(n[1])]]  = -1.0;
            s.a_mat[[idx(innoy), idx(n[2])]]  = -1.0;
            s.a_mat[[idx(innoy), idx(n[3])]]  = -2.0;
            s.a_mat[[idx(innoy), idx(n[4])]]  = -1.0;
            s.a_mat[[idx(innoy), idx(n[14])]] = -1.0;
            s.a_mat[[idx(innoy), idx(n[19])]] = -1.0;
            s.a_mat[[idx(innoy), idx(n[20])]] = -1.0;
            s.a_mat[[idx(innoy), idx(n[22])]] = -1.0;

            for icol in 0..ntot {
                s.a_mat[[idx(inclx), icol]] = 0.0;
            }
            s.a_mat[[idx(inclx), idx(n[15])]] = -1.0;
            s.a_mat[[idx(inclx), idx(n[16])]] = -1.0;
            s.a_mat[[idx(inclx), idx(n[17])]] = -2.0;
            s.a_mat[[idx(inclx), idx(n[18])]] = -1.0;
            s.a_mat[[idx(inclx), idx(n[19])]] = -1.0;
            s.a_mat[[idx(inclx), idx(n[21])]] = -1.0;
            s.a_mat[[idx(inclx), idx(n[27])]] = -1.0;
            s.a_mat[[idx(inclx), idx(n[28])]] = -2.0;

            if s.lbrom {
                for icol in 0..ntot {
                    s.a_mat[[idx(inbrx), icol]] = 0.0;
                }
                s.a_mat[[idx(inbrx), idx(n[11])]] = -1.0;
                s.a_mat[[idx(inbrx), idx(n[12])]] = -1.0;
                s.a_mat[[idx(inbrx), idx(n[13])]] = -1.0;
                s.a_mat[[idx(inbrx), idx(n[22])]] = -1.0;
                s.a_mat[[idx(inbrx), idx(n[23])]] = -1.0;
            }
        }

        DVector::zeros(ntot) // placeholder — caller uses A for Jacobian path
    }
}

// ── REACT — accumulate Jacobian entries ──────────────────────────────────────

/// Update A for one reaction. nl1, nl2, nl3 are reactants (1-based, may be negative);
/// np1, np2 are products (1-based, positive or zero).
/// Negative nl means: use |nl| as column key, but don't add self-loss for that row.
/// Fortran: SUBROUTINE REACT(RRATE, ML1, ML2, ML3, MP1, MP2)
fn react(a: &mut ndarray::Array2<f64>, rrate: f64, nl1: i32, nl2: i32, nl3: i32, np1: i32, np2: i32) {
    if nl1 == 0 { return; }
    let kol1 = (nl1.unsigned_abs() as usize).saturating_sub(1);
    if nl1 > 0 { a[[idx(nl1 as usize), kol1]] -= rrate; }
    if nl2 > 0 { a[[idx(nl2 as usize), kol1]] -= rrate; }
    if nl3 > 0 { a[[idx(nl3 as usize), kol1]] -= rrate; }
    if np1 > 0 { a[[idx(np1 as usize), kol1]] += rrate; }
    if np2 > 0 { a[[idx(np2 as usize), kol1]] += rrate; }

    if nl2 == 0 { return; }
    let kol2 = (nl2.unsigned_abs() as usize).saturating_sub(1);
    if nl1 > 0 { a[[idx(nl1 as usize), kol2]] -= rrate; }
    if nl2 > 0 { a[[idx(nl2 as usize), kol2]] -= rrate; }
    if nl3 > 0 { a[[idx(nl3 as usize), kol2]] -= rrate; }
    if np1 > 0 { a[[idx(np1 as usize), kol2]] += rrate; }
    if np2 > 0 { a[[idx(np2 as usize), kol2]] += rrate; }

    if nl3 == 0 { return; }
    let kol3 = (nl3.unsigned_abs() as usize).saturating_sub(1);
    if nl1 > 0 { a[[idx(nl1 as usize), kol3]] -= rrate; }
    if nl2 > 0 { a[[idx(nl2 as usize), kol3]] -= rrate; }
    if nl3 > 0 { a[[idx(nl3 as usize), kol3]] -= rrate; }
    if np1 > 0 { a[[idx(np1 as usize), kol3]] += rrate; }
    if np2 > 0 { a[[idx(np2 as usize), kol3]] += rrate; }
}
