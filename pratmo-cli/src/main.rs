use anyhow::Result;
use clap::Parser;
use pratmo_core::{
    ctm::ctmlfq,
    diurnal::diurn,
    path::tpath,
    reader::{FortranReader, ModelReader},
    state::ModelState,
};
use std::io::BufWriter;

/// PRATMO: Stratospheric Photochemical Box Model (Rust port of Fortran v6.0)
#[derive(Parser, Debug)]
#[command(version, about)]
struct Args {
    /// Input file directory (must contain fort10_cam06.x, fort11_jpl09.x, etc.)
    #[arg(short, long, default_value = ".")]
    input_dir: std::path::PathBuf,

    /// Output file path (default: fort01.x in input_dir)
    #[arg(short, long)]
    output: Option<std::path::PathBuf>,
}

fn main() -> Result<()> {
    let args = Args::parse();

    println!(" Photochemical Box-Model  version 6.0 (Prather) [Rust port]");

    let mut state = ModelState::new();
    state.cinpdir = args.input_dir.to_string_lossy().into_owned();

    let mut reader = FortranReader::new(&args.input_dir);
    reader.read_all(&mut state)?;

    // Open output files (Fortran units 7, 8, 9)
    let out_dir = args.output.as_deref()
        .and_then(|p| p.parent())
        .unwrap_or(&args.input_dir);
    if let Ok(f) = std::fs::File::create(out_dir.join("fort07.x")) {
        state.out_unit7 = Some(BufWriter::new(f));
    }
    if let Ok(f) = std::fs::File::create(out_dir.join("fort08.x")) {
        state.out_unit8 = Some(BufWriter::new(f));
    }
    if let Ok(f) = std::fs::File::create(out_dir.join("fort09.x")) {
        state.out_unit9 = Some(BufWriter::new(f));
    }

    // Mode dispatch: mirrors batmo.f logic on ND216
    // Fortran: IF(LPRTJV) GOTO 1 → if LPRTJV true, skip to end (J-value print only mode)
    if !state.lprtjv {
        if state.nd216 > 0 {
            ctmlfq(&mut state)?;
        } else if state.nd216 == 0 {
            diurn(&mut state)?;
            tpath(&mut state)?;
        } else {
            eprintln!("DERIVS mode (nd216={}) not yet implemented", state.nd216);
        }
    }

    println!(" RAXLOOP={} RADCOUNT={}", state.raxloop, state.radcount);
    Ok(())
}
