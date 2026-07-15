use std::{
    path::{Path, PathBuf},
    sync::OnceLock,
    time::{SystemTime, UNIX_EPOCH},
};

/// Copy the tracked legacy inputs to an isolated writable directory.
///
/// The CTM driver writes `boxout.dat` beside its inputs, so tests must not run
/// directly against the source fixture tree. Keeping the copy behind a
/// `OnceLock` also makes parallel tests share one completed model run safely.
pub fn input_dir() -> &'static Path {
    static INPUT_DIR: OnceLock<PathBuf> = OnceLock::new();

    INPUT_DIR
        .get_or_init(|| {
            let source = Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests")
                .join("fixtures")
                .join("legacy_inputs");
            let nonce = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .expect("system clock predates Unix epoch")
                .as_nanos();
            let destination = std::env::temp_dir().join(format!(
                "pratmo-legacy-inputs-{}-{nonce}",
                std::process::id()
            ));
            std::fs::create_dir_all(&destination).expect("cannot create test input directory");

            for entry in std::fs::read_dir(&source).expect("cannot read tracked test fixtures") {
                let entry = entry.expect("cannot read tracked test fixture entry");
                if entry
                    .file_type()
                    .expect("cannot stat test fixture")
                    .is_file()
                {
                    std::fs::copy(entry.path(), destination.join(entry.file_name()))
                        .expect("cannot copy tracked test fixture");
                }
            }

            destination
        })
        .as_path()
}
