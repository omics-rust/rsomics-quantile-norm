use std::process::Command;

fn ours() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_rsomics-quantile-norm"))
}

fn path(rel: &str) -> String {
    format!("{}/{}", env!("CARGO_MANIFEST_DIR"), rel)
}

fn cells(s: &str) -> Vec<Vec<String>> {
    s.trim()
        .lines()
        .map(|l| l.split('\t').map(str::to_string).collect())
        .collect()
}

fn run_ours(counts: &str) -> String {
    let out = Command::new(ours()).arg(counts).output().unwrap();
    assert!(
        out.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn diff_values(mine: &str, theirs: &str, eps: f64) {
    let m = cells(mine);
    let t = cells(theirs);
    assert_eq!(m[0], t[0], "header mismatch");
    assert_eq!(m.len(), t.len(), "row count mismatch");
    let mut max_dev = 0.0f64;
    for (r, (a, b)) in m[1..].iter().zip(t[1..].iter()).enumerate() {
        assert_eq!(a[0], b[0], "row {r} gene id mismatch");
        for (c, (av, bv)) in a[1..].iter().zip(b[1..].iter()).enumerate() {
            let x: f64 = av.parse().unwrap();
            let y: f64 = bv.parse().unwrap();
            let dev = (x - y).abs() / (1.0 + y.abs());
            max_dev = max_dev.max(dev);
            assert!(
                dev < eps,
                "row {r} col {c}: ours={x} oracle={y} reldev={dev}"
            );
        }
    }
    eprintln!("max relative deviation = {max_dev:e}");
}

// Always-on: diff ours against the committed golden (captured once from
// limma::normalizeQuantiles). Runs in CI with no R installed.
#[test]
fn matches_committed_golden() {
    let mine = run_ours(&path("tests/golden/counts.tsv"));
    let golden = std::fs::read_to_string(path("tests/golden/normalized.tsv")).unwrap();
    diff_values(&mine, &golden, 1e-6);
}

// Live differential vs limma in the staged conda env. Loud-skips when the
// r-bioc env is absent (CI), so the committed golden remains the CI gate.
#[test]
fn matches_limma_oracle() {
    if !r_bioc_available() {
        eprintln!("SKIP matches_limma_oracle: `conda run -n r-bioc Rscript` unavailable");
        return;
    }
    let counts = path("tests/golden/counts.tsv");
    let oracle = path("tests/quantile_norm_oracle.R");
    let ref_out = Command::new("conda")
        .args(["run", "-n", "r-bioc", "Rscript", &oracle, &counts])
        .output()
        .unwrap();
    assert!(
        ref_out.status.success(),
        "oracle failed: {}",
        String::from_utf8_lossy(&ref_out.stderr)
    );
    let theirs = String::from_utf8(ref_out.stdout).unwrap();
    let mine = run_ours(&counts);
    diff_values(&mine, &theirs, 1e-6);
}

fn r_bioc_available() -> bool {
    Command::new("conda")
        .args(["run", "-n", "r-bioc", "Rscript", "-e", "library(limma)"])
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}
