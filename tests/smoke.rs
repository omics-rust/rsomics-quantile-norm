use std::io::Write;
use std::process::Command;

fn ours() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_rsomics-quantile-norm"))
}

#[test]
fn normalizes_columns_to_a_common_distribution() {
    let dir = std::env::temp_dir().join("rsomics-qn-smoke");
    std::fs::create_dir_all(&dir).unwrap();
    let counts = dir.join("counts.tsv");
    let mut f = std::fs::File::create(&counts).unwrap();
    writeln!(f, "gene\ts1\ts2\ts3").unwrap();
    writeln!(f, "g1\t5\t4\t3").unwrap();
    writeln!(f, "g2\t2\t1\t4").unwrap();
    writeln!(f, "g3\t3\t4\t6").unwrap();
    writeln!(f, "g4\t4\t2\t8").unwrap();
    drop(f);

    let out = Command::new(ours()).arg(&counts).output().unwrap();
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let lines: Vec<&str> = s.trim().lines().collect();
    assert_eq!(lines[0], "gene\ts1\ts2\ts3");
    assert_eq!(lines.len(), 5);

    // After quantile normalization every column has the same multiset of
    // values (the sorted reference quantiles), so column sums are equal.
    let n = lines[0].split('\t').count() - 1;
    let mut sums = vec![0.0f64; n];
    for line in &lines[1..] {
        for (j, v) in line.split('\t').skip(1).enumerate() {
            sums[j] += v.parse::<f64>().unwrap();
        }
    }
    for j in 1..n {
        assert!(
            (sums[j] - sums[0]).abs() < 1e-3,
            "column sums differ: {sums:?}"
        );
    }
}
