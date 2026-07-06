use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub struct Matrix {
    pub header: String,
    pub genes: Vec<String>,
    pub n_genes: usize,
    pub n_samples: usize,
    /// Column-major: column `j` occupies `data[j*n_genes .. (j+1)*n_genes]`.
    pub data: Vec<f64>,
}

pub fn quantile_normalize(input: &Path, output: &mut dyn Write) -> Result<u64> {
    let m = read_matrix(input)?;
    let normed = normalize_matrix(&m);
    write_matrix(&m, &normed, output)?;
    Ok(m.n_genes as u64)
}

fn read_matrix(path: &Path) -> Result<Matrix> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut reader = BufReader::new(file);

    let mut header = String::new();
    if reader.read_line(&mut header).map_err(RsomicsError::Io)? == 0 {
        return Err(RsomicsError::InvalidInput("empty count matrix".into()));
    }
    while header.ends_with('\n') || header.ends_with('\r') {
        header.pop();
    }
    let n_samples = header.split('\t').count() - 1;
    if n_samples < 2 {
        return Err(RsomicsError::InvalidInput(
            "need at least 2 sample columns to quantile-normalize".into(),
        ));
    }

    // Read row-major into per-sample columns directly.
    let mut genes: Vec<String> = Vec::new();
    let mut columns: Vec<Vec<f64>> = vec![Vec::new(); n_samples];

    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line).map_err(RsomicsError::Io)? == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(['\n', '\r']);
        if trimmed.is_empty() {
            continue;
        }
        let mut fields = trimmed.split('\t');
        let gene = fields
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("row with no gene id".into()))?;
        genes.push(gene.to_string());
        for (j, col) in columns.iter_mut().enumerate() {
            let v = fields.next().ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "gene {gene}: expected {n_samples} values, column {j} missing"
                ))
            })?;
            col.push(parse_f64(v, gene)?);
        }
    }

    let n_genes = genes.len();
    let mut data = Vec::with_capacity(n_genes * n_samples);
    for col in columns {
        data.extend_from_slice(&col);
    }

    Ok(Matrix {
        header,
        genes,
        n_genes,
        n_samples,
        data,
    })
}

fn parse_f64(s: &str, gene: &str) -> Result<f64> {
    // read.delim maps the literal "NA" to a missing value; limma treats it as
    // NA-aware. We carry missing as NaN through the normalization.
    if s == "NA" {
        return Ok(f64::NAN);
    }
    s.parse::<f64>().map_err(|_| {
        RsomicsError::InvalidInput(format!("gene {gene}: cannot parse value '{s}' as a number"))
    })
}

/// limma `normalizeQuantiles(A, ties=TRUE)`, NA-aware.
///
/// Each column's non-missing values are sorted and, if the column has any NA,
/// stretched by linear interpolation onto the full `n`-point grid; the
/// across-column mean of those curves is the reference `m`. Every non-missing
/// value is then replaced by `m` interpolated at its average-tie rank taken
/// over the column's non-missing entries (`(rank-1)/(nobs-1)`), so equal values
/// map to the mean of the reference quantiles they collectively occupy. NA
/// positions stay NA.
pub fn normalize_matrix(m: &Matrix) -> Vec<f64> {
    let n = m.n_genes;
    let s = m.n_samples;
    if n == 0 {
        return Vec::new();
    }

    let mut ref_curve = vec![0.0f64; n];
    let mut sorted: Vec<f64> = Vec::with_capacity(n);
    for j in 0..s {
        let col = &m.data[j * n..(j + 1) * n];
        sorted.clear();
        sorted.extend(col.iter().copied().filter(|v| !v.is_nan()));
        sorted.sort_unstable_by(f64::total_cmp);
        let nobs = sorted.len();
        if nobs == n {
            for (acc, &v) in ref_curve.iter_mut().zip(sorted.iter()) {
                *acc += v;
            }
        } else if nobs == 0 {
            for acc in ref_curve.iter_mut() {
                *acc += f64::NAN;
            }
        } else {
            let span = (nobs - 1) as f64 / (n - 1) as f64;
            for (p, acc) in ref_curve.iter_mut().enumerate() {
                *acc += interp_uniform(&sorted, p as f64 * span);
            }
        }
    }
    let inv_s = 1.0 / s as f64;
    for v in &mut ref_curve {
        *v *= inv_s;
    }

    let mut out = vec![0.0f64; n * s];
    let mut avg_rank = vec![0.0f64; n];
    let mut order: Vec<u32> = Vec::with_capacity(n);
    for j in 0..s {
        let col = &m.data[j * n..(j + 1) * n];
        let nobs = average_ranks(col, &mut order, &mut avg_rank);
        let denom = (nobs as f64) - 1.0;
        let scale = (n - 1) as f64;
        let dst = &mut out[j * n..(j + 1) * n];
        for i in 0..n {
            if col[i].is_nan() {
                dst[i] = f64::NAN;
            } else {
                let quantile = (avg_rank[i] - 1.0) / denom;
                dst[i] = interp_uniform(&ref_curve, quantile * scale);
            }
        }
    }
    out
}

/// R `rank(x, ties.method="average")` over the non-NA entries of `col`. Writes
/// the 1-based average rank into `avg_rank` for each non-missing position
/// (missing positions are left untouched) and returns the non-missing count.
fn average_ranks(col: &[f64], order: &mut Vec<u32>, avg_rank: &mut [f64]) -> usize {
    order.clear();
    order.extend((0..col.len() as u32).filter(|&i| !col[i as usize].is_nan()));
    order.sort_unstable_by(|&a, &b| col[a as usize].total_cmp(&col[b as usize]));
    let nobs = order.len();
    let mut i = 0;
    while i < nobs {
        let vi = col[order[i] as usize];
        let mut k = i + 1;
        while k < nobs && col[order[k] as usize] == vi {
            k += 1;
        }
        let mean_rank = (i + k + 1) as f64 * 0.5;
        for &idx in &order[i..k] {
            avg_rank[idx as usize] = mean_rank;
        }
        i = k;
    }
    nobs
}

/// Linear interpolation of `arr` at a 0-based fractional index, clamped at the
/// ends (R `approx` / numpy `interp` edge behaviour). `arr` must be non-empty.
fn interp_uniform(arr: &[f64], fpos: f64) -> f64 {
    if arr.len() == 1 || fpos <= 0.0 {
        return arr[0];
    }
    let last = (arr.len() - 1) as f64;
    if fpos >= last {
        return arr[arr.len() - 1];
    }
    let lo = fpos.floor();
    let lo_i = lo as usize;
    let frac = fpos - lo;
    arr[lo_i] + frac * (arr[lo_i + 1] - arr[lo_i])
}

fn write_matrix(m: &Matrix, normed: &[f64], output: &mut dyn Write) -> Result<()> {
    let n = m.n_genes;
    let s = m.n_samples;
    let mut out = BufWriter::new(output);
    out.write_all(m.header.as_bytes())
        .map_err(RsomicsError::Io)?;
    out.write_all(b"\n").map_err(RsomicsError::Io)?;

    let mut buf: Vec<u8> = Vec::with_capacity(24);
    for i in 0..n {
        out.write_all(m.genes[i].as_bytes())
            .map_err(RsomicsError::Io)?;
        for j in 0..s {
            out.write_all(b"\t").map_err(RsomicsError::Io)?;
            buf.clear();
            fmt6(&mut buf, normed[j * n + i]);
            out.write_all(&buf).map_err(RsomicsError::Io)?;
        }
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

fn fmt6(buf: &mut Vec<u8>, v: f64) {
    use std::io::Write as _;
    if v.is_nan() {
        buf.extend_from_slice(b"NA");
    } else {
        write!(buf, "{v:.6}").unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mat(cols: &[&[f64]]) -> Matrix {
        let s = cols.len();
        let n = cols[0].len();
        let mut data = Vec::new();
        for c in cols {
            data.extend_from_slice(c);
        }
        Matrix {
            header: String::new(),
            genes: (0..n).map(|i| format!("g{i}")).collect(),
            n_genes: n,
            n_samples: s,
            data,
        }
    }

    fn col(out: &[f64], j: usize, n: usize) -> Vec<f64> {
        out[j * n..(j + 1) * n].to_vec()
    }

    #[test]
    fn matches_limma_three_col_example() {
        // limma normalizeQuantiles on cols [5,2,3,4] [4,1,4,2] [3,4,6,8].
        let m = mat(&[
            &[5.0, 2.0, 3.0, 4.0],
            &[4.0, 1.0, 4.0, 2.0],
            &[3.0, 4.0, 6.0, 8.0],
        ]);
        let out = normalize_matrix(&m);
        let approx = |a: &[f64], b: &[f64]| {
            for (x, y) in a.iter().zip(b) {
                assert!((x - y).abs() < 1e-9, "{x} vs {y}");
            }
        };
        approx(&col(&out, 0, 4), &[5.666666667, 2.0, 3.0, 4.666666667]);
        approx(&col(&out, 1, 4), &[5.166666667, 2.0, 5.166666667, 3.0]);
        approx(&col(&out, 2, 4), &[2.0, 3.0, 4.666666667, 5.666666667]);
    }

    #[test]
    fn average_rank_handles_ties() {
        let col = [4.0, 1.0, 4.0, 2.0];
        let mut order = Vec::new();
        let mut r = vec![0.0f64; 4];
        let nobs = average_ranks(&col, &mut order, &mut r);
        assert_eq!(nobs, 4);
        assert_eq!(r, vec![3.5, 1.0, 3.5, 2.0]);
    }

    #[test]
    fn identical_columns_are_unchanged() {
        let m = mat(&[&[3.0, 1.0, 2.0], &[3.0, 1.0, 2.0]]);
        let out = normalize_matrix(&m);
        assert_eq!(col(&out, 0, 3), vec![3.0, 1.0, 2.0]);
        assert_eq!(col(&out, 1, 3), vec![3.0, 1.0, 2.0]);
    }

    #[test]
    fn matches_limma_na_aware_example() {
        // limma normalizeQuantiles on cols [NaN,3,5] [2,8,1]:
        // g1: NaN,3.0  g2: 2.0,6.5  g3: 6.5,2.0 (NA excluded from rank/reference).
        let nan = f64::NAN;
        let m = mat(&[&[nan, 3.0, 5.0], &[2.0, 8.0, 1.0]]);
        let out = normalize_matrix(&m);
        assert!(out[0].is_nan());
        assert_eq!(&col(&out, 0, 3)[1..], &[2.0, 6.5]);
        assert_eq!(col(&out, 1, 3), vec![3.0, 6.5, 2.0]);
    }
}
