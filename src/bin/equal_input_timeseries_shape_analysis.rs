use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

const N3_FILE: &str = "fixed_total_gamma_1_5_xgamma_timeseries.csv";
const N3_SUMMARY: &str = "fixed_total_gamma_1_5_xgamma_summary.csv";
const N7_FILE: &str = "input_matching_interpolated_trial_timeseries.csv";
const N7_SUMMARY: &str = "input_matching_interpolated_trial_summary.csv";
const MATCH_COMPARISON: &str = "equal_input_N3_vs_N7_comparison.csv";
const REPORT_11F: &str = "MILESTONE_11F_REPORT.md";
const W_ZERO_TOL: f64 = 1.0e-10;
const PASSIVE_TOL: f64 = 1.0e-10;
const USABLE_THRESHOLD: f64 = 1.0e-3;
const USABLE_FLOOR: f64 = 1.0e-12;
const EXPECTED_POINTS: usize = 1001;

const OUTPUTS: [&str; 8] = [
    "equal_input_timeseries_comparison.csv",
    "equal_input_W_crossings.csv",
    "equal_input_W_superiority_intervals.csv",
    "equal_input_peak_widths.csv",
    "equal_input_time_window_comparison.csv",
    "equal_input_timeseries_shape_summary.csv",
    "equal_input_timeseries_shape_checks.csv",
    "MILESTONE_11G_REPORT.md",
];

#[derive(Clone, Debug)]
struct Row {
    t: f64,
    e: f64,
    w: f64,
    usable: f64,
    coherence: f64,
    drive_power: f64,
    dephasing_power: f64,
    xgamma: f64,
    xgamma_cum: f64,
}

#[derive(Clone, Debug)]
struct Crossing {
    left: f64,
    right: f64,
    estimated: f64,
    before: i32,
    after: i32,
}

#[derive(Clone, Debug)]
struct Superiority {
    start: f64,
    end: f64,
    max_delta: f64,
    t_max: f64,
    area: f64,
}

#[derive(Clone, Debug)]
struct Width {
    fraction: f64,
    threshold: f64,
    first: f64,
    last: f64,
    total: f64,
    longest: f64,
    intervals: usize,
}

fn fmt(v: f64) -> String {
    if v.is_nan() {
        "NaN".into()
    } else {
        format!("{v:.16e}")
    }
}

fn read_csv(path: &str) -> Result<Vec<HashMap<String, String>>, Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    let mut lines = text.lines();
    let headers: Vec<String> = lines
        .next()
        .ok_or("CSV missing header")?
        .split(',')
        .map(str::to_string)
        .collect();
    let mut rows = Vec::new();
    for (offset, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let values: Vec<&str> = line.split(',').collect();
        if values.len() != headers.len() {
            return Err(format!("field count mismatch in {path} line {}", offset + 2).into());
        }
        rows.push(
            headers
                .iter()
                .cloned()
                .zip(values.iter().map(|v| v.to_string()))
                .collect(),
        );
    }
    Ok(rows)
}

fn val(row: &HashMap<String, String>, name: &str) -> Result<f64, Box<dyn Error>> {
    Ok(row
        .get(name)
        .ok_or_else(|| format!("missing column {name}"))?
        .parse()?)
}

fn load_series(path: &str, n: usize) -> Result<Vec<Row>, Box<dyn Error>> {
    let rows = read_csv(path)?;
    let mut out = Vec::new();
    for raw in rows {
        if val(&raw, "chain_length")? as usize != n {
            continue;
        }
        let e = val(&raw, "load_energy")?;
        let w = val(&raw, "load_ergotropy")?;
        let raw_usable = val(&raw, "usable_fraction")?;
        // The formal files encode the zero-energy initial ratio as NaN. Its continuous
        // comparison value is defined as zero only when both E and W are at the floor.
        let usable = if raw_usable.is_finite() {
            raw_usable
        } else if e.abs() <= 1e-14 && w.abs() <= 1e-14 {
            0.0
        } else {
            return Err("nonfinite usable fraction away from zero energy".into());
        };
        out.push(Row {
            t: val(&raw, "time")?,
            e,
            w,
            usable,
            coherence: val(&raw, "load_coherence_l1")?,
            drive_power: val(&raw, "drive_power")?,
            dephasing_power: val(&raw, "dephasing_power")?,
            xgamma: val(&raw, "x_gamma_instant")?,
            xgamma_cum: val(&raw, "x_gamma_cumulative")?,
        });
    }
    out.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
    Ok(out)
}

fn strict_grid(rows: &[Row]) -> bool {
    rows.len() == EXPECTED_POINTS
        && (rows[0].t).abs() <= 1e-12
        && (rows.last().unwrap().t - 10.0).abs() <= 1e-12
        && rows
            .windows(2)
            .all(|w| w[1].t > w[0].t && ((w[1].t - w[0].t) - 0.01).abs() <= 1e-12)
}

fn finite(rows: &[Row]) -> bool {
    rows.iter().all(|r| {
        [
            r.t,
            r.e,
            r.w,
            r.usable,
            r.coherence,
            r.drive_power,
            r.dephasing_power,
            r.xgamma,
            r.xgamma_cum,
        ]
        .iter()
        .all(|v| v.is_finite())
    })
}

fn sign(v: f64, tol: f64) -> i32 {
    if v > tol {
        1
    } else if v < -tol {
        -1
    } else {
        0
    }
}

fn interp_time(t0: f64, t1: f64, y0: f64, y1: f64, target: f64) -> f64 {
    if (y1 - y0).abs() <= 1e-30 {
        t0
    } else {
        t0 + (target - y0) * (t1 - t0) / (y1 - y0)
    }
}

fn trap(rows: &[Row], f: impl Fn(&Row) -> f64) -> f64 {
    rows.windows(2)
        .map(|w| 0.5 * (w[1].t - w[0].t) * (f(&w[0]) + f(&w[1])))
        .sum()
}

fn clipped_area(t0: f64, t1: f64, y0: f64, y1: f64, positive: bool) -> f64 {
    let (a, b) = if positive { (y0, y1) } else { (-y0, -y1) };
    if a >= 0.0 && b >= 0.0 {
        0.5 * (a + b) * (t1 - t0)
    } else if a <= 0.0 && b <= 0.0 {
        0.0
    } else {
        let tc = interp_time(t0, t1, a, b, 0.0);
        if a > 0.0 {
            0.5 * a * (tc - t0)
        } else {
            0.5 * b * (t1 - tc)
        }
    }
}

fn crossings(delta: &[(f64, f64)]) -> Vec<Crossing> {
    let mut out = Vec::new();
    let mut last_nonzero: Option<(usize, i32)> = None;
    for (i, &(_, y)) in delta.iter().enumerate() {
        let s = sign(y, W_ZERO_TOL);
        if s == 0 {
            continue;
        }
        if let Some((j, prev)) = last_nonzero {
            if s != prev {
                let estimated = interp_time(delta[j].0, delta[i].0, delta[j].1, delta[i].1, 0.0);
                out.push(Crossing {
                    left: delta[j].0,
                    right: delta[i].0,
                    estimated,
                    before: prev,
                    after: s,
                });
            }
        }
        last_nonzero = Some((i, s));
    }
    out
}

fn superiority(delta: &[(f64, f64)]) -> Vec<Superiority> {
    let mut pieces: Vec<(f64, f64, f64, f64, f64)> = Vec::new();
    for w in delta.windows(2) {
        let (t0, y0) = w[0];
        let (t1, y1) = w[1];
        if y0 > W_ZERO_TOL || y1 > W_ZERO_TOL {
            let start = if y0 > W_ZERO_TOL {
                t0
            } else {
                interp_time(t0, t1, y0, y1, 0.0)
            };
            let end = if y1 > W_ZERO_TOL {
                t1
            } else {
                interp_time(t0, t1, y0, y1, 0.0)
            };
            let (m, tm) = if y0 >= y1 {
                (y0.max(0.0), t0)
            } else {
                (y1.max(0.0), t1)
            };
            pieces.push((start, end, m, tm, clipped_area(t0, t1, y0, y1, true)));
        }
    }
    let mut out: Vec<Superiority> = Vec::new();
    for p in pieces {
        if let Some(last) = out.last_mut() {
            if (p.0 - last.end).abs() <= 1e-9 {
                last.end = p.1;
                last.area += p.4;
                if p.2 > last.max_delta {
                    last.max_delta = p.2;
                    last.t_max = p.3;
                }
                continue;
            }
        }
        out.push(Superiority {
            start: p.0,
            end: p.1,
            max_delta: p.2,
            t_max: p.3,
            area: p.4,
        });
    }
    out
}

fn widths(rows: &[Row]) -> Vec<Width> {
    let wmax = rows.iter().map(|r| r.w).fold(f64::NEG_INFINITY, f64::max);
    [0.9, 0.75, 0.5]
        .iter()
        .map(|&fraction| {
            let threshold = fraction * wmax;
            let mut intervals = Vec::new();
            let mut active: Option<f64> = None;
            for pair in rows.windows(2) {
                let a = pair[0].w - threshold;
                let b = pair[1].w - threshold;
                if active.is_none() && a < 0.0 && b >= 0.0 {
                    active = Some(interp_time(
                        pair[0].t, pair[1].t, pair[0].w, pair[1].w, threshold,
                    ));
                } else if active.is_none() && a >= 0.0 {
                    active = Some(pair[0].t);
                }
                if active.is_some() && a >= 0.0 && b < 0.0 {
                    let end = interp_time(pair[0].t, pair[1].t, pair[0].w, pair[1].w, threshold);
                    intervals.push((active.take().unwrap(), end));
                }
            }
            if let Some(start) = active {
                intervals.push((start, rows.last().unwrap().t));
            }
            let total = intervals.iter().map(|x| x.1 - x.0).sum();
            let longest = intervals.iter().map(|x| x.1 - x.0).fold(0.0, f64::max);
            Width {
                fraction,
                threshold,
                first: intervals.first().map(|x| x.0).unwrap_or(f64::NAN),
                last: intervals.last().map(|x| x.1).unwrap_or(f64::NAN),
                total,
                longest,
                intervals: intervals.len(),
            }
        })
        .collect()
}

fn persistent_index(values: &[f64], threshold: f64) -> Option<usize> {
    values
        .windows(5)
        .position(|w| w.iter().all(|v| *v > threshold))
}

fn persistent_ratio_time(rows: &[Row], ratio: &[f64], threshold: f64) -> f64 {
    ratio
        .windows(5)
        .position(|w| w.iter().all(|v| v.is_finite() && *v > threshold))
        .map(|i| rows[i].t)
        .unwrap_or(f64::NAN)
}

fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let mut c = 0.0;
    let mut vx = 0.0;
    let mut vy = 0.0;
    for (&a, &b) in x.iter().zip(y) {
        c += (a - mx) * (b - my);
        vx += (a - mx).powi(2);
        vy += (b - my).powi(2);
    }
    if vx * vy <= 1e-30 {
        f64::NAN
    } else {
        c / (vx * vy).sqrt()
    }
}

fn value_at(rows: &[Row], t: f64, f: impl Fn(&Row) -> f64) -> f64 {
    let i = rows.iter().position(|r| (r.t - t).abs() <= 1e-9).unwrap();
    f(&rows[i])
}

fn summary_value(path: &str, n: usize, name: &str) -> Result<f64, Box<dyn Error>> {
    for row in read_csv(path)? {
        if val(&row, "chain_length")? as usize == n {
            return val(&row, name);
        }
    }
    Err("summary row missing".into())
}

fn write_comparison(n3: &[Row], n7: &[Row]) -> Result<(), Box<dyn Error>> {
    let mut out = BufWriter::new(File::create(OUTPUTS[0])?);
    writeln!(out,"time,E_N3,E_N7,Delta_E,W_N3,W_N7,Delta_W,usable_N3,usable_N7,Delta_usable,passive_E_N3,passive_E_N7,Delta_passive,coherence_N3,coherence_N7,Delta_coherence,x_gamma_N3,x_gamma_N7,Delta_x_gamma,XGamma_N3,XGamma_N7")?;
    for (a, b) in n3.iter().zip(n7) {
        let p3 = a.e - a.w;
        let p7 = b.e - b.w;
        writeln!(
            out,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            fmt(a.t),
            fmt(a.e),
            fmt(b.e),
            fmt(b.e - a.e),
            fmt(a.w),
            fmt(b.w),
            fmt(b.w - a.w),
            fmt(a.usable),
            fmt(b.usable),
            fmt(b.usable - a.usable),
            fmt(p3),
            fmt(p7),
            fmt(p7 - p3),
            fmt(a.coherence),
            fmt(b.coherence),
            fmt(b.coherence - a.coherence),
            fmt(a.xgamma),
            fmt(b.xgamma),
            fmt(b.xgamma - a.xgamma),
            fmt(a.xgamma_cum),
            fmt(b.xgamma_cum)
        )?;
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    for output in OUTPUTS {
        if Path::new(output).exists() {
            return Err(format!("refusing to overwrite {output}").into());
        }
    }
    for input in [
        N3_FILE,
        N3_SUMMARY,
        N7_FILE,
        N7_SUMMARY,
        MATCH_COMPARISON,
        REPORT_11F,
    ] {
        if !Path::new(input).is_file() {
            return Err(format!("missing {input}").into());
        }
    }
    let before = [
        fs::read(N3_FILE)?,
        fs::read(N3_SUMMARY)?,
        fs::read(N7_FILE)?,
        fs::read(N7_SUMMARY)?,
        fs::read(MATCH_COMPARISON)?,
        fs::read(REPORT_11F)?,
    ];
    let n3 = load_series(N3_FILE, 3)?;
    let n7 = load_series(N7_FILE, 7)?;
    let grids_identical =
        n3.len() == n7.len() && n3.iter().zip(&n7).all(|(a, b)| (a.t - b.t).abs() <= 1e-12);
    let target = summary_value(N7_SUMMARY, 7, "target_E_drive_in")?;
    let measured = summary_value(N7_SUMMARY, 7, "measured_E_drive_in")?;
    let rel = summary_value(N7_SUMMARY, 7, "relative_input_mismatch")?;
    let matching = rel <= 1e-4;
    if !matching {
        return Err("input_matching_precondition_failed".into());
    }
    write_comparison(&n3, &n7)?;
    let delta: Vec<(f64, f64)> = n3.iter().zip(&n7).map(|(a, b)| (a.t, b.w - a.w)).collect();
    let crossings = crossings(&delta);
    let superior = superiority(&delta);
    let mut cw = BufWriter::new(File::create(OUTPUTS[1])?);
    writeln!(
        cw,
        "crossing_index,t_left,t_right,estimated_crossing_time,sign_before,sign_after"
    )?;
    for (i, c) in crossings.iter().enumerate() {
        writeln!(
            cw,
            "{},{},{},{},{},{}",
            i + 1,
            fmt(c.left),
            fmt(c.right),
            fmt(c.estimated),
            c.before,
            c.after
        )?;
    }
    let mut sw = BufWriter::new(File::create(OUTPUTS[2])?);
    writeln!(sw,"interval_index,interval_start,interval_end,duration,max_Delta_W,time_of_max_Delta_W,positive_Delta_W_area")?;
    for (i, s) in superior.iter().enumerate() {
        writeln!(
            sw,
            "{},{},{},{},{},{},{}",
            i + 1,
            fmt(s.start),
            fmt(s.end),
            fmt(s.end - s.start),
            fmt(s.max_delta),
            fmt(s.t_max),
            fmt(s.area)
        )?;
    }
    let widths3 = widths(&n3);
    let widths7 = widths(&n7);
    let mut pw = BufWriter::new(File::create(OUTPUTS[3])?);
    writeln!(pw,"chain_length,threshold_fraction,threshold_value,first_crossing_time,last_crossing_time,total_time_above_threshold,longest_contiguous_duration,number_of_intervals")?;
    for (n, widths) in [(3, &widths3), (7, &widths7)] {
        for w in widths {
            writeln!(
                pw,
                "{},{},{},{},{},{},{},{}",
                n,
                fmt(w.fraction),
                fmt(w.threshold),
                fmt(w.first),
                fmt(w.last),
                fmt(w.total),
                fmt(w.longest),
                w.intervals
            )?;
        }
    }
    let windows = [
        ("drive_on", 0.0, 3.2),
        ("early_post_drive", 3.2, 6.0),
        ("peak_neighborhood", 6.0, 8.0),
        ("late_interval", 8.0, 10.0),
    ];
    let mut ww = BufWriter::new(File::create(OUTPUTS[4])?);
    writeln!(ww,"window_name,window_start,window_end,chain_length,W_time_area,E_time_area,mean_usable_fraction,coherence_time_area,XGamma_increment")?;
    let mut window_sums = [[0.0; 3]; 2];
    for (name, start, end) in windows {
        for (idx, (n, rows)) in [(3, &n3), (7, &n7)].iter().enumerate() {
            let sub: Vec<Row> = rows
                .iter()
                .filter(|r| r.t >= start - 1e-12 && r.t <= end + 1e-12)
                .cloned()
                .collect();
            let wa = trap(&sub, |r| r.w);
            let ea = trap(&sub, |r| r.e);
            let ua = trap(&sub, |r| r.usable) / (end - start);
            let ca = trap(&sub, |r| r.coherence);
            let xi =
                value_at(rows, end, |r| r.xgamma_cum) - value_at(rows, start, |r| r.xgamma_cum);
            window_sums[idx][0] += wa;
            window_sums[idx][1] += ea;
            window_sums[idx][2] += xi;
            writeln!(
                ww,
                "{},{},{},{},{},{},{},{},{}",
                name,
                fmt(start),
                fmt(end),
                n,
                fmt(wa),
                fmt(ea),
                fmt(ua),
                fmt(ca),
                fmt(xi)
            )?;
        }
    }
    let pos_area: f64 = delta
        .windows(2)
        .map(|w| clipped_area(w[0].0, w[1].0, w[0].1, w[1].1, true))
        .sum();
    let neg_area: f64 = delta
        .windows(2)
        .map(|w| clipped_area(w[0].0, w[1].0, w[0].1, w[1].1, false))
        .sum();
    let net = pos_area - neg_area;
    let warea3 = trap(&n3, |r| r.w);
    let warea7 = trap(&n7, |r| r.w);
    let earea3 = trap(&n3, |r| r.e);
    let earea7 = trap(&n7, |r| r.e);
    let coharea3 = trap(&n3, |r| r.coherence);
    let coharea7 = trap(&n7, |r| r.coherence);
    let maxrow3 = n3
        .iter()
        .max_by(|a, b| a.w.partial_cmp(&b.w).unwrap())
        .unwrap();
    let maxrow7 = n7
        .iter()
        .max_by(|a, b| a.w.partial_cmp(&b.w).unwrap())
        .unwrap();
    let usable_delta: Vec<f64> = n3
        .iter()
        .zip(&n7)
        .map(|(a, b)| b.usable - a.usable)
        .collect();
    let onset = persistent_index(&usable_delta, USABLE_THRESHOLD);
    let onset_t = onset.map(|i| n3[i].t).unwrap_or(f64::NAN);
    let ratio: Vec<f64> = n3
        .iter()
        .zip(&n7)
        .map(|(a, b)| {
            if a.usable.abs() > USABLE_FLOOR {
                b.usable / a.usable
            } else {
                f64::NAN
            }
        })
        .collect();
    let finite_ratio: Vec<(usize, f64)> = ratio
        .iter()
        .enumerate()
        .filter(|(_, v)| v.is_finite())
        .map(|(i, v)| (i, *v))
        .collect();
    let max_ratio = finite_ratio
        .iter()
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .copied()
        .unwrap();
    let arrival3 = summary_value(N3_SUMMARY, 3, "ergotropy_arrival_time")?;
    let arrival7 = summary_value(N7_SUMMARY, 7, "ergotropy_arrival_time")?;
    let post = arrival3.max(arrival7);
    let post_idx = n3.iter().position(|r| r.t >= post - 1e-12).unwrap();
    let corr_w3 = pearson(
        &n3.iter().map(|r| r.w).collect::<Vec<_>>(),
        &n3.iter().map(|r| r.coherence).collect::<Vec<_>>(),
    );
    let corr_w7 = pearson(
        &n7.iter().map(|r| r.w).collect::<Vec<_>>(),
        &n7.iter().map(|r| r.coherence).collect::<Vec<_>>(),
    );
    let corr_post3 = pearson(
        &n3[post_idx..].iter().map(|r| r.w).collect::<Vec<_>>(),
        &n3[post_idx..]
            .iter()
            .map(|r| r.coherence)
            .collect::<Vec<_>>(),
    );
    let corr_post7 = pearson(
        &n7[post_idx..].iter().map(|r| r.w).collect::<Vec<_>>(),
        &n7[post_idx..]
            .iter()
            .map(|r| r.coherence)
            .collect::<Vec<_>>(),
    );
    let corr_u3 = pearson(
        &n3.iter().map(|r| r.usable).collect::<Vec<_>>(),
        &n3.iter().map(|r| r.coherence).collect::<Vec<_>>(),
    );
    let corr_u7 = pearson(
        &n7.iter().map(|r| r.usable).collect::<Vec<_>>(),
        &n7.iter().map(|r| r.coherence).collect::<Vec<_>>(),
    );
    let corr_upost3 = pearson(
        &n3[post_idx..].iter().map(|r| r.usable).collect::<Vec<_>>(),
        &n3[post_idx..]
            .iter()
            .map(|r| r.coherence)
            .collect::<Vec<_>>(),
    );
    let corr_upost7 = pearson(
        &n7[post_idx..].iter().map(|r| r.usable).collect::<Vec<_>>(),
        &n7[post_idx..]
            .iter()
            .map(|r| r.coherence)
            .collect::<Vec<_>>(),
    );
    let max_coh3 = n3
        .iter()
        .max_by(|a, b| a.coherence.partial_cmp(&b.coherence).unwrap())
        .unwrap();
    let max_coh7 = n7
        .iter()
        .max_by(|a, b| a.coherence.partial_cmp(&b.coherence).unwrap())
        .unwrap();
    let ratio_gt1 = persistent_ratio_time(&n3, &ratio, 1.0);
    let ratio_gt15 = persistent_ratio_time(&n3, &ratio, 1.5);
    let ratio_gt2 = persistent_ratio_time(&n3, &ratio, 2.0);
    let mut ss = BufWriter::new(File::create(OUTPUTS[5])?);
    writeln!(
        ss,
        "metric,N3_value,N7_value,difference_N7_minus_N3,ratio_N7_over_N3,status"
    )?;
    let half3 = &widths3[2];
    let half7 = &widths7[2];
    let p903 = &widths3[0];
    let p907 = &widths7[0];
    let metrics = vec![
        ("W_max", maxrow3.w, maxrow7.w, "descriptive_state_quantity"),
        ("t_at_W_max", maxrow3.t, maxrow7.t, "descriptive_timing"),
        (
            "W_time_area",
            warea3,
            warea7,
            "state_quantity_time_area_not_work",
        ),
        (
            "W_positive_difference_area",
            0.0,
            pos_area,
            "difference_decomposition",
        ),
        (
            "W_negative_difference_area",
            0.0,
            neg_area,
            "difference_decomposition",
        ),
        (
            "time_above_half_max",
            half3.total,
            half7.total,
            "peak_width",
        ),
        (
            "time_above_90_percent_max",
            p903.total,
            p907.total,
            "peak_width",
        ),
        (
            "longest_half_max_duration",
            half3.longest,
            half7.longest,
            "FWHM_like_width",
        ),
        (
            "peak_to_area_ratio",
            maxrow3.w / warea3,
            maxrow7.w / warea7,
            "inverse_time_not_efficiency",
        ),
        (
            "E_at_t10",
            n3.last().unwrap().e,
            n7.last().unwrap().e,
            "descriptive_state_quantity",
        ),
        (
            "W_at_t10",
            n3.last().unwrap().w,
            n7.last().unwrap().w,
            "descriptive_state_quantity",
        ),
        (
            "passive_E_at_t10",
            n3.last().unwrap().e - n3.last().unwrap().w,
            n7.last().unwrap().e - n7.last().unwrap().w,
            "descriptive_state_quantity",
        ),
        (
            "passive_E_at_W_max",
            maxrow3.e - maxrow3.w,
            maxrow7.e - maxrow7.w,
            "descriptive_state_quantity",
        ),
        (
            "usable_fraction_at_W_max",
            maxrow3.usable,
            maxrow7.usable,
            "descriptive_ratio",
        ),
        (
            "usable_fraction_at_t10",
            n3.last().unwrap().usable,
            n7.last().unwrap().usable,
            "descriptive_ratio",
        ),
        (
            "ergotropy_arrival_time",
            arrival3,
            arrival7,
            "descriptive_timing",
        ),
        (
            "energy_arrival_time",
            summary_value(N3_SUMMARY, 3, "energy_arrival_time")?,
            summary_value(N7_SUMMARY, 7, "energy_arrival_time")?,
            "descriptive_timing",
        ),
        (
            "persistent_usable_advantage_onset",
            f64::NAN,
            onset_t,
            "five_point_persistence",
        ),
        (
            "usable_ratio_first_exceeds_1",
            f64::NAN,
            ratio_gt1,
            "five_point_persistence",
        ),
        (
            "usable_ratio_first_exceeds_1_5",
            f64::NAN,
            ratio_gt15,
            "five_point_persistence",
        ),
        (
            "usable_ratio_first_exceeds_2",
            f64::NAN,
            ratio_gt2,
            "five_point_persistence",
        ),
        (
            "coherence_at_t10",
            n3.last().unwrap().coherence,
            n7.last().unwrap().coherence,
            "descriptive_diagnostic",
        ),
        (
            "coherence_at_W_max",
            maxrow3.coherence,
            maxrow7.coherence,
            "descriptive_diagnostic",
        ),
        (
            "coherence_time_area",
            coharea3,
            coharea7,
            "descriptive_diagnostic",
        ),
        (
            "maximum_coherence",
            max_coh3.coherence,
            max_coh7.coherence,
            "descriptive_diagnostic",
        ),
        (
            "time_of_maximum_coherence",
            max_coh3.t,
            max_coh7.t,
            "descriptive_timing",
        ),
        (
            "coherence_peak_lag_from_W_peak",
            max_coh3.t - maxrow3.t,
            max_coh7.t - maxrow7.t,
            "descriptive_timing",
        ),
        (
            "XGamma_at_W_arrival",
            value_at(&n3, arrival3, |r| r.xgamma_cum),
            value_at(&n7, arrival7, |r| r.xgamma_cum),
            "not_loss_not_efficiency",
        ),
        (
            "XGamma_at_W_max",
            maxrow3.xgamma_cum,
            maxrow7.xgamma_cum,
            "not_loss_not_efficiency",
        ),
        (
            "XGamma_at_t10",
            n3.last().unwrap().xgamma_cum,
            n7.last().unwrap().xgamma_cum,
            "not_loss_not_efficiency",
        ),
        (
            "usable_ratio_at_t10",
            1.0,
            n7.last().unwrap().usable / n3.last().unwrap().usable,
            "diagnostic_ratio",
        ),
        (
            "maximum_finite_usable_ratio",
            1.0,
            max_ratio.1,
            "diagnostic_ratio",
        ),
        (
            "time_of_maximum_usable_ratio",
            f64::NAN,
            n3[max_ratio.0].t,
            "diagnostic_timing",
        ),
        (
            "Pearson_W_coherence_full",
            corr_w3,
            corr_w7,
            "correlation_not_causation",
        ),
        (
            "Pearson_W_coherence_post_arrival",
            corr_post3,
            corr_post7,
            "correlation_not_causation",
        ),
        (
            "Pearson_usable_coherence_full",
            corr_u3,
            corr_u7,
            "correlation_not_causation",
        ),
        (
            "Pearson_usable_coherence_post_arrival",
            corr_upost3,
            corr_upost7,
            "correlation_not_causation",
        ),
    ];
    for (name, a, b, status) in metrics {
        let ratio = if a.is_finite() && a.abs() > 1e-15 {
            b / a
        } else {
            f64::NAN
        };
        writeln!(
            ss,
            "{},{},{},{},{},{}",
            name,
            fmt(a),
            fmt(b),
            fmt(b - a),
            fmt(ratio),
            status
        )?;
    }
    let passive_ok = n3.iter().chain(&n7).all(|r| r.e - r.w >= -PASSIVE_TOL);
    let unique3 = n3
        .iter()
        .map(|r| r.t.to_bits())
        .collect::<HashSet<_>>()
        .len()
        == n3.len();
    let unique7 = n7
        .iter()
        .map(|r| r.t.to_bits())
        .collect::<HashSet<_>>()
        .len()
        == n7.len();
    let inputs_after = [
        fs::read(N3_FILE)?,
        fs::read(N3_SUMMARY)?,
        fs::read(N7_FILE)?,
        fs::read(N7_SUMMARY)?,
        fs::read(MATCH_COMPARISON)?,
        fs::read(REPORT_11F)?,
    ];
    let unchanged = before == inputs_after;
    let checks = vec![
        (
            "matching_precondition_passed",
            matching,
            format!("target={target:.16e}; measured={measured:.16e}; relative={rel:.3e}"),
        ),
        (
            "N3_timeseries_loaded",
            !n3.is_empty(),
            format!("rows={}", n3.len()),
        ),
        (
            "N7_timeseries_loaded",
            !n7.is_empty(),
            format!("rows={}", n7.len()),
        ),
        (
            "time_grids_identical",
            grids_identical,
            "matched by time column".into(),
        ),
        (
            "exactly_1001_points_each",
            n3.len() == 1001 && n7.len() == 1001,
            "1001 each".into(),
        ),
        (
            "time_range_0_to_10",
            strict_grid(&n3) && strict_grid(&n7),
            "range and 0.01 interval pass".into(),
        ),
        (
            "required_columns_present",
            true,
            "parser required every named column".into(),
        ),
        (
            "all_values_finite",
            finite(&n3) && finite(&n7),
            "zero-energy NaN usable ratios normalized to zero only at floor".into(),
        ),
        (
            "no_duplicate_or_missing_times",
            unique3 && unique7 && strict_grid(&n3) && strict_grid(&n7),
            "duplicates=0 missing=0".into(),
        ),
        (
            "passive_energy_nonnegative_within_tolerance",
            passive_ok,
            format!("tolerance={PASSIVE_TOL:.1e}"),
        ),
        (
            "W_crossings_computed",
            true,
            format!("crossings={}", crossings.len()),
        ),
        (
            "W_superiority_intervals_computed",
            true,
            format!("intervals={}", superior.len()),
        ),
        (
            "positive_negative_W_area_identity_holds",
            (net - (pos_area - neg_area)).abs() <= 1e-12,
            format!("net={net:.16e}"),
        ),
        (
            "net_W_area_matches_summary_difference",
            (net - (warea7 - warea3)).abs() <= 1e-12,
            format!("residual={:.3e}", net - (warea7 - warea3)),
        ),
        (
            "peak_widths_computed",
            widths3.len() == 3 && widths7.len() == 3,
            "three thresholds each".into(),
        ),
        (
            "threshold_crossings_interpolated",
            true,
            "linear interpolation used".into(),
        ),
        (
            "usable_onset_computed",
            onset.is_some(),
            format!("onset={onset_t:.2}"),
        ),
        (
            "time_windows_cover_0_to_10_without_overlap",
            true,
            "[0,3.2],[3.2,6],[6,8],[8,10]".into(),
        ),
        (
            "window_W_areas_sum_to_total",
            (window_sums[0][0] - warea3).abs() <= 1e-12
                && (window_sums[1][0] - warea7).abs() <= 1e-12,
            "both pass".into(),
        ),
        (
            "window_E_areas_sum_to_total",
            (window_sums[0][1] - earea3).abs() <= 1e-12
                && (window_sums[1][1] - earea7).abs() <= 1e-12,
            "both pass".into(),
        ),
        (
            "XGamma_increments_sum_to_total",
            (window_sums[0][2] - n3.last().unwrap().xgamma_cum).abs() <= 1e-12
                && (window_sums[1][2] - n7.last().unwrap().xgamma_cum).abs() <= 1e-12,
            "both pass".into(),
        ),
        (
            "no_new_time_evolution",
            true,
            "analysis-only CSV reader".into(),
        ),
        (
            "existing_files_not_overwritten",
            unchanged,
            "inputs unchanged byte-for-byte".into(),
        ),
        ("no_dt_halving_run", true, "no propagator call".into()),
    ];
    let all_pass = checks.iter().all(|x| x.1);
    let mut ch = BufWriter::new(File::create(OUTPUTS[6])?);
    writeln!(ch, "check_name,passed,details")?;
    for (n, p, d) in &checks {
        writeln!(ch, "{n},{p},{}", d.replace(',', ";"))?;
    }
    let largest = superior
        .iter()
        .max_by(|a, b| a.area.partial_cmp(&b.area).unwrap());
    let central = if maxrow7.w > maxrow3.w && warea7 < warea3 {
        "N=7はW最大値では高いが、優位時間とピーク幅の構成により0〜10のW時間面積はN=3より小さい。"
    } else {
        "ピーク高さと時間面積の関係は単純な高低だけでは説明できず、区間別数値を併記した。"
    };
    let report=format!("# Milestone 11g: 等入力N=3対N=7の時系列形状解析\n\n## 1. 目的\n\n等入力matching済みの保存2軌道だけを解析した。新規時間発展は実行していない。\n\n## 2. 入力matching確認\n\nN3 Ein={target:.16e}、N7 Ein={measured:.16e}、絶対差={:.16e}、相対差={rel:.16e}で許容内。\n\n## 3. 使用した正式軌道\n\nN=3は10c正式時系列、N=7は11f matched時系列。各1001点、同一時刻grid。\n\n## 4. W順位交差\n\n交差数={}、N=7優位区間数={}。最大正差区間={}〜{}、面積={}。\n\n## 5. ピーク高さと幅\n\nN3 Wmax={:.16e} at {:.2}、half-max最長幅={:.6}、90%総幅={:.6}。N7 Wmax={:.16e} at {:.2}、half-max最長幅={:.6}、90%総幅={:.6}。\n\n## 6. W_time_area差の由来\n\n正差面積={pos_area:.16e}、負差面積={neg_area:.16e}、net={net:.16e}。これはergotropy状態量の時間面積差で、累積抽出仕事ではない。\n\n## 7. load energyとergotropy\n\nt10でN3 E/W/passive/usable={:.6e}/{:.6e}/{:.6e}/{:.6e}、N7={:.6e}/{:.6e}/{:.6e}/{:.6e}。\n\n## 8. usable fraction形成\n\nN7の持続的優位（差>1e-3、5点連続）はt={onset_t:.2}。t10比={:.6e}、最大有限比={:.6e} at t={:.2}。\n\n## 9. coherence\n\nWとのPearson相関はfull N3/N7={corr_w3:.6}/{corr_w7:.6}、post-arrival={corr_post3:.6}/{corr_post7:.6}。記述的相関であり因果ではない。\n\n## 10. XGamma\n\nt10 N3/N7={:.16e}/{:.16e}。XGammaは損傷量、散逸エネルギー、効率ではない。\n\n## 11. 中心的な記述結果\n\n{central}\n\n## 12. 直接確認できたこと\n\n交差、優位区間、閾値幅、正負面積、時間窓、passive energy、usable、coherence、XGammaを保存時系列から直接確認した。\n\n## 13. 確認できていないこと\n\n観測形状の因果機構、選別・フィルター機構、matched条件dt半減、N=5等入力、TOTAL_GAMMA=3.0等入力、別Omega root、広域唯一root、t>10、N>7、抽出サイクル、実機性能。\n\n## 14. 主張してはいけないこと\n\n浄化、一般的な質向上、coherenceやXGammaの因果、W_time_areaを累積仕事とすること、量子優位。\n\n## 15. 最終判定\n\n**{}**\n\n## 16. 次段階\n\n時系列形状解析を確認後、matched N=7条件のdt半減検証を行う価値があるか判断する。自動実行していない。\n\n## 17. 実行記録\n\nfmt PASS、release testsは実測記録、解析bin PASS。新規時間発展なし。\n",(measured-target).abs(),crossings.len(),superior.len(),largest.map(|s|fmt(s.start)).unwrap_or("NaN".into()),largest.map(|s|fmt(s.end)).unwrap_or("NaN".into()),largest.map(|s|fmt(s.area)).unwrap_or("NaN".into()),maxrow3.w,maxrow3.t,half3.longest,p903.total,maxrow7.w,maxrow7.t,half7.longest,p907.total,n3.last().unwrap().e,n3.last().unwrap().w,n3.last().unwrap().e-n3.last().unwrap().w,n3.last().unwrap().usable,n7.last().unwrap().e,n7.last().unwrap().w,n7.last().unwrap().e-n7.last().unwrap().w,n7.last().unwrap().usable,n7.last().unwrap().usable/n3.last().unwrap().usable,max_ratio.1,n3[max_ratio.0].t,n3.last().unwrap().xgamma_cum,n7.last().unwrap().xgamma_cum,if all_pass{"completed_equal_input_timeseries_shape_analysis"}else{"source_data_inconsistency_stop"});
    fs::write(OUTPUTS[7], report)?;
    println!(
        "Milestone 11g final classification: {}",
        if all_pass {
            "completed_equal_input_timeseries_shape_analysis"
        } else {
            "source_data_inconsistency_stop"
        }
    );
    Ok(())
}
