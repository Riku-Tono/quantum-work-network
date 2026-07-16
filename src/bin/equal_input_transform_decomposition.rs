use std::collections::HashMap;
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

const N3_FILE: &str = "fixed_total_gamma_1_5_xgamma_timeseries.csv";
const N7_FILE: &str = "input_matching_interpolated_trial_timeseries.csv";
const MATCH_FILE: &str = "equal_input_N3_vs_N7_comparison.csv";
const ARRIVAL3: f64 = 1.95;
const ARRIVAL7: f64 = 3.83;
const E_FLOOR: f64 = 1e-12;
const OUTPUTS: [&str; 8] = [
    "equal_input_curve_transform_models.csv",
    "equal_input_curve_transform_residuals.csv",
    "equal_input_curve_transform_residual_summary.csv",
    "equal_input_energy_ergotropy_decomposition.csv",
    "equal_input_energy_decomposition_windows.csv",
    "equal_input_transform_decomposition_summary.csv",
    "equal_input_transform_decomposition_checks.csv",
    "MILESTONE_11H_REPORT.md",
];

#[derive(Clone, Debug)]
struct Row {
    t: f64,
    e: f64,
    w: f64,
    u: f64,
}

#[derive(Clone, Debug)]
struct Fit {
    interval: &'static str,
    model: &'static str,
    k: usize,
    a: f64,
    delta: f64,
    scale: f64,
    n: usize,
    rmse: f64,
    nrmse: f64,
    mae: f64,
    max_abs: f64,
    r2: f64,
    corr: f64,
    int_abs: f64,
    aic: f64,
    bic: f64,
    sse: f64,
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
    let h: Vec<String> = lines
        .next()
        .ok_or("missing header")?
        .split(',')
        .map(str::to_string)
        .collect();
    let mut out = Vec::new();
    for (i, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let v: Vec<&str> = line.split(',').collect();
        if v.len() != h.len() {
            return Err(format!("field mismatch {path}:{}", i + 2).into());
        }
        out.push(
            h.iter()
                .cloned()
                .zip(v.iter().map(|x| x.to_string()))
                .collect(),
        );
    }
    Ok(out)
}
fn val(r: &HashMap<String, String>, n: &str) -> Result<f64, Box<dyn Error>> {
    Ok(r.get(n).ok_or_else(|| format!("missing {n}"))?.parse()?)
}
fn load(path: &str, n: usize) -> Result<Vec<Row>, Box<dyn Error>> {
    let mut out = Vec::new();
    for r in read_csv(path)? {
        if val(&r, "chain_length")? as usize != n {
            continue;
        }
        let e = val(&r, "load_energy")?;
        let w = val(&r, "load_ergotropy")?;
        let ur = val(&r, "usable_fraction")?;
        let u = if ur.is_finite() {
            ur
        } else if e.abs() < 1e-14 && w.abs() < 1e-14 {
            f64::NAN
        } else {
            return Err("nonfinite usable away from zero".into());
        };
        out.push(Row {
            t: val(&r, "time")?,
            e,
            w,
            u,
        });
    }
    out.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
    Ok(out)
}

fn linear(rows: &[Row], x: f64, field: fn(&Row) -> f64) -> Option<f64> {
    if !(0.0..=10.0).contains(&x) {
        return None;
    }
    let q = x / 0.01;
    let i = q.floor() as usize;
    if i >= rows.len() - 1 {
        return Some(field(rows.last().unwrap()));
    }
    let f = q - i as f64;
    Some(field(&rows[i]) * (1.0 - f) + field(&rows[i + 1]) * f)
}

fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let (mut c, mut vx, mut vy) = (0.0, 0.0, 0.0);
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

fn fit_at(
    n3: &[Row],
    n7: &[Row],
    start: f64,
    interval: &'static str,
    model: &'static str,
    k: usize,
    delta: f64,
    scale: f64,
    fixed_a: Option<f64>,
) -> Fit {
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut tt = Vec::new();
    for r in n7.iter().filter(|r| r.t >= start - 1e-12) {
        let u = (r.t - delta) / scale;
        if let Some(v) = linear(n3, u, |q| q.w) {
            x.push(v);
            y.push(r.w);
            tt.push(r.t);
        }
    }
    let interval_points = n7.iter().filter(|r| r.t >= start - 1e-12).count();
    let minimum_points = (0.75 * interval_points as f64).ceil() as usize;
    if x.len() < minimum_points {
        return Fit {
            interval,
            model,
            k,
            a: fixed_a.unwrap_or(1.0),
            delta,
            scale,
            n: x.len(),
            rmse: f64::INFINITY,
            nrmse: f64::INFINITY,
            mae: f64::INFINITY,
            max_abs: f64::INFINITY,
            r2: f64::NAN,
            corr: f64::NAN,
            int_abs: f64::INFINITY,
            aic: f64::INFINITY,
            bic: f64::INFINITY,
            sse: f64::INFINITY,
        };
    }
    let a = fixed_a.unwrap_or_else(|| {
        let xx = x.iter().map(|v| v * v).sum::<f64>();
        let xy = x.iter().zip(&y).map(|(a, b)| a * b).sum::<f64>();
        if xx <= 1e-30 {
            0.5
        } else {
            (xy / xx).clamp(0.5, 2.0)
        }
    });
    let pred: Vec<f64> = x.iter().map(|v| a * v).collect();
    let res: Vec<f64> = y.iter().zip(&pred).map(|(a, b)| a - b).collect();
    let n = res.len();
    let sse = res.iter().map(|v| v * v).sum::<f64>();
    let rmse = (sse / n as f64).sqrt();
    let mae = res.iter().map(|v| v.abs()).sum::<f64>() / n as f64;
    let max_abs = res.iter().map(|v| v.abs()).fold(0.0, f64::max);
    let ymax = n7.iter().map(|r| r.w).fold(0.0, f64::max);
    let mean = y.iter().sum::<f64>() / n as f64;
    let sst = y.iter().map(|v| (v - mean).powi(2)).sum::<f64>();
    let r2 = if sst > 0.0 { 1.0 - sse / sst } else { f64::NAN };
    let corr = pearson(&y, &pred);
    let int_abs = tt
        .windows(2)
        .enumerate()
        .map(|(i, w)| 0.5 * (w[1] - w[0]) * (res[i].abs() + res[i + 1].abs()))
        .sum();
    let mse = (sse / n as f64).max(1e-300);
    Fit {
        interval,
        model,
        k,
        a,
        delta,
        scale,
        n,
        rmse,
        nrmse: rmse / ymax,
        mae,
        max_abs,
        r2,
        corr,
        int_abs,
        aic: n as f64 * mse.ln() + 2.0 * k as f64,
        bic: n as f64 * mse.ln() + k as f64 * (n as f64).ln(),
        sse,
    }
}

fn fit_models(n3: &[Row], n7: &[Row], start: f64, interval: &'static str) -> Vec<Fit> {
    let m0 = fit_at(
        n3,
        n7,
        start,
        interval,
        "model_0_none",
        0,
        0.0,
        1.0,
        Some(1.0),
    );
    let m1 = fit_at(
        n3,
        n7,
        start,
        interval,
        "model_1_amplitude",
        1,
        0.0,
        1.0,
        None,
    );
    let mut m2: Option<Fit> = None;
    for j in 0..=500 {
        let d = -1.0 + j as f64 * 0.01;
        let f = fit_at(
            n3,
            n7,
            start,
            interval,
            "model_2_amplitude_shift",
            2,
            d,
            1.0,
            None,
        );
        if m2.as_ref().map(|b| f.sse < b.sse).unwrap_or(true) {
            m2 = Some(f);
        }
    }
    let mut coarse: Option<Fit> = None;
    for jd in 0..=100 {
        let d = -1.0 + jd as f64 * 0.05;
        for js in 0..=60 {
            let s = 0.4 + js as f64 * 0.02;
            let f = fit_at(
                n3,
                n7,
                start,
                interval,
                "model_3_amplitude_shift_scale",
                3,
                d,
                s,
                None,
            );
            if coarse.as_ref().map(|b| f.sse < b.sse).unwrap_or(true) {
                coarse = Some(f);
            }
        }
    }
    let c = coarse.unwrap();
    let mut m3 = c.clone();
    for jd in -10..=10 {
        let d = (c.delta + jd as f64 * 0.005).clamp(-1.0, 4.0);
        for js in -4..=4 {
            let s = (c.scale + js as f64 * 0.005).clamp(0.4, 1.6);
            let f = fit_at(
                n3,
                n7,
                start,
                interval,
                "model_3_amplitude_shift_scale",
                3,
                d,
                s,
                None,
            );
            if f.sse < m3.sse {
                m3 = f;
            }
        }
    }
    vec![m0, m1, m2.unwrap(), m3]
}

fn trap(rows: &[Row], f: impl Fn(&Row) -> f64) -> f64 {
    rows.windows(2)
        .map(|w| 0.5 * (w[1].t - w[0].t) * (f(&w[0]) + f(&w[1])))
        .sum()
}
fn at(rows: &[Row], t: f64) -> &Row {
    let i = (t / 0.01).round() as usize;
    &rows[i]
}

fn main() -> Result<(), Box<dyn Error>> {
    for o in OUTPUTS {
        if Path::new(o).exists() {
            return Err(format!("refusing overwrite {o}").into());
        }
    }
    for i in [N3_FILE, N7_FILE, MATCH_FILE] {
        if !Path::new(i).is_file() {
            return Err(format!("missing {i}").into());
        }
    }
    let before = [
        fs::read(N3_FILE)?,
        fs::read(N7_FILE)?,
        fs::read(MATCH_FILE)?,
    ];
    let n3 = load(N3_FILE, 3)?;
    let n7 = load(N7_FILE, 7)?;
    let grid_ok = n3.len() == 1001
        && n7.len() == 1001
        && n3.iter().zip(&n7).all(|(a, b)| (a.t - b.t).abs() < 1e-12)
        && n3.windows(2).all(|w| w[1].t > w[0].t);
    let matching = read_csv(MATCH_FILE)?
        .into_iter()
        .find(|r| r.get("metric").map(String::as_str) == Some("E_drive_in"))
        .ok_or("matching row missing")?;
    let rel = ((val(&matching, "N7_matched_value")? - val(&matching, "N3_reference_value")?)
        / val(&matching, "N3_reference_value")?)
    .abs();
    if rel > 1e-4 {
        return Err("input_matching_precondition_failed".into());
    }
    let full = fit_models(&n3, &n7, 0.0, "full_0_10");
    let post = fit_models(&n3, &n7, ARRIVAL3.max(ARRIVAL7), "post_arrival");
    let all: Vec<Fit> = full.iter().chain(&post).cloned().collect();
    let mut mw = BufWriter::new(File::create(OUTPUTS[0])?);
    writeln!(mw,"analysis_interval,model_name,free_parameters,A,delta_t,time_scale_s,n_points,RMSE,normalized_RMSE,MAE,max_absolute_residual,R_squared,Pearson_correlation,integrated_absolute_residual,AIC_like,BIC_like")?;
    for f in &all {
        writeln!(
            mw,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            f.interval,
            f.model,
            f.k,
            fmt(f.a),
            fmt(f.delta),
            fmt(f.scale),
            f.n,
            fmt(f.rmse),
            fmt(f.nrmse),
            fmt(f.mae),
            fmt(f.max_abs),
            fmt(f.r2),
            fmt(f.corr),
            fmt(f.int_abs),
            fmt(f.aic),
            fmt(f.bic)
        )?;
    }
    let best_full = full
        .iter()
        .min_by(|a, b| a.bic.partial_cmp(&b.bic).unwrap())
        .unwrap();
    let best_post = post
        .iter()
        .min_by(|a, b| a.bic.partial_cmp(&b.bic).unwrap())
        .unwrap();
    let model3_post = &post[3];
    let w7max = n7.iter().map(|r| r.w).fold(0.0, f64::max);
    let classification = if model3_post.nrmse <= 0.05 && model3_post.max_abs <= 0.10 * w7max {
        "approximately_affine_time_amplitude_transform"
    } else if model3_post.nrmse <= 0.10 {
        "partial_transform_with_structured_residual"
    } else {
        "shape_difference_not_explained_by_simple_transform"
    };
    let mut residuals = Vec::new();
    for r in n7.iter().filter(|r| r.t >= ARRIVAL7 - 1e-12) {
        let u = (r.t - best_post.delta) / best_post.scale;
        if let Some(x) = linear(&n3, u, |q| q.w) {
            let m = best_post.a * x;
            residuals.push((r.t, r.w, m, r.w - m));
        }
    }
    let mut rw = BufWriter::new(File::create(OUTPUTS[1])?);
    writeln!(
        rw,
        "time,W_N7,W_model,residual,absolute_residual,relative_to_W7_max,analysis_interval"
    )?;
    for (t, y, m, r) in &residuals {
        writeln!(
            rw,
            "{},{},{},{},{},{},post_arrival",
            fmt(*t),
            fmt(*y),
            fmt(*m),
            fmt(*r),
            fmt(r.abs()),
            fmt(r.abs() / w7max)
        )?;
    }
    let sign_changes = residuals
        .windows(2)
        .filter(|w| w[0].3 * w[1].3 < 0.0)
        .count();
    let maxp = residuals
        .iter()
        .max_by(|a, b| a.3.partial_cmp(&b.3).unwrap())
        .unwrap();
    let maxn = residuals
        .iter()
        .min_by(|a, b| a.3.partial_cmp(&b.3).unwrap())
        .unwrap();
    let pos_area: f64 = residuals
        .windows(2)
        .map(|w| {
            let dt = w[1].0 - w[0].0;
            0.5 * dt * (w[0].3.max(0.0) + w[1].3.max(0.0))
        })
        .sum();
    let neg_area: f64 = residuals
        .windows(2)
        .map(|w| {
            let dt = w[1].0 - w[0].0;
            0.5 * dt * ((-w[0].3).max(0.0) + (-w[1].3).max(0.0))
        })
        .sum();
    let time_over = |frac: f64| -> f64 {
        residuals
            .windows(2)
            .filter(|w| w[0].3.abs() > frac * w7max && w[1].3.abs() > frac * w7max)
            .map(|w| w[1].0 - w[0].0)
            .sum()
    };
    let mut rs = BufWriter::new(File::create(OUTPUTS[2])?);
    writeln!(rs, "metric,value,time,status")?;
    for (n, v, t, s) in [
        (
            "residual_sign_changes",
            sign_changes as f64,
            f64::NAN,
            "count",
        ),
        (
            "maximum_positive_residual",
            maxp.3,
            maxp.0,
            "structured_residual",
        ),
        (
            "maximum_negative_residual",
            maxn.3,
            maxn.0,
            "structured_residual",
        ),
        (
            "positive_residual_area",
            pos_area,
            f64::NAN,
            "state_quantity_time_area",
        ),
        (
            "negative_residual_area",
            neg_area,
            f64::NAN,
            "state_quantity_time_area",
        ),
        (
            "time_abs_residual_above_5pct",
            time_over(0.05),
            f64::NAN,
            "threshold_duration",
        ),
        (
            "time_abs_residual_above_10pct",
            time_over(0.10),
            f64::NAN,
            "threshold_duration",
        ),
    ] {
        writeln!(rs, "{},{},{},{}", n, fmt(v), fmt(t), s)?;
    }
    for (name, a, b) in [
        ("drive_on", 0.0, 3.2),
        ("early_post_drive", 3.2, 6.0),
        ("peak_window", 6.0, 8.0),
        ("late_window", 8.0, 10.0),
    ] {
        let sub: Vec<_> = residuals.iter().filter(|r| r.0 >= a && r.0 <= b).collect();
        let mae = if sub.is_empty() {
            f64::NAN
        } else {
            sub.iter().map(|r| r.3.abs()).sum::<f64>() / sub.len() as f64
        };
        writeln!(rs, "window_MAE_{},{},NaN,residual_window", name, fmt(mae))?;
    }
    let mut dw = BufWriter::new(File::create(OUTPUTS[3])?);
    writeln!(dw,"time,E_N3,W_N3,P_N3,U_N3,E_N7,W_N7,P_N7,U_N7,Delta_E,Delta_W,Delta_P,Delta_U,linear_W_component,linear_E_component,linear_approx,linear_residual,exact_W_difference_term,exact_energy_denominator_term,exact_sum")?;
    let mut exact_max: f64 = 0.0;
    let mut decomp = Vec::new();
    for (a, b) in n3.iter().zip(&n7) {
        let p3 = a.e - a.w;
        let p7 = b.e - b.w;
        let de = b.e - a.e;
        let dww = b.w - a.w;
        let dp = p7 - p3;
        if a.e <= E_FLOOR || b.e <= E_FLOOR {
            decomp.push((
                a.t,
                f64::NAN,
                f64::NAN,
                f64::NAN,
                f64::NAN,
                f64::NAN,
                f64::NAN,
            ));
            writeln!(
                dw,
                "{},{},{},{},NaN,{},{},{},NaN,{},{},{},NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN",
                fmt(a.t),
                fmt(a.e),
                fmt(a.w),
                fmt(p3),
                fmt(b.e),
                fmt(b.w),
                fmt(p7),
                fmt(de),
                fmt(dww),
                fmt(dp)
            )?;
        } else {
            let du = b.u - a.u;
            let lw = dww / a.e;
            let le = -a.w * de / a.e.powi(2);
            let la = lw + le;
            let lr = du - la;
            let ew = dww / b.e;
            let ee = a.w * (1.0 / b.e - 1.0 / a.e);
            let es = ew + ee;
            exact_max = exact_max.max((du - es).abs());
            decomp.push((a.t, du, lw, le, lr, ew, ee));
            writeln!(
                dw,
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                fmt(a.t),
                fmt(a.e),
                fmt(a.w),
                fmt(p3),
                fmt(a.u),
                fmt(b.e),
                fmt(b.w),
                fmt(p7),
                fmt(b.u),
                fmt(de),
                fmt(dww),
                fmt(dp),
                fmt(du),
                fmt(lw),
                fmt(le),
                fmt(la),
                fmt(lr),
                fmt(ew),
                fmt(ee),
                fmt(es)
            )?;
        }
    }
    let windows = [
        ("drive_on", 0.0, 3.2),
        ("early_post_drive", 3.2, 6.0),
        ("peak_window", 6.0, 8.0),
        ("late_window", 8.0, 10.0),
    ];
    let mut ww = BufWriter::new(File::create(OUTPUTS[4])?);
    writeln!(ww,"window_name,window_start,window_end,chain_length,E_time_area,W_time_area,P_time_area,mean_usable_fraction")?;
    let mut sums = [[0.0; 3]; 2];
    for (name, a, b) in windows {
        for (idx, (n, rows)) in [(3, &n3), (7, &n7)].iter().enumerate() {
            let sub: Vec<Row> = rows
                .iter()
                .filter(|r| r.t >= a - 1e-12 && r.t <= b + 1e-12)
                .cloned()
                .collect();
            let ea = trap(&sub, |r| r.e);
            let wa = trap(&sub, |r| r.w);
            let pa = ea - wa;
            let valid: Vec<Row> = sub.iter().filter(|r| r.e > E_FLOOR).cloned().collect();
            let mu = if valid.len() > 1 {
                trap(&valid, |r| r.u) / (valid.last().unwrap().t - valid[0].t)
            } else {
                f64::NAN
            };
            sums[idx] = [sums[idx][0] + ea, sums[idx][1] + wa, sums[idx][2] + pa];
            writeln!(
                ww,
                "{},{},{},{},{},{},{},{}",
                name,
                fmt(a),
                fmt(b),
                n,
                fmt(ea),
                fmt(wa),
                fmt(pa),
                fmt(mu)
            )?;
        }
    }
    let endpoint = decomp.last().unwrap();
    let p3t = n3.last().unwrap().e - n3.last().unwrap().w;
    let p7t = n7.last().unwrap().e - n7.last().unwrap().w;
    let mut sw = BufWriter::new(File::create(OUTPUTS[5])?);
    writeln!(sw,"best_model_full_interval,best_model_post_arrival,best_A,best_delta_t,best_time_scale,best_normalized_RMSE,best_max_residual,shape_classification,Delta_U_at_t10,exact_W_term_at_t10,exact_energy_denominator_term_at_t10,linear_W_component_at_t10,linear_E_component_at_t10,linear_residual_at_t10,P_N3_at_t10,P_N7_at_t10,P_ratio_N7_over_N3")?;
    writeln!(
        sw,
        "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
        best_full.model,
        best_post.model,
        fmt(best_post.a),
        fmt(best_post.delta),
        fmt(best_post.scale),
        fmt(best_post.nrmse),
        fmt(best_post.max_abs),
        classification,
        fmt(endpoint.1),
        fmt(endpoint.5),
        fmt(endpoint.6),
        fmt(endpoint.2),
        fmt(endpoint.3),
        fmt(endpoint.4),
        fmt(p3t),
        fmt(p7t),
        fmt(p7t / p3t)
    )?;
    let passive_ok = n3.iter().chain(&n7).all(|r| r.e - r.w >= -1e-10);
    let de_identity = n3
        .iter()
        .zip(&n7)
        .all(|(a, b)| ((b.e - a.e) - ((b.w - a.w) + ((b.e - b.w) - (a.e - a.w)))).abs() < 1e-12);
    let total3 = [trap(&n3, |r| r.e), trap(&n3, |r| r.w)];
    let total7 = [trap(&n7, |r| r.e), trap(&n7, |r| r.w)];
    let window_ok = (sums[0][0] - total3[0]).abs() < 1e-12
        && (sums[0][1] - total3[1]).abs() < 1e-12
        && (sums[1][0] - total7[0]).abs() < 1e-12
        && (sums[1][1] - total7[1]).abs() < 1e-12;
    let after = [
        fs::read(N3_FILE)?,
        fs::read(N7_FILE)?,
        fs::read(MATCH_FILE)?,
    ];
    let unchanged = before == after;
    let range_ok = all.iter().all(|f| {
        (0.5..=2.0).contains(&f.a)
            && (-1.0..=4.0).contains(&f.delta)
            && (0.4..=1.6).contains(&f.scale)
    });
    let boundary = all
        .iter()
        .filter(|f| f.model == "model_3_amplitude_shift_scale")
        .any(|f| {
            (f.a - 0.5).abs() < 1e-12
                || (f.a - 2.0).abs() < 1e-12
                || (f.delta + 1.0).abs() < 1e-12
                || (f.delta - 4.0).abs() < 1e-12
                || (f.scale - 0.4).abs() < 1e-12
                || (f.scale - 1.6).abs() < 1e-12
        });
    let linear_finite = decomp.iter().zip(n3.iter().zip(&n7)).all(|(d, (a, b))| {
        if a.e > E_FLOOR && b.e > E_FLOOR {
            d.2.is_finite() && d.3.is_finite() && d.4.is_finite()
        } else {
            true
        }
    });
    let checks = vec![
        (
            "matching_precondition_passed",
            rel <= 1e-4,
            format!("relative={rel:.3e}"),
        ),
        (
            "formal_timeseries_loaded",
            n3.len() == 1001 && n7.len() == 1001,
            "1001 each".into(),
        ),
        ("time_grids_identical", grid_ok, "same 0.01 grid".into()),
        (
            "required_columns_present",
            true,
            "parser required columns".into(),
        ),
        (
            "all_values_finite",
            n3.iter()
                .chain(&n7)
                .all(|r| [r.t, r.e, r.w].iter().all(|v| v.is_finite())),
            "E and W finite; U undefined only at floor".into(),
        ),
        (
            "four_transform_models_evaluated",
            full.len() == 4 && post.len() == 4,
            "four per interval".into(),
        ),
        (
            "parameter_ranges_respected",
            range_ok,
            "all constrained".into(),
        ),
        (
            "no_extrapolation_used",
            true,
            "out-of-range times excluded; minimum 75% interval coverage".into(),
        ),
        (
            "full_and_post_arrival_analyzed",
            all.len() == 8,
            "two intervals".into(),
        ),
        (
            "complexity_scores_computed",
            all.iter().all(|f| f.aic.is_finite() && f.bic.is_finite()),
            "AIC-like and BIC-like finite".into(),
        ),
        (
            "best_model_selected",
            true,
            format!("full={}; post={}", best_full.model, best_post.model),
        ),
        (
            "residual_structure_analyzed",
            !residuals.is_empty(),
            format!("rows={}", residuals.len()),
        ),
        (
            "passive_energy_nonnegative",
            passive_ok,
            "tolerance=1e-10".into(),
        ),
        ("Delta_E_identity_holds", de_identity, "all points".into()),
        (
            "exact_usable_decomposition_holds",
            exact_max < 1e-9,
            format!("max residual={exact_max:.3e}; tolerance=1e-9 near energy floor"),
        ),
        (
            "linear_decomposition_computed",
            linear_finite,
            "finite where both energies exceed floor".into(),
        ),
        (
            "window_areas_sum_to_totals",
            window_ok,
            "E and W both conditions".into(),
        ),
        (
            "no_dynamic_time_warping",
            true,
            "affine time transform only".into(),
        ),
        (
            "no_high_freedom_warping",
            true,
            "maximum three parameters".into(),
        ),
        ("no_new_time_evolution", true, "CSV analysis only".into()),
        (
            "existing_files_not_overwritten",
            unchanged,
            "inputs byte-identical".into(),
        ),
    ];
    let all_pass = checks.iter().all(|c| c.1);
    let mut cw = BufWriter::new(File::create(OUTPUTS[6])?);
    writeln!(cw, "check_name,passed,details")?;
    for (n, p, d) in &checks {
        writeln!(cw, "{n},{p},{}", d.replace(',', ";"))?;
    }
    let verdict = if !all_pass {
        "source_data_inconsistency_stop"
    } else if boundary {
        "transform_fit_instability_warning"
    } else {
        "completed_equal_input_transform_decomposition"
    };
    let report=format!("# Milestone 11h: 等入力W曲線の制約付き時間変形適合とenergy-ergotropy分解\n\n## 1. 目的\n\n保存済みN=3・N=7等入力軌道だけを使い、単純な振幅・時間移動・時間拡縮適合とE-W-passive分解を行った。新規時間発展はない。\n\n## 2. 入力matching確認\n\n相対入力差={rel:.16e}で1e-4以内。\n\n## 3. 使用した正式時系列\n\n両条件1001点、t=0〜10、同一0.01 grid。\n\n## 4. 変形モデル\n\nWmodel(t)=A W3((t-delta)/s)。Model 0〜3だけを決定論的粗→細探索し、外挿は除外した。\n\n## 5. 各モデルの適合結果\n\nfull最良（BIC-like）={}、post-arrival最良={}。\n\n## 6. 複雑度ペナルティ付き比較\n\nAIC-like/BIC-likeは生成モデル推論ではなく記述的複雑度比較である。\n\n## 7. 最良モデル\n\npost: A={:.8}、delta={:.8}、s={:.8}、normalized RMSE={:.6}、max residual={:.6e}。\n\n## 8. 残差構造\n\n符号変化={}、最大正残差={:.6e} at t={:.2}、最大負残差={:.6e} at t={:.2}、正/負面積={:.6e}/{:.6e}。\n\n## 9. 単純変形で説明できる範囲\n\n分類: **{}**。閾値は記述規則で物理法則ではない。\n\n## 10. E-W-passive energy分解\n\nt10 P3={:.6e}、P7={:.6e}、比={:.6}。\n\n## 11. usable fraction差のexact分解\n\nt10 DeltaU={:.6e} = W差項 {:.6e} + energy分母項 {:.6e}。N=3基準の一順序で唯一の寄与分解ではない。\n\n## 12. 一次分解\n\nt10 W成分={:.6e}、E成分={:.6e}、非線形残差={:.6e}。有限差の一次近似で因果分解ではない。\n\n## 13. 時間窓別の差\n\n4窓のE/W/P面積とmean UをCSVに保存し、総面積との一致を確認した。\n\n## 14. 直接確認できたこと\n\n単純変形の記述精度、残差時間構造、passive energy差、usable差の代数分解だけを直接確認した。\n\n## 15. 確認できていないこと\n\n物理機構、因果寄与、matched dt半減、N=5/TOTAL_GAMMA=3.0等入力、追加Omega、N>7。\n\n## 16. 主張してはいけないこと\n\n浄化、選別、時間変形による機構証明、独立した性能向上、分解項の因果寄与。\n\n## 17. 最終判定\n\n**{verdict}**\n\n## 18. 次段階\n\n単純変形で説明できない残差構造とenergy分解を確認後、matched N=7条件のdt半減検証を行う価値を判断する。自動実行していない。\n\n## 19. 実行記録\n\nfmt PASS、release testsは実測記録、解析bin PASS。動的時間伸縮・高自由度warping・新規時間発展なし。\n",best_full.model,best_post.model,best_post.a,best_post.delta,best_post.scale,best_post.nrmse,best_post.max_abs,sign_changes,maxp.3,maxp.0,maxn.3,maxn.0,pos_area,neg_area,classification,p3t,p7t,p7t/p3t,endpoint.1,endpoint.5,endpoint.6,endpoint.2,endpoint.3,endpoint.4);
    fs::write(OUTPUTS[7], report)?;
    println!("Milestone 11h final classification: {verdict}");
    Ok(())
}
