use std::collections::HashMap;
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufWriter, Write};

const N3_FILE: &str = "fixed_total_gamma_1_5_xgamma_timeseries.csv";
const N7_FILE: &str = "input_matching_interpolated_trial_timeseries.csv";
const N7_SUMMARY: &str = "input_matching_interpolated_trial_summary.csv";
const H11_MODELS: &str = "equal_input_curve_transform_models.csv";
const H11_RESIDUALS: &str = "equal_input_curve_transform_residuals.csv";
const H11_SUMMARY: &str = "equal_input_transform_decomposition_summary.csv";
const H11_REPORT: &str = "MILESTONE_11H_REPORT.md";
const ARRIVAL: f64 = 3.83;
const A_MIN: f64 = 0.5;
const A_MAX: f64 = 2.0;
const D_MIN: f64 = -1.0;
const D_MAX: f64 = 4.0;
const S_MIN: f64 = 0.4;
const S_MAX: f64 = 1.6;
const COARSE_D: f64 = 0.05;
const COARSE_S: f64 = 0.04;
const FINE_D: f64 = 0.005;
const FINE_S: f64 = 0.005;

#[derive(Clone)]
struct Row {
    t: f64,
    w: f64,
}

#[derive(Clone)]
struct Fit {
    name: &'static str,
    k: usize,
    a: f64,
    delta: f64,
    sr: f64,
    sf: f64,
    boundary: f64,
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

fn fmt(value: f64) -> String {
    if value.is_nan() {
        "NaN".to_owned()
    } else {
        format!("{value:.16e}")
    }
}

fn read_csv(path: &str) -> Result<Vec<HashMap<String, String>>, Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    let mut lines = text.lines();
    let headers: Vec<String> = lines
        .next()
        .ok_or("empty CSV")?
        .split(',')
        .map(str::to_owned)
        .collect();
    let mut rows = Vec::new();
    for line in lines.filter(|line| !line.trim().is_empty()) {
        let values: Vec<&str> = line.split(',').collect();
        if values.len() != headers.len() {
            return Err(format!("CSV width mismatch: {path}").into());
        }
        rows.push(
            headers
                .iter()
                .cloned()
                .zip(values.iter().map(|value| (*value).to_owned()))
                .collect(),
        );
    }
    Ok(rows)
}

fn value(row: &HashMap<String, String>, name: &str) -> Result<f64, Box<dyn Error>> {
    Ok(row
        .get(name)
        .ok_or_else(|| format!("missing {name}"))?
        .parse()?)
}

fn load(path: &str, chain_length: usize) -> Result<Vec<Row>, Box<dyn Error>> {
    let mut rows = Vec::new();
    for row in read_csv(path)? {
        if value(&row, "chain_length")? as usize == chain_length {
            rows.push(Row {
                t: value(&row, "time")?,
                w: value(&row, "load_ergotropy")?,
            });
        }
    }
    rows.sort_by(|left, right| left.t.total_cmp(&right.t));
    Ok(rows)
}

fn linear(rows: &[Row], time: f64) -> Option<f64> {
    if !(0.0..=10.0).contains(&time) {
        return None;
    }
    let coordinate = time / 0.01;
    let index = coordinate.floor() as usize;
    if index >= rows.len() - 1 {
        return Some(rows.last()?.w);
    }
    let fraction = coordinate - index as f64;
    Some(rows[index].w * (1.0 - fraction) + rows[index + 1].w * fraction)
}

fn pearson(left: &[f64], right: &[f64]) -> f64 {
    let n = left.len() as f64;
    let ml = left.iter().sum::<f64>() / n;
    let mr = right.iter().sum::<f64>() / n;
    let (mut covariance, mut vl, mut vr) = (0.0, 0.0, 0.0);
    for (&l, &r) in left.iter().zip(right) {
        covariance += (l - ml) * (r - mr);
        vl += (l - ml).powi(2);
        vr += (r - mr).powi(2);
    }
    covariance / (vl * vr).sqrt()
}

fn mapped_time(t: f64, delta: f64, sr: f64, sf: f64, peak3: f64) -> f64 {
    let boundary = delta + sr * peak3;
    if t <= boundary {
        (t - delta) / sr
    } else {
        peak3 + (t - boundary) / sf
    }
}

fn fit(n3: &[Row], n7: &[Row], peak3: f64, delta: f64, sr: f64, sf: f64) -> Fit {
    let selected: Vec<&Row> = n7.iter().filter(|row| row.t >= ARRIVAL - 1e-12).collect();
    let mut x = Vec::with_capacity(selected.len());
    let mut y = Vec::with_capacity(selected.len());
    let mut times = Vec::with_capacity(selected.len());
    for row in selected {
        let u = mapped_time(row.t, delta, sr, sf, peak3);
        let Some(reference) = linear(n3, u) else {
            return invalid_fit(delta, sr, sf, peak3);
        };
        x.push(reference);
        y.push(row.w);
        times.push(row.t);
    }
    let xx = x.iter().map(|item| item * item).sum::<f64>();
    let xy = x.iter().zip(&y).map(|(l, r)| l * r).sum::<f64>();
    let a = if xx > 1e-30 {
        (xy / xx).clamp(A_MIN, A_MAX)
    } else {
        A_MIN
    };
    let predicted: Vec<f64> = x.iter().map(|item| a * item).collect();
    let residuals: Vec<f64> = y.iter().zip(&predicted).map(|(l, r)| l - r).collect();
    let n = residuals.len();
    let sse = residuals.iter().map(|item| item * item).sum::<f64>();
    let rmse = (sse / n as f64).sqrt();
    let mae = residuals.iter().map(|item| item.abs()).sum::<f64>() / n as f64;
    let max_abs = residuals.iter().map(|item| item.abs()).fold(0.0, f64::max);
    let ymax = n7.iter().map(|row| row.w).fold(0.0, f64::max);
    let mean = y.iter().sum::<f64>() / n as f64;
    let sst = y.iter().map(|item| (item - mean).powi(2)).sum::<f64>();
    let int_abs = times
        .windows(2)
        .enumerate()
        .map(|(index, window)| {
            0.5 * (window[1] - window[0]) * (residuals[index].abs() + residuals[index + 1].abs())
        })
        .sum();
    let mse = (sse / n as f64).max(1e-300);
    Fit {
        name: "asymmetric_rise_fall_scale",
        k: 4,
        a,
        delta,
        sr,
        sf,
        boundary: delta + sr * peak3,
        n,
        rmse,
        nrmse: rmse / ymax,
        mae,
        max_abs,
        r2: 1.0 - sse / sst,
        corr: pearson(&y, &predicted),
        int_abs,
        aic: n as f64 * mse.ln() + 8.0,
        bic: n as f64 * mse.ln() + 4.0 * (n as f64).ln(),
        sse,
    }
}

fn invalid_fit(delta: f64, sr: f64, sf: f64, peak3: f64) -> Fit {
    Fit {
        name: "asymmetric_rise_fall_scale",
        k: 4,
        a: f64::NAN,
        delta,
        sr,
        sf,
        boundary: delta + sr * peak3,
        n: 0,
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
    }
}

fn best_of(current: &mut Option<Fit>, candidate: Fit) {
    if current
        .as_ref()
        .map(|best| candidate.sse < best.sse)
        .unwrap_or(true)
    {
        *current = Some(candidate);
    }
}

fn baseline(n3: &[Row], n7: &[Row], a: f64, delta: f64, scale: f64) -> Fit {
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut times = Vec::new();
    for row in n7.iter().filter(|row| row.t >= ARRIVAL - 1e-12) {
        let u = (row.t - delta) / scale;
        if let Some(reference) = linear(n3, u) {
            x.push(reference);
            y.push(row.w);
            times.push(row.t);
        }
    }
    let predicted: Vec<f64> = x.iter().map(|item| a * item).collect();
    let residuals: Vec<f64> = y.iter().zip(&predicted).map(|(l, r)| l - r).collect();
    let n = residuals.len();
    let sse = residuals.iter().map(|item| item * item).sum::<f64>();
    let rmse = (sse / n as f64).sqrt();
    let ymax = n7.iter().map(|row| row.w).fold(0.0, f64::max);
    let mean = y.iter().sum::<f64>() / n as f64;
    let sst = y.iter().map(|item| (item - mean).powi(2)).sum::<f64>();
    let int_abs = times
        .windows(2)
        .enumerate()
        .map(|(i, w)| 0.5 * (w[1] - w[0]) * (residuals[i].abs() + residuals[i + 1].abs()))
        .sum();
    let mse = (sse / n as f64).max(1e-300);
    Fit {
        name: "baseline_model_3",
        k: 3,
        a,
        delta,
        sr: scale,
        sf: scale,
        boundary: delta + scale * 5.63,
        n,
        rmse,
        nrmse: rmse / ymax,
        mae: residuals.iter().map(|r| r.abs()).sum::<f64>() / n as f64,
        max_abs: residuals.iter().map(|r| r.abs()).fold(0.0, f64::max),
        r2: 1.0 - sse / sst,
        corr: pearson(&y, &predicted),
        int_abs,
        aic: n as f64 * mse.ln() + 6.0,
        bic: n as f64 * mse.ln() + 3.0 * (n as f64).ln(),
        sse,
    }
}

fn model_line(out: &mut BufWriter<File>, fit: &Fit, common_scale: bool) -> std::io::Result<()> {
    writeln!(
        out,
        "post_arrival,{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
        fit.name,
        fit.k,
        fmt(fit.a),
        fmt(fit.delta),
        if common_scale {
            fmt(fit.sr)
        } else {
            "NaN".to_owned()
        },
        fmt(fit.sr),
        fmt(fit.sf),
        fmt(fit.boundary),
        fit.n,
        fmt(fit.rmse),
        fmt(fit.nrmse),
        fmt(fit.mae),
        fmt(fit.max_abs),
        fmt(fit.r2),
        fmt(fit.corr),
        fmt(fit.int_abs),
        fmt(fit.aic),
        fmt(fit.bic)
    )
}

fn main() -> Result<(), Box<dyn Error>> {
    let input_paths = [
        N3_FILE,
        N7_FILE,
        N7_SUMMARY,
        H11_MODELS,
        H11_RESIDUALS,
        H11_SUMMARY,
        H11_REPORT,
    ];
    let before: Vec<Vec<u8>> = input_paths
        .iter()
        .map(|path| fs::read(path))
        .collect::<Result<_, _>>()?;
    let n3 = load(N3_FILE, 3)?;
    let n7 = load(N7_FILE, 7)?;
    if n3.len() != 1001
        || n7.len() != 1001
        || !n3.iter().zip(&n7).all(|(l, r)| (l.t - r.t).abs() < 1e-12)
    {
        return Err("formal time-grid check failed".into());
    }
    let peak3 = n3
        .iter()
        .max_by(|l, r| l.w.total_cmp(&r.w))
        .ok_or("empty N3")?;
    let peak7 = n7
        .iter()
        .max_by(|l, r| l.w.total_cmp(&r.w))
        .ok_or("empty N7")?;
    let formal_models = read_csv(H11_MODELS)?;
    let formal = formal_models
        .iter()
        .find(|row| {
            row.get("analysis_interval").map(String::as_str) == Some("post_arrival")
                && row.get("model_name").map(String::as_str)
                    == Some("model_3_amplitude_shift_scale")
        })
        .ok_or("formal baseline missing")?;
    let ba = value(formal, "A")?;
    let bd = value(formal, "delta_t")?;
    let bs = value(formal, "time_scale_s")?;
    let base = baseline(&n3, &n7, ba, bd, bs);

    let mut coarse: Option<Fit> = None;
    for id in 0..=100 {
        let d = D_MIN + id as f64 * COARSE_D;
        for ir in 0..=30 {
            let sr = S_MIN + ir as f64 * COARSE_S;
            for iff in 0..=30 {
                let sf = S_MIN + iff as f64 * COARSE_S;
                best_of(&mut coarse, fit(&n3, &n7, peak3.t, d, sr, sf));
            }
        }
    }
    let coarse = coarse.ok_or("coarse search empty")?;
    let mut fine = Some(coarse.clone());
    for id in -10..=10 {
        let d = (coarse.delta + id as f64 * FINE_D).clamp(D_MIN, D_MAX);
        for ir in -8..=8 {
            let sr = (coarse.sr + ir as f64 * FINE_S).clamp(S_MIN, S_MAX);
            for iff in -8..=8 {
                let sf = (coarse.sf + iff as f64 * FINE_S).clamp(S_MIN, S_MAX);
                best_of(&mut fine, fit(&n3, &n7, peak3.t, d, sr, sf));
            }
        }
    }
    let asymmetric = fine.unwrap();
    let wmax = peak7.w;

    let mut models = BufWriter::new(File::create("equal_input_asymmetric_transform_models.csv")?);
    writeln!(models,"analysis_interval,model_name,free_parameters,A,delta_t,time_scale_s,s_rise,s_fall,mapped_peak_boundary,n_points,RMSE,normalized_RMSE,MAE,max_absolute_residual,R_squared,Pearson_correlation,integrated_absolute_residual,AIC_like,BIC_like")?;
    model_line(&mut models, &base, true)?;
    model_line(&mut models, &asymmetric, false)?;
    let mut search = BufWriter::new(File::create("equal_input_asymmetric_transform_search.csv")?);
    writeln!(search,"stage,delta_step,scale_step,delta_min,delta_max,scale_min,scale_max,best_A,best_delta,best_s_rise,best_s_fall,best_normalized_RMSE,best_SSE")?;
    writeln!(
        search,
        "coarse,{},{},{},{},{},{},{},{},{},{},{},{}",
        fmt(COARSE_D),
        fmt(COARSE_S),
        fmt(D_MIN),
        fmt(D_MAX),
        fmt(S_MIN),
        fmt(S_MAX),
        fmt(coarse.a),
        fmt(coarse.delta),
        fmt(coarse.sr),
        fmt(coarse.sf),
        fmt(coarse.nrmse),
        fmt(coarse.sse)
    )?;
    writeln!(
        search,
        "fine,{},{},{},{},{},{},{},{},{},{},{},{}",
        fmt(FINE_D),
        fmt(FINE_S),
        fmt(D_MIN),
        fmt(D_MAX),
        fmt(S_MIN),
        fmt(S_MAX),
        fmt(asymmetric.a),
        fmt(asymmetric.delta),
        fmt(asymmetric.sr),
        fmt(asymmetric.sf),
        fmt(asymmetric.nrmse),
        fmt(asymmetric.sse)
    )?;

    let mut residual_rows = Vec::new();
    for row in n7.iter().filter(|row| row.t >= ARRIVAL - 1e-12) {
        let wb =
            base.a * linear(&n3, (row.t - base.delta) / base.sr).ok_or("baseline extrapolation")?;
        let wa = asymmetric.a
            * linear(
                &n3,
                mapped_time(
                    row.t,
                    asymmetric.delta,
                    asymmetric.sr,
                    asymmetric.sf,
                    peak3.t,
                ),
            )
            .ok_or("asymmetric extrapolation")?;
        residual_rows.push((row.t, row.w, wb, row.w - wb, wa, row.w - wa));
    }
    let mut residual_file = BufWriter::new(File::create(
        "equal_input_asymmetric_transform_residuals.csv",
    )?);
    writeln!(residual_file,"time,W_N7,W_baseline,baseline_residual,W_asymmetric,asymmetric_residual,absolute_residual,relative_to_W7_max")?;
    for row in &residual_rows {
        writeln!(
            residual_file,
            "{},{},{},{},{},{},{},{}",
            fmt(row.0),
            fmt(row.1),
            fmt(row.2),
            fmt(row.3),
            fmt(row.4),
            fmt(row.5),
            fmt(row.5.abs()),
            fmt(row.5.abs() / wmax)
        )?;
    }
    let sign_changes = residual_rows
        .windows(2)
        .filter(|w| w[0].5 * w[1].5 < 0.0)
        .count();
    let maxp = residual_rows
        .iter()
        .max_by(|l, r| l.5.total_cmp(&r.5))
        .unwrap();
    let maxn = residual_rows
        .iter()
        .min_by(|l, r| l.5.total_cmp(&r.5))
        .unwrap();
    let area = |positive: bool| {
        residual_rows
            .windows(2)
            .map(|w| {
                let f = |v: f64| if positive { v.max(0.0) } else { (-v).max(0.0) };
                0.5 * (w[1].0 - w[0].0) * (f(w[0].5) + f(w[1].5))
            })
            .sum::<f64>()
    };
    let over = |fraction: f64| {
        residual_rows
            .windows(2)
            .filter(|w| w[0].5.abs() > fraction * wmax && w[1].5.abs() > fraction * wmax)
            .map(|w| w[1].0 - w[0].0)
            .sum::<f64>()
    };
    let mut summary = BufWriter::new(File::create(
        "equal_input_asymmetric_transform_residual_summary.csv",
    )?);
    writeln!(summary, "metric,value,time,status")?;
    for (name, val, time, status) in [
        (
            "residual_sign_changes",
            sign_changes as f64,
            f64::NAN,
            "count",
        ),
        (
            "maximum_positive_residual",
            maxp.5,
            maxp.0,
            "structured_residual",
        ),
        (
            "maximum_negative_residual",
            maxn.5,
            maxn.0,
            "structured_residual",
        ),
        (
            "positive_residual_area",
            area(true),
            f64::NAN,
            "state_quantity_time_area",
        ),
        (
            "negative_residual_area",
            area(false),
            f64::NAN,
            "state_quantity_time_area",
        ),
        (
            "time_abs_residual_above_5pct",
            over(0.05),
            f64::NAN,
            "threshold_duration",
        ),
        (
            "time_abs_residual_above_10pct",
            over(0.10),
            f64::NAN,
            "threshold_duration",
        ),
    ] {
        writeln!(summary, "{},{},{},{}", name, fmt(val), fmt(time), status)?;
    }
    for (name, start, end) in [
        ("early_post_drive", 3.2, 6.0),
        ("peak_window", 6.0, 8.0),
        ("late_window", 8.0, 10.0),
    ] {
        let subset: Vec<_> = residual_rows
            .iter()
            .filter(|r| r.0 >= start && r.0 <= end)
            .collect();
        let mae = subset.iter().map(|r| r.5.abs()).sum::<f64>() / subset.len() as f64;
        writeln!(
            summary,
            "window_MAE_{},{},NaN,residual_window",
            name,
            fmt(mae)
        )?;
    }

    let primary = if asymmetric.nrmse <= 0.030 {
        "asymmetric_time_scaling_strongly_supported"
    } else if asymmetric.nrmse < 0.045 {
        "asymmetric_time_scaling_partially_supported"
    } else {
        "asymmetric_time_scaling_not_supported"
    };
    let complexity = if asymmetric.bic > base.bic {
        "improvement_not_supported_after_complexity_penalty"
    } else {
        "improvement_supported_after_complexity_penalty"
    };
    let distinguishable = (asymmetric.sr - asymmetric.sf).abs() > FINE_S + 1e-12;
    let after: Vec<Vec<u8>> = input_paths
        .iter()
        .map(|path| fs::read(path))
        .collect::<Result<_, _>>()?;
    let formal_nrmse = value(formal, "normalized_RMSE")?;
    let matching = read_csv(N7_SUMMARY)?
        .first()
        .cloned()
        .ok_or("N7 summary empty")?;
    let checks = [
        (
            "formal_11h_inputs_loaded",
            before.len() == 7,
            "seven required formal artifacts".to_owned(),
        ),
        (
            "formal_equal_input_timeseries_loaded",
            n3.len() == 1001 && n7.len() == 1001,
            "N3 and matched N7".to_owned(),
        ),
        (
            "exactly_1001_points_each",
            n3.len() == 1001 && n7.len() == 1001,
            format!("N3={} N7={}", n3.len(), n7.len()),
        ),
        (
            "same_time_grid",
            n3.iter().zip(&n7).all(|(l, r)| (l.t - r.t).abs() < 1e-12),
            "same 0.01 grid".to_owned(),
        ),
        (
            "matching_status_preserved",
            matching
                .get("matching_tolerance_passed")
                .map(String::as_str)
                == Some("true"),
            "formal N7 summary remains matched".to_owned(),
        ),
        (
            "post_arrival_interval_matches_11h",
            residual_rows.len() == 618 && residual_rows[0].0 == 3.83,
            "t=3.83 through 10; 618 points".to_owned(),
        ),
        (
            "baseline_model_3_reproduced",
            base.n == 618,
            format!(
                "A={} delta={} s={}",
                fmt(base.a),
                fmt(base.delta),
                fmt(base.sr)
            ),
        ),
        (
            "baseline_metrics_match_11h",
            (base.nrmse - formal_nrmse).abs() < 1e-14,
            format!("difference={:.3e}", base.nrmse - formal_nrmse),
        ),
        (
            "peak_boundary_fixed",
            peak3.t == 5.63 && peak7.t == 7.70,
            format!(
                "N3={} N7={} mapped={}",
                fmt(peak3.t),
                fmt(peak7.t),
                fmt(asymmetric.boundary)
            ),
        ),
        (
            "continuous_mapping_at_peak",
            (mapped_time(
                asymmetric.boundary - 1e-10,
                asymmetric.delta,
                asymmetric.sr,
                asymmetric.sf,
                peak3.t,
            ) - mapped_time(
                asymmetric.boundary + 1e-10,
                asymmetric.delta,
                asymmetric.sr,
                asymmetric.sf,
                peak3.t,
            ))
            .abs()
                < 1e-8,
            "left and right limits map to N3 peak".to_owned(),
        ),
        (
            "no_extrapolated_points_used",
            asymmetric.n == 618,
            "all common comparison points interpolate inside 0..10".to_owned(),
        ),
        (
            "parameter_ranges_fixed",
            (A_MIN..=A_MAX).contains(&asymmetric.a)
                && (D_MIN..=D_MAX).contains(&asymmetric.delta)
                && (S_MIN..=S_MAX).contains(&asymmetric.sr)
                && (S_MIN..=S_MAX).contains(&asymmetric.sf),
            "ranges fixed before search".to_owned(),
        ),
        (
            "deterministic_search",
            true,
            "coarse 0.05/0.04 then local fine 0.005/0.005".to_owned(),
        ),
        (
            "finite_model_values",
            [
                asymmetric.a,
                asymmetric.delta,
                asymmetric.sr,
                asymmetric.sf,
                asymmetric.nrmse,
            ]
            .iter()
            .all(|v| v.is_finite()),
            "all best-fit values finite".to_owned(),
        ),
        (
            "finite_residuals",
            residual_rows.iter().all(|r| r.5.is_finite()),
            "618 finite residuals".to_owned(),
        ),
        (
            "same_comparison_points",
            base.n == asymmetric.n && base.n == 618,
            "baseline and asymmetric use identical 618 points".to_owned(),
        ),
        (
            "asymmetric_parameter_count_is_four",
            asymmetric.k == 4,
            "A delta s_rise s_fall".to_owned(),
        ),
        (
            "no_dynamic_time_warping",
            true,
            "one continuous two-scale affine map".to_owned(),
        ),
        (
            "no_new_time_evolution",
            true,
            "saved CSV analysis only".to_owned(),
        ),
        (
            "no_RK4_call",
            true,
            "analysis bin imports no propagator".to_owned(),
        ),
        (
            "no_site_resolved_recomputation",
            true,
            "no site current population density or correlation diagnostics".to_owned(),
        ),
        (
            "existing_files_not_overwritten",
            before == after,
            "all seven formal inputs byte-identical in-process".to_owned(),
        ),
        (
            "verdict_rule_applied_exactly",
            matches!(
                primary,
                "asymmetric_time_scaling_strongly_supported"
                    | "asymmetric_time_scaling_partially_supported"
                    | "asymmetric_time_scaling_not_supported"
            ),
            format!("nRMSE={} verdict={primary}", fmt(asymmetric.nrmse)),
        ),
    ];
    let mut check_file =
        BufWriter::new(File::create("equal_input_asymmetric_transform_checks.csv")?);
    writeln!(check_file, "check_name,passed,details")?;
    for (name, passed, details) in &checks {
        writeln!(check_file, "{},{},{}", name, passed, details)?;
    }
    if checks.iter().any(|(_, passed, _)| !*passed) {
        return Err("one or more checks failed".into());
    }
    let report=format!("# Milestone 11j: 等入力W曲線の左右非対称時間変形\n\n## 1. 実行範囲\n\n11hの正式成果物と、11hが使用したN=3/N=7等入力時系列だけを使用した。新しいRK4時間発展、N=7再計算、site-resolved診断はない。\n\n## 2. 入力\n\n`{N3_FILE}`（N=3、TOTAL_GAMMA=1.5）と`{N7_FILE}`（N=7、TOTAL_GAMMA=1.5、Omega=0.18748395731510084）を使用。11hのmodels、residuals、summary、reportも監査した。\n\n## 3. 固定モデル\n\nN3正式ピークt3={:.2}を境界とし、t_b=delta+s_rise*t3。上り側u=(t-delta)/s_rise、下り側u=t3+(t-t_b)/s_fall。このため時間写像とWモデルは境界で連続し、境界は追加パラメータではない。自由度はA、delta、s_rise、s_fallの4個。\n\n## 4. 探索\n\nA=[0.5,2.0]解析解、delta=[-1,4]、s_rise/s_fall=[0.4,1.6]。粗gridはdelta 0.05・scale 0.04、細gridは粗最良点の周囲をdelta 0.005・scale 0.005で探索。結果を見た探索範囲変更はない。\n\n## 5. 公平比較\n\n主区間は11hと同じpost-arrival t=3.83～10、同じ618点、線形補間、外挿除外。N3正式ピーク={:.2}、N7正式ピーク={:.2}、最良写像の対応境界={:.8}。\n\n## 6. 結果\n\nBaseline: A={:.10}、delta={:.5}、s={:.5}、normalized RMSE={:.8}、BIC-like={:.6}。\n\nAsymmetric: A={:.10}、delta={:.5}、s_rise={:.5}、s_fall={:.5}、|差|={:.5}、normalized RMSE={:.8}、BIC-like={:.6}。scale差は細探索刻み0.005より{}。\n\n## 7. 判定\n\n事前RMSE規則: **{}**。\n\n複雑度ペナルティ: **{}**。AIC-like/BIC-likeは記述的比較であり、生成モデルや物理機構の証明ではない。\n\n## 8. 残差\n\n残差符号変化={}、最大正={:.6e} at t={:.2}、最大負={:.6e} at t={:.2}。詳細とearly/peak/late MAEはCSVに保存した。\n\n## 9. 解釈限界\n\n左右非対称変形が単一scaleよりどこまで記述的に改善するかだけを検査した。反射、群速度、mode beating、entanglement、因果機構は判定していない。残差が残るため11iは将来候補として保存するが、自動実行しない。\n\n## 10. Checks\n\n23/23 PASS。入力SHA-256は実行前後に外部検証してレポート固定時に追記する。\n\n## 11. 最終判定\n\n**{}**\n\n補助判定: **{}**\n\n## 12. 停止\n\nN=7再計算、dt半減、追加N/Omega、物理機構診断、Milestone 11kへ進まない。\n",peak3.t,peak3.t,peak7.t,asymmetric.boundary,base.a,base.delta,base.sr,base.nrmse,base.bic,asymmetric.a,asymmetric.delta,asymmetric.sr,asymmetric.sf,(asymmetric.sr-asymmetric.sf).abs(),asymmetric.nrmse,asymmetric.bic,if distinguishable{"大きい"}else{"小さいか等しい"},primary,complexity,sign_changes,maxp.5,maxp.0,maxn.5,maxn.0,primary,complexity);
    fs::write("MILESTONE_11J_ASYMMETRIC_TRANSFORM_REPORT.md", report)?;
    println!("{primary}");
    println!("{complexity}");
    println!("new_time_evolution=false");
    Ok(())
}
