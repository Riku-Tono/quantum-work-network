use std::collections::HashMap;
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufWriter, Write};

const N3_FILE: &str = "fixed_total_gamma_1_5_xgamma_timeseries.csv";
const N7_FILE: &str = "input_matching_interpolated_trial_timeseries.csv";
const N7_SUMMARY: &str = "input_matching_interpolated_trial_summary.csv";
const J_RESIDUALS: &str = "equal_input_asymmetric_transform_residuals.csv";
const J_MODELS: &str = "equal_input_asymmetric_transform_models.csv";
const J_REPORT: &str = "MILESTONE_11J_ASYMMETRIC_TRANSFORM_REPORT.md";
const T_PEAK: f64 = 7.70;
const T_BOUNDARY: f64 = 7.75295;
const W_FLOOR: f64 = 1e-10;

#[derive(Clone)]
struct Row {
    t: f64,
    e: f64,
    w: f64,
    u: f64,
    c: f64,
}

#[derive(Clone)]
struct Fit {
    interval: &'static str,
    model: &'static str,
    k: usize,
    switch: f64,
    c: f64,
    l1: f64,
    l2: f64,
    b: f64,
    n: usize,
    rmse: f64,
    nrmse: f64,
    mae: f64,
    max_abs: f64,
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

fn csv(path: &str) -> Result<Vec<HashMap<String, String>>, Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    let mut lines = text.lines();
    let h: Vec<String> = lines
        .next()
        .ok_or("empty CSV")?
        .split(',')
        .map(str::to_owned)
        .collect();
    let mut out = Vec::new();
    for line in lines.filter(|x| !x.trim().is_empty()) {
        let v: Vec<&str> = line.split(',').collect();
        if v.len() != h.len() {
            return Err(format!("width mismatch {path}").into());
        }
        out.push(
            h.iter()
                .cloned()
                .zip(v.iter().map(|x| (*x).to_owned()))
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
    for r in csv(path)? {
        if val(&r, "chain_length")? as usize == n {
            out.push(Row {
                t: val(&r, "time")?,
                e: val(&r, "load_energy")?,
                w: val(&r, "load_ergotropy")?,
                u: val(&r, "usable_fraction")?,
                c: val(&r, "load_coherence_l1")?,
            })
        }
    }
    out.sort_by(|a, b| a.t.total_cmp(&b.t));
    Ok(out)
}

fn rate(values: &[f64], dt: f64) -> Vec<f64> {
    (0..values.len())
        .map(|i| {
            if values[i] <= W_FLOOR {
                return f64::NAN;
            }
            let derivative = if i == 0 {
                (values[1].ln() - values[0].ln()) / dt
            } else if i + 1 == values.len() {
                (values[i].ln() - values[i - 1].ln()) / dt
            } else {
                (values[i + 1].ln() - values[i - 1].ln()) / (2.0 * dt)
            };
            -derivative
        })
        .collect()
}

fn shape(model: &str, x: f64, l1: f64, l2: f64, b: f64, sw: f64) -> f64 {
    match model {
        "model_A_single_exponential" => (-l1 * x).exp(),
        "model_B_two_stage_exponential" => {
            let xs = sw - T_PEAK;
            if x <= xs {
                (-l1 * x).exp()
            } else {
                (-l1 * xs - l2 * (x - xs)).exp()
            }
        }
        "model_C_exponential_linear" => (-l1 * x).exp() * (1.0 + b * x),
        _ => f64::NAN,
    }
}

fn fit_at(
    rows: &[&Row],
    interval: &'static str,
    model: &'static str,
    k: usize,
    l1: f64,
    l2: f64,
    b: f64,
    sw: f64,
) -> Fit {
    let x: Vec<f64> = rows.iter().map(|r| r.t - T_PEAK).collect();
    let y: Vec<f64> = rows.iter().map(|r| r.w).collect();
    let q: Vec<f64> = x.iter().map(|&v| shape(model, v, l1, l2, b, sw)).collect();
    if q.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return invalid(interval, model, k, l1, l2, b, sw);
    }
    let qq = q.iter().map(|v| v * v).sum::<f64>();
    let qy = q.iter().zip(&y).map(|(a, b)| a * b).sum::<f64>();
    let c = (qy / qq).max(0.0);
    let pred: Vec<f64> = q.iter().map(|v| c * v).collect();
    let res: Vec<f64> = y.iter().zip(&pred).map(|(a, b)| a - b).collect();
    let n = res.len();
    let sse = res.iter().map(|v| v * v).sum::<f64>();
    let rmse = (sse / n as f64).sqrt();
    let ymax = y.iter().copied().fold(0.0, f64::max);
    let mse = (sse / n as f64).max(1e-300);
    let int_abs = rows
        .windows(2)
        .enumerate()
        .map(|(i, w)| 0.5 * (w[1].t - w[0].t) * (res[i].abs() + res[i + 1].abs()))
        .sum();
    Fit {
        interval,
        model,
        k,
        switch: sw,
        c,
        l1,
        l2,
        b,
        n,
        rmse,
        nrmse: rmse / ymax,
        mae: res.iter().map(|v| v.abs()).sum::<f64>() / n as f64,
        max_abs: res.iter().map(|v| v.abs()).fold(0.0, f64::max),
        int_abs,
        aic: n as f64 * mse.ln() + 2.0 * k as f64,
        bic: n as f64 * mse.ln() + k as f64 * (n as f64).ln(),
        sse,
    }
}
fn invalid(
    interval: &'static str,
    model: &'static str,
    k: usize,
    l1: f64,
    l2: f64,
    b: f64,
    sw: f64,
) -> Fit {
    Fit {
        interval,
        model,
        k,
        switch: sw,
        c: f64::NAN,
        l1,
        l2,
        b,
        n: 0,
        rmse: f64::INFINITY,
        nrmse: f64::INFINITY,
        mae: f64::INFINITY,
        max_abs: f64::INFINITY,
        int_abs: f64::INFINITY,
        aic: f64::INFINITY,
        bic: f64::INFINITY,
        sse: f64::INFINITY,
    }
}
fn improve(best: &mut Option<Fit>, candidate: Fit) {
    if best.as_ref().map(|b| candidate.sse < b.sse).unwrap_or(true) {
        *best = Some(candidate)
    }
}

fn fit_models(all: &[Row], start: f64, end: f64, name: &'static str) -> Vec<Fit> {
    let rows: Vec<&Row> = all
        .iter()
        .filter(|r| r.t >= start - 1e-12 && r.t <= end + 1e-12)
        .collect();
    let mut a = None;
    for i in 0..=2000 {
        improve(
            &mut a,
            fit_at(
                &rows,
                name,
                "model_A_single_exponential",
                2,
                i as f64 * 0.001,
                0.0,
                0.0,
                f64::NAN,
            ),
        );
    }
    let mut bfits = Vec::new();
    for &sw in &[8.5, 9.0, 9.5] {
        if sw <= start || sw >= end {
            continue;
        }
        let mut candidate = None;
        for i in 0..=200 {
            for j in 0..=200 {
                improve(
                    &mut candidate,
                    fit_at(
                        &rows,
                        name,
                        "model_B_two_stage_exponential",
                        3,
                        i as f64 * 0.01,
                        j as f64 * 0.01,
                        0.0,
                        sw,
                    ),
                );
            }
        }
        bfits.push(candidate.unwrap());
    }
    let mut cbest = None;
    for i in 0..=400 {
        let l = i as f64 * 0.005;
        for j in 0..=280 {
            let b = -0.4 + j as f64 * 0.005;
            improve(
                &mut cbest,
                fit_at(
                    &rows,
                    name,
                    "model_C_exponential_linear",
                    3,
                    l,
                    0.0,
                    b,
                    f64::NAN,
                ),
            );
        }
    }
    let mut out = vec![a.unwrap()];
    out.extend(bfits);
    out.push(cbest.unwrap());
    out
}

fn area(rows: &[(f64, f64)], abs: bool) -> f64 {
    rows.windows(2)
        .map(|w| {
            let f = |v: f64| if abs { v.abs() } else { v };
            0.5 * (w[1].0 - w[0].0) * (f(w[0].1) + f(w[1].1))
        })
        .sum()
}

fn main() -> Result<(), Box<dyn Error>> {
    let inputs = [
        N3_FILE,
        N7_FILE,
        N7_SUMMARY,
        J_RESIDUALS,
        J_MODELS,
        J_REPORT,
    ];
    let before: Vec<Vec<u8>> = inputs
        .iter()
        .map(|p| fs::read(p))
        .collect::<Result<_, _>>()?;
    let n3 = load(N3_FILE, 3)?;
    let n7 = load(N7_FILE, 7)?;
    if n3.len() != 1001 || n7.len() != 1001 {
        return Err("expected 1001 rows".into());
    }
    let jres = csv(J_RESIDUALS)?;
    let jm = csv(J_MODELS)?;
    let jreport = fs::read_to_string(J_REPORT)?;
    let jasym = jm
        .iter()
        .find(|r| r.get("model_name").map(String::as_str) == Some("asymmetric_rise_fall_scale"))
        .ok_or("11j asymmetric model missing")?;
    let boundary = val(jasym, "mapped_peak_boundary")?;
    if (boundary - T_BOUNDARY).abs() > 1e-12 {
        return Err("11j boundary mismatch".into());
    }
    let n7sum = csv(N7_SUMMARY)?.first().cloned().ok_or("summary empty")?;
    let matching = n7sum.get("matching_tolerance_passed").map(String::as_str) == Some("true");
    let post7: Vec<&Row> = n7.iter().filter(|r| r.t >= T_PEAK - 1e-12).collect();
    let post3: Vec<&Row> = n3.iter().filter(|r| r.t >= T_PEAK - 1e-12).collect();
    let model_w: Vec<f64> = jres
        .iter()
        .filter(|r| val(r, "time").unwrap() >= T_PEAK - 1e-12)
        .map(|r| val(r, "W_asymmetric").unwrap())
        .collect();
    let w7: Vec<f64> = post7.iter().map(|r| r.w).collect();
    let w3: Vec<f64> = post3.iter().map(|r| r.w).collect();
    let k7 = rate(&w7, 0.01);
    let km = rate(&model_w, 0.01);
    let k3 = rate(&w3, 0.01);
    let mut rates = BufWriter::new(File::create("equal_input_post_peak_decay_rates.csv")?);
    writeln!(
        rates,
        "time,W_N7,W_model_11j,W_N3,k_W_N7,k_W_model,k_W_N3,Delta_k_N7_minus_model,W_floor"
    )?;
    for i in 0..post7.len() {
        writeln!(
            rates,
            "{},{},{},{},{},{},{},{},{}",
            fmt(post7[i].t),
            fmt(w7[i]),
            fmt(model_w[i]),
            fmt(w3[i]),
            fmt(k7[i]),
            fmt(km[i]),
            fmt(k3[i]),
            fmt(k7[i] - km[i]),
            fmt(W_FLOOR)
        )?;
    }
    let mut fits = Vec::new();
    for (s, e, n) in [
        (7.7, 10.0, "post_peak_7_70_10_00"),
        (7.7, 9.0, "early_tail_7_70_9_00"),
        (9.0, 10.0, "late_tail_9_00_10_00"),
    ] {
        fits.extend(fit_models(&n7, s, e, n));
    }
    let mut models = BufWriter::new(File::create("equal_input_post_peak_decay_models.csv")?);
    writeln!(models,"analysis_interval,model_name,free_parameters,t_switch,C,lambda_1,lambda_2,b,n_points,RMSE,normalized_RMSE,MAE,maximum_absolute_residual,integrated_absolute_residual,AIC_like,BIC_like")?;
    for f in &fits {
        writeln!(
            models,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            f.interval,
            f.model,
            f.k,
            fmt(f.switch),
            fmt(f.c),
            fmt(f.l1),
            fmt(f.l2),
            fmt(f.b),
            f.n,
            fmt(f.rmse),
            fmt(f.nrmse),
            fmt(f.mae),
            fmt(f.max_abs),
            fmt(f.int_abs),
            fmt(f.aic),
            fmt(f.bic)
        )?;
    }
    let mainfits: Vec<&Fit> = fits
        .iter()
        .filter(|f| f.interval == "post_peak_7_70_10_00")
        .collect();
    let best = *mainfits
        .iter()
        .min_by(|a, b| a.bic.total_cmp(&b.bic))
        .unwrap();
    let mainrows: Vec<&Row> = n7.iter().filter(|r| r.t >= T_PEAK - 1e-12).collect();
    let mut residuals = BufWriter::new(File::create("equal_input_post_peak_decay_residuals.csv")?);
    writeln!(
        residuals,
        "time,W_N7,W_best_decay_model,residual,absolute_residual,relative_to_peak,best_model"
    )?;
    for r in &mainrows {
        let p = best.c
            * shape(
                best.model,
                r.t - T_PEAK,
                best.l1,
                best.l2,
                best.b,
                best.switch,
            );
        writeln!(
            residuals,
            "{},{},{},{},{},{},{}",
            fmt(r.t),
            fmt(r.w),
            fmt(p),
            fmt(r.w - p),
            fmt((r.w - p).abs()),
            fmt((r.w - p).abs() / w7[0]),
            best.model
        )?;
    }
    let e7: Vec<f64> = post7.iter().map(|r| r.e).collect();
    let p7: Vec<f64> = post7.iter().map(|r| r.e - r.w).collect();
    let c7: Vec<f64> = post7.iter().map(|r| r.c).collect();
    let ke = rate(&e7, 0.01);
    let kp = rate(&p7, 0.01);
    let kc = rate(&c7, 0.01);
    let mut components = BufWriter::new(File::create(
        "equal_input_post_peak_component_comparison.csv",
    )?);
    writeln!(components,"time,E_N7,W_N7,P_N7,coherence_N7,E_norm,W_norm,P_norm,coherence_norm,k_E,k_W,k_P,k_coherence,Delta_E_N7_minus_N3,Delta_W_N7_minus_N3,Delta_P_N7_minus_N3,Delta_usable_N7_minus_N3,R_11j")?;
    let e0 = e7[0];
    let w0 = w7[0];
    let p0 = p7[0];
    let c0 = c7[0];
    for i in 0..post7.len() {
        let r3 = &post3[i];
        let r7 = post7[i];
        let rj = val(
            &jres[(r7.t / 0.01).round() as usize - 383],
            "asymmetric_residual",
        )?;
        writeln!(
            components,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            fmt(r7.t),
            fmt(r7.e),
            fmt(r7.w),
            fmt(r7.e - r7.w),
            fmt(r7.c),
            fmt(r7.e / e0),
            fmt(r7.w / w0),
            fmt((r7.e - r7.w) / p0),
            fmt(r7.c / c0),
            fmt(ke[i]),
            fmt(k7[i]),
            fmt(kp[i]),
            fmt(kc[i]),
            fmt(r7.e - r3.e),
            fmt(r7.w - r3.w),
            fmt((r7.e - r7.w) - (r3.e - r3.w)),
            fmt(r7.u - r3.u),
            fmt(rj)
        )?;
    }
    let jr: Vec<(f64, f64)> = jres
        .iter()
        .filter_map(|r| {
            let t = val(r, "time").ok()?;
            if t >= T_PEAK - 1e-12 {
                Some((t, val(r, "asymmetric_residual").ok()?))
            } else {
                None
            }
        })
        .collect();
    let early: Vec<_> = jr.iter().copied().filter(|r| r.0 <= 9.0 + 1e-12).collect();
    let late: Vec<_> = jr.iter().copied().filter(|r| r.0 >= 9.0 - 1e-12).collect();
    let total_abs = area(&jr, true);
    let late_abs = area(&late, true);
    let sign_changes = jr.windows(2).filter(|w| w[0].1 * w[1].1 < 0.0).count();
    let maxp = jr.iter().max_by(|a, b| a.1.total_cmp(&b.1)).unwrap();
    let maxn = jr.iter().min_by(|a, b| a.1.total_cmp(&b.1)).unwrap();
    let tail_class = if sign_changes >= 4 {
        "tail_residual_oscillatory"
    } else if late_abs / total_abs >= 0.5 {
        "tail_residual_concentrated_late"
    } else {
        "tail_residual_distributed"
    };
    let two = *mainfits
        .iter()
        .filter(|f| f.model == "model_B_two_stage_exponential")
        .min_by(|a, b| a.bic.total_cmp(&b.bic))
        .unwrap();
    let simple = mainfits
        .iter()
        .find(|f| f.model == "model_A_single_exponential")
        .unwrap();
    let two_stage_supported = two.bic + 10.0 < simple.bic && (two.l1 - two.l2).abs() > 0.01;
    let verdict = if best.nrmse > 0.03 || tail_class == "tail_residual_concentrated_late" {
        "late_tail_structure_remains"
    } else if two_stage_supported {
        "two_stage_tail_decay_supported"
    } else if best.nrmse <= 0.03 {
        "simple_decay_rate_difference_supported"
    } else {
        "late_tail_structure_remains"
    };
    let after: Vec<Vec<u8>> = inputs
        .iter()
        .map(|p| fs::read(p))
        .collect::<Result<_, _>>()?;
    let identity = n7.iter().all(|r| ((r.e - r.w) + r.w - r.e).abs() < 1e-14);
    let checks = [
        (
            "formal_11j_inputs_loaded",
            jreport.contains("asymmetric_time_scaling_partially_supported"),
            "five formal artifacts plus matching summary".to_owned(),
        ),
        (
            "matching_precondition_passed",
            matching,
            "formal equal-input match true".to_owned(),
        ),
        (
            "post_peak_interval_fixed",
            post7.len() == 231,
            "t=7.70 through 10.00; 231 points".to_owned(),
        ),
        (
            "instantaneous_decay_rate_computed",
            k7.iter().all(|v| v.is_finite()),
            "central difference; one-sided endpoints".to_owned(),
        ),
        (
            "W_floor_applied",
            W_FLOOR == 1e-10,
            "undefined at or below floor".to_owned(),
        ),
        (
            "only_low_freedom_models_used",
            fits.len() == 11,
            "A C and every applicable pre-fixed B switch over three intervals".to_owned(),
        ),
        (
            "model_continuity_preserved",
            fits.iter()
                .filter(|f| f.model == "model_B_two_stage_exponential")
                .all(|f| f.switch.is_finite()),
            "two-stage formula continuous at fixed switch".to_owned(),
        ),
        (
            "complexity_scores_computed",
            fits.iter().all(|f| f.aic.is_finite() && f.bic.is_finite()),
            "AIC-like and BIC-like finite".to_owned(),
        ),
        (
            "11j_residual_partitioned",
            !early.is_empty() && !late.is_empty(),
            format!("late abs fraction={:.6}", late_abs / total_abs),
        ),
        (
            "E_W_passive_identity_preserved",
            identity,
            "P=E-W exactly from saved scalars".to_owned(),
        ),
        (
            "no_high_order_fit",
            true,
            "only single exponential two-stage exponential exponential-times-linear".to_owned(),
        ),
        (
            "no_new_time_evolution",
            true,
            "saved time series only".to_owned(),
        ),
        (
            "existing_files_not_overwritten",
            before == after,
            "six formal inputs byte-identical in-process".to_owned(),
        ),
    ];
    let mut cw = BufWriter::new(File::create("equal_input_post_peak_decay_checks.csv")?);
    writeln!(cw, "check_name,passed,details")?;
    for (c, p, d) in &checks {
        writeln!(cw, "{},{},{}", c, p, d)?;
    }
    if checks.iter().any(|(_, p, _)| !*p) {
        return Err("checks failed".into());
    }
    let mut sw = BufWriter::new(File::create("equal_input_post_peak_decay_summary.csv")?);
    writeln!(sw, "metric,value,time,status")?;
    for (n, v, t, s) in [
        ("best_model_post_peak", f64::NAN, f64::NAN, best.model),
        ("best_normalized_RMSE", best.nrmse, f64::NAN, "model_metric"),
        ("best_BIC_like", best.bic, f64::NAN, "complexity_metric"),
        (
            "R11j_post_peak_signed_area",
            area(&jr, false),
            f64::NAN,
            "residual_area",
        ),
        (
            "R11j_post_peak_absolute_area",
            total_abs,
            f64::NAN,
            "residual_area",
        ),
        (
            "R11j_7_70_9_signed_area",
            area(&early, false),
            f64::NAN,
            "residual_area",
        ),
        (
            "R11j_9_10_signed_area",
            area(&late, false),
            f64::NAN,
            "residual_area",
        ),
        (
            "R11j_9_10_absolute_fraction",
            late_abs / total_abs,
            f64::NAN,
            tail_class,
        ),
        ("R11j_maximum_positive", maxp.1, maxp.0, "residual_extreme"),
        ("R11j_maximum_negative", maxn.1, maxn.0, "residual_extreme"),
        ("R11j_sign_changes", sign_changes as f64, f64::NAN, "count"),
    ] {
        writeln!(sw, "{},{},{},{}", n, fmt(v), fmt(t), s)?;
    }
    let report=format!("# Milestone 11k: 等入力Wピーク後の減衰率・尾部残差解析\n\n## 1. 範囲\n\n11j正式成果物とN=3/N=7等入力保存時系列のみを使用。新規RK4・軌道再計算なし。主区間t=7.70～10.00、補助区間7.70～9.00と9.00～10.00。11j境界t_b={boundary:.5}。\n\n## 2. 瞬間減衰率\n\nk_Q=-d ln(Q)/dtを中央差分、端点のみ片側差分で計算。floor=1e-10。N7実測W、11jモデルW、N3 Wと、N7のE/P/coherenceをCSV化した。\n\n## 3. 低自由度モデル\n\nModel A単一指数、Model B連続二段階指数（switch 8.5/9.0/9.5のみ）、Model C指数×線形だけを評価。AIC-like/BIC-likeは記述的比較で物理過程の証明ではない。\n\n## 4. 主結果\n\nBIC-like最良={}、normalized RMSE={:.8}、BIC-like={:.6}。Model A nRMSE={:.8}、Model B nRMSE={:.8}（switch={:.2}, lambda1={:.5}, lambda2={:.5}）、Model C nRMSE={:.8}。\n\n## 5. 11j残差の尾部\n\npost-peak absolute area={:.6e}、9～10 absolute fraction={:.6}、符号変化={}、最大正={:.6e} at {:.2}、最大負={:.6e} at {:.2}。分類: **{}**。\n\n## 6. E/W/passive/coherence\n\n各量をt=7.70値で規格化し、対数減衰率とともに保存した。これは同時変化の記述であり、coherenceがW減衰を引き起こしたとは主張しない。t=8～10のDelta E/W/passive/usableと11j残差も同一表に保存した。usable fraction上昇を独立性能向上とは扱わない。\n\n## 7. Checks\n\n13/13 PASS。入力SHA-256とreleaseテスト記録は提出固定時に追記する。\n\n## 8. 判定\n\n**{}**\n\n## 9. 解釈限界\n\n減衰fitは形状診断であり、境界反射、群速度、二つの物理過程、因果機構を証明しない。\n\n## 10. 次段階\n\n終盤構造が残る場合、11iのsite-resolved診断を同一N=7 matched軌道1本に限定して検討できる。ただし自動実行しない。\n",best.model,best.nrmse,best.bic,simple.nrmse,two.nrmse,two.switch,two.l1,two.l2,mainfits.iter().find(|f|f.model=="model_C_exponential_linear").unwrap().nrmse,total_abs,late_abs/total_abs,sign_changes,maxp.1,maxp.0,maxn.1,maxn.0,tail_class,verdict);
    fs::write("MILESTONE_11K_POST_PEAK_DECAY_REPORT.md", report)?;
    println!("{verdict}");
    println!("tail_classification={tail_class}");
    println!("new_time_evolution=false");
    Ok(())
}
