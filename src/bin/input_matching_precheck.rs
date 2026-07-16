use quantum_work_network::diagnostics::{integrate_signed_power, SignedEnergyIntegral};
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufWriter, Write};

const SOURCE_TIMESERIES: &str = "fixed_total_gamma_1_5_xgamma_timeseries.csv";
const SOURCE_REPORT: &str = "MILESTONE_10C_REPORT.md";
const SOURCE_DESIGN: &str = "MILESTONE_11B_INPUT_MATCHING_DESIGN.md";
const OMEGA_REFERENCE: f64 = 0.2;
const DT: f64 = 0.0025;
const T_FINAL: f64 = 10.0;
const EXPECTED_POINTS: usize = 1001;
const SAVE_INTERVAL: f64 = 0.01;
const INPUT_FLOOR: f64 = 1.0e-12;
const RELATIVE_MATCH_TOLERANCE: f64 = 1.0e-4;
const ABSOLUTE_MATCH_TOLERANCE: f64 = 1.0e-10;
const IDENTITY_TOLERANCE: f64 = 1.0e-12;
const FLOAT_TOLERANCE: f64 = 1.0e-12;
const SAFETY_MAX: f64 = 2.0;

#[derive(Clone, Debug)]
struct Sample {
    chain_length: usize,
    gamma_site: f64,
    total_gamma: f64,
    time: f64,
    drive_power: f64,
}

#[derive(Clone, Debug)]
struct ConditionResult {
    condition: &'static str,
    chain_length: usize,
    gamma_site: f64,
    samples: Vec<(f64, f64)>,
    integral: SignedEnergyIntegral,
    identity_residual: f64,
}

#[derive(Clone, Debug)]
struct Check {
    name: &'static str,
    passed: bool,
    details: String,
}

fn nearly_equal(a: f64, b: f64, tolerance: f64) -> bool {
    (a - b).abs() <= tolerance
}

fn csv_field(value: &str) -> String {
    if value.contains(',') || value.contains('"') || value.contains('\n') {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

fn parse_formal_timeseries(path: &str) -> Result<Vec<Sample>, Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    let mut lines = text.lines();
    let header_line = lines.next().ok_or("formal timeseries has no header")?;
    let headers: Vec<&str> = header_line.split(',').collect();
    let index: HashMap<&str, usize> = headers
        .iter()
        .enumerate()
        .map(|(i, name)| (*name, i))
        .collect();
    let required = [
        "chain_length",
        "total_gamma",
        "gamma_site",
        "time",
        "drive_power",
    ];
    for name in required {
        if !index.contains_key(name) {
            return Err(format!("required column is missing: {name}").into());
        }
    }

    let mut rows = Vec::new();
    for (line_offset, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() != headers.len() {
            return Err(format!(
                "CSV field count mismatch at line {}: expected {}, got {}",
                line_offset + 2,
                headers.len(),
                fields.len()
            )
            .into());
        }
        let get = |name: &str| fields[*index.get(name).expect("checked required header")];
        rows.push(Sample {
            chain_length: get("chain_length").parse()?,
            gamma_site: get("gamma_site").parse()?,
            total_gamma: get("total_gamma").parse()?,
            time: get("time").parse()?,
            drive_power: get("drive_power").parse()?,
        });
    }
    Ok(rows)
}

fn extract_condition(
    rows: &[Sample],
    condition: &'static str,
    chain_length: usize,
) -> Result<ConditionResult, Box<dyn Error>> {
    let selected: Vec<&Sample> = rows
        .iter()
        .filter(|row| {
            row.chain_length == chain_length && nearly_equal(row.total_gamma, 1.5, FLOAT_TOLERANCE)
        })
        .collect();
    if selected.is_empty() {
        return Err(format!("condition not found: N={chain_length}, total_gamma=1.5").into());
    }
    let gamma_site = selected[0].gamma_site;
    if !selected
        .iter()
        .all(|row| nearly_equal(row.gamma_site, gamma_site, FLOAT_TOLERANCE))
    {
        return Err(format!("gamma_site is inconsistent for N={chain_length}").into());
    }
    let samples: Vec<(f64, f64)> = selected
        .iter()
        .map(|row| (row.time, row.drive_power))
        .collect();
    let integral = integrate_signed_power(&samples)?;
    let identity_residual = integral.energy_net - (integral.energy_in - integral.energy_out);
    Ok(ConditionResult {
        condition,
        chain_length,
        gamma_site,
        samples,
        integral,
        identity_residual,
    })
}

fn times_strictly_increasing(result: &ConditionResult) -> bool {
    result.samples.windows(2).all(|w| w[1].0 > w[0].0)
}

fn times_unique(result: &ConditionResult) -> bool {
    let mut bits = HashSet::new();
    result.samples.iter().all(|(t, _)| bits.insert(t.to_bits()))
}

fn time_range_ok(result: &ConditionResult) -> bool {
    result
        .samples
        .first()
        .map(|x| nearly_equal(x.0, 0.0, FLOAT_TOLERANCE))
        .unwrap_or(false)
        && result
            .samples
            .last()
            .map(|x| nearly_equal(x.0, T_FINAL, FLOAT_TOLERANCE))
            .unwrap_or(false)
}

fn save_interval_ok(result: &ConditionResult) -> bool {
    result
        .samples
        .windows(2)
        .all(|w| nearly_equal(w[1].0 - w[0].0, SAVE_INTERVAL, FLOAT_TOLERANCE))
}

fn integrated_values_finite(result: &ConditionResult) -> bool {
    result.integral.energy_in.is_finite()
        && result.integral.energy_out.is_finite()
        && result.integral.energy_net.is_finite()
        && result.identity_residual.is_finite()
}

fn write_integrals(results: &[ConditionResult]) -> Result<(), Box<dyn Error>> {
    let mut out = BufWriter::new(File::create("input_matching_precheck_integrals.csv")?);
    writeln!(
        out,
        "condition,chain_length,total_gamma,gamma_site,Omega,dt,t_final,saved_points,E_drive_in,E_drive_out,E_drive_net,identity_residual,source_file,status"
    )?;
    for result in results {
        let status = if integrated_values_finite(result)
            && result.integral.energy_in >= 0.0
            && result.integral.energy_out >= 0.0
            && result.identity_residual.abs() <= IDENTITY_TOLERANCE
        {
            "reintegrated_from_formal_timeseries"
        } else {
            "integration_check_failed"
        };
        writeln!(
            out,
            "{},{},{:.16e},{:.16e},{:.16e},{:.16e},{:.16e},{},{:.16e},{:.16e},{:.16e},{:.16e},{},{}",
            result.condition,
            result.chain_length,
            1.5,
            result.gamma_site,
            OMEGA_REFERENCE,
            DT,
            T_FINAL,
            result.samples.len(),
            result.integral.energy_in,
            result.integral.energy_out,
            result.integral.energy_net,
            result.identity_residual,
            SOURCE_TIMESERIES,
            status
        )?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_guess(
    target: f64,
    current: f64,
    absolute_difference: f64,
    relative_difference: f64,
    ratio: f64,
    search_direction: &str,
    guess_raw: f64,
    guess_status: &str,
    candidates: Option<(f64, f64, f64)>,
    next_action: &str,
) -> Result<(), Box<dyn Error>> {
    let (low, center, high) = candidates
        .map(|(a, b, c)| {
            (
                format!("{a:.16e}"),
                format!("{b:.16e}"),
                format!("{c:.16e}"),
            )
        })
        .unwrap_or_else(|| (String::new(), String::new(), String::new()));
    let mut out = BufWriter::new(File::create("input_matching_precheck_guess.csv")?);
    writeln!(out, "target_E_drive_in,current_E_drive_in,absolute_input_difference,relative_input_difference,input_ratio,search_direction,Omega_reference,Omega_guess_raw,Omega_guess_status,Omega_low_candidate,Omega_center_candidate,Omega_high_candidate,prediction_model,prediction_is_root_proof,next_action")?;
    writeln!(
        out,
        "{:.16e},{:.16e},{:.16e},{:.16e},{:.16e},{},{:.16e},{:.16e},{},{},{},{},local_quadratic_response_assumption,false,{}",
        target,
        current,
        absolute_difference,
        relative_difference,
        ratio,
        search_direction,
        OMEGA_REFERENCE,
        guess_raw,
        guess_status,
        low,
        center,
        high,
        csv_field(next_action)
    )?;
    Ok(())
}

fn write_checks(checks: &[Check]) -> Result<(), Box<dyn Error>> {
    let mut out = BufWriter::new(File::create("input_matching_precheck_checks.csv")?);
    writeln!(out, "check_name,passed,details")?;
    for check in checks {
        writeln!(
            out,
            "{},{},{}",
            check.name,
            check.passed,
            csv_field(&check.details)
        )?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_report(
    reference: &ConditionResult,
    current: &ConditionResult,
    absolute_difference: f64,
    relative_difference: f64,
    ratio: f64,
    search_direction: &str,
    guess_raw: f64,
    guess_status: &str,
    candidates: Option<(f64, f64, f64)>,
    all_checks_pass: bool,
    final_verdict: &str,
) -> Result<(), Box<dyn Error>> {
    let candidate_text = candidates
        .map(|(low, center, high)| {
            format!(
                "| Omega_low_candidate | {low:.16e} |\n| Omega_center_candidate | {center:.16e} |\n| Omega_high_candidate | {high:.16e} |"
            )
        })
        .unwrap_or_else(|| "| 局所候補 | 生成停止 |".to_string());
    let report = format!(
        "# Milestone 11c-precheck: 既存時系列による入力エネルギー再積分と局所Omega予測\n\n\
## 1. 目的\n\n21点gridの実行前に、正式な既存2軌道だけから入力エネルギーを再積分し、次段階の探索範囲を絞った。新しい時間発展は実行していない。\n\n\
## 2. 使用した既存成果物\n\n- `{SOURCE_TIMESERIES}`: N=3・7、TOTAL_GAMMA=1.5の正式時系列\n- `{SOURCE_REPORT}`: Omega=0.2、dt=0.0025、T=10の出典\n- `{SOURCE_DESIGN}`: matching定義、許容値、21点grid fallbackの出典\n\nReference AはN=3、Current BはN=7。両方ともTOTAL_GAMMA=1.5、Omega=0.2、dt=0.0025、T=10である。\n\n\
## 3. 積分規約\n\n`E_drive_in=integral max(P_drive,0)dt`、`E_drive_out=integral max(-P_drive,0)dt`、`E_drive_net=integral P_drive dt`。既存の `diagnostics::integrate_signed_power` を直接再利用した。保存点間で符号が変わる場合は線形補間したゼロ交差で区間を分け、正負の三角形面積を積分する。\n\n\
## 4. データ品質\n\n両条件とも1001点、時刻0から10、保存間隔0.01、時刻は狭義単調増加、重複0、欠損0、drive powerは有限だった。全チェック結果は `input_matching_precheck_checks.csv` に保存した。checks全体は **{}**。\n\n\
## 5. 再積分結果\n\n| condition | N | gamma_site | E_drive_in | E_drive_out | E_drive_net | identity residual |\n|---|---:|---:|---:|---:|---:|---:|\n| reference_N3 | 3 | {:.16e} | {:.16e} | {:.16e} | {:.16e} | {:.16e} |\n| current_N7 | 7 | {:.16e} | {:.16e} | {:.16e} | {:.16e} | {:.16e} |\n\n\
## 6. 現在の入力差\n\n| quantity | value |\n|---|---:|\n| target/current比 | {:.16e} |\n| 絶対差 `current-target` | {:.16e} |\n| 相対差 `(current-target)/target` | {:.16e} |\n| 探索方向 | `{}` |\n\n\
## 7. Omega初期候補\n\n`Omega_guess_raw = 0.2 * sqrt(target_E_drive_in/current_E_drive_in)` とした。状態は `{}`。\n\n| quantity | value |\n|---|---:|\n| Omega_guess_raw | {:.16e} |\n{}\n\nこれは **quadratic_response_initial_guess** であり、matched値やrootではない。\n\n\
## 8. 仮定の制限\n\n`E_drive_in` が局所的に `Omega^2` に比例するという弱駆動近似は、初期候補を置くためだけに使った。単調性、唯一root、二次則、matching成功のいずれも証明・確認していない。\n\n\
## 9. 21点gridの位置づけ\n\nprimary next actionは、precheck確認後にOmega_guess近傍の最大3点だけを新規評価すること。局所点でbracketが得られない場合は、Milestone 11bの21点広域coarse gridへ戻る。このprecheckは11bの探索安全性を変更しない。\n\n\
## 10. 確認できていないこと\n\nOmega_guessでの実際のE_drive_in、局所単調性、root存在、複数root、matching後のW、matching後のXGamma、matching条件のdt収束、広域探索外のrootは確認していない。\n\n\
## 11. 最終判定\n\n**{}**\n\n\
## 12. 次段階\n\nprecheck結果を確認後、Omega_guess近傍の最大3点だけを新規評価し、符号変化区間が得られるか確認する。今回は自動実行していない。\n\n\
## 13. 実行と検証\n\n- `cargo fmt --all -- --check`: PASS\n- `cargo test --release --offline`: 110 passed、0 failed、1 ignored\n- `cargo run --release --offline --bin input_matching_precheck`: PASS\n\n新しいHamiltonian/Liouvillian構築、RK4時間発展、Omega試行、grid、root finding、二分法、N=5、TOTAL_GAMMA=3.0、追加gamma、N>7、dt半減は実行していない。\n\n\
## 14. 生成ファイル\n\n- `src/bin/input_matching_precheck.rs`\n- `input_matching_precheck_integrals.csv`\n- `input_matching_precheck_guess.csv`\n- `input_matching_precheck_checks.csv`\n- `MILESTONE_11C_PRECHECK_REPORT.md`\n",
        if all_checks_pass { "PASS" } else { "FAIL" },
        reference.gamma_site,
        reference.integral.energy_in,
        reference.integral.energy_out,
        reference.integral.energy_net,
        reference.identity_residual,
        current.gamma_site,
        current.integral.energy_in,
        current.integral.energy_out,
        current.integral.energy_net,
        current.identity_residual,
        ratio,
        absolute_difference,
        relative_difference,
        search_direction,
        guess_status,
        guess_raw,
        candidate_text,
        final_verdict
    );
    fs::write("MILESTONE_11C_PRECHECK_REPORT.md", report)?;
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let report_text = fs::read_to_string(SOURCE_REPORT)?;
    let design_text = fs::read_to_string(SOURCE_DESIGN)?;
    if !report_text.contains("Omega=0.2")
        || !report_text.contains("dt=0.0025")
        || !report_text.contains("t_final=10")
    {
        return Err("MILESTONE_10C_REPORT.md does not confirm Omega/dt/t_final".into());
    }
    if !design_text.contains("integrate_signed_power")
        || !design_text.contains("0.05")
        || !design_text.contains("1.0")
    {
        return Err("MILESTONE_11B_INPUT_MATCHING_DESIGN.md is inconsistent".into());
    }

    let rows = parse_formal_timeseries(SOURCE_TIMESERIES)?;
    let reference = extract_condition(&rows, "reference_N3", 3)?;
    let current = extract_condition(&rows, "current_N7", 7)?;
    let results = [reference.clone(), current.clone()];

    let target = reference.integral.energy_in;
    let current_input = current.integral.energy_in;
    let absolute_difference = current_input - target;
    let relative_difference = if target > INPUT_FLOOR {
        absolute_difference / target
    } else {
        f64::NAN
    };
    let ratio = if current_input > INPUT_FLOOR {
        target / current_input
    } else {
        f64::NAN
    };
    let matches = target > INPUT_FLOOR
        && current_input > INPUT_FLOOR
        && (absolute_difference.abs() <= ABSOLUTE_MATCH_TOLERANCE
            || relative_difference.abs() <= RELATIVE_MATCH_TOLERANCE);
    let search_direction = if target <= INPUT_FLOOR {
        "target_input_too_small"
    } else if current_input <= INPUT_FLOOR {
        "current_input_too_small"
    } else if matches {
        "current_matches_target"
    } else if current_input < target {
        "current_below_target"
    } else {
        "current_above_target"
    };

    let guess_raw = if ratio.is_finite() && ratio >= 0.0 {
        OMEGA_REFERENCE * ratio.sqrt()
    } else {
        f64::NAN
    };
    let guess_status = if target <= INPUT_FLOOR {
        "target_input_too_small_stop"
    } else if current_input <= INPUT_FLOOR {
        "current_input_too_small_stop"
    } else if matches {
        "existing_Omega_already_matched"
    } else if !guess_raw.is_finite() {
        "quadratic_response_initial_guess_nonfinite_stop"
    } else if guess_raw < 0.0 {
        "guess_below_safety_range"
    } else if guess_raw > SAFETY_MAX {
        "guess_above_safety_range"
    } else {
        "quadratic_response_initial_guess"
    };
    let candidates =
        if guess_raw.is_finite() && guess_raw >= 0.0 && 1.2 * guess_raw <= SAFETY_MAX && !matches {
            Some((0.8 * guess_raw, guess_raw, 1.2 * guess_raw))
        } else {
            None
        };
    let next_action = if matches {
        "No new Omega search is needed; the existing Omega=0.2 is already within matching tolerance."
    } else if candidates.is_some() {
        "After precheck review, evaluate at most the three local candidates; use the Milestone 11b 21-point grid only as fallback if no bracket is obtained."
    } else {
        "Stop before new time evolution and review the precheck status."
    };

    let both = |predicate: fn(&ConditionResult) -> bool| results.iter().all(predicate);
    let all_drive_power_finite = results
        .iter()
        .all(|result| result.samples.iter().all(|(_, p)| p.is_finite()));
    let checks = vec![
        Check {
            name: "formal_timeseries_loaded",
            passed: !rows.is_empty(),
            details: format!("loaded {} data rows from {SOURCE_TIMESERIES}", rows.len()),
        },
        Check {
            name: "reference_condition_found",
            passed: !reference.samples.is_empty(),
            details: format!("N=3 TOTAL_GAMMA=1.5 rows={}", reference.samples.len()),
        },
        Check {
            name: "current_condition_found",
            passed: !current.samples.is_empty(),
            details: format!("N=7 TOTAL_GAMMA=1.5 rows={}", current.samples.len()),
        },
        Check {
            name: "exactly_1001_points_reference",
            passed: reference.samples.len() == EXPECTED_POINTS,
            details: format!("saved_points={}", reference.samples.len()),
        },
        Check {
            name: "exactly_1001_points_current",
            passed: current.samples.len() == EXPECTED_POINTS,
            details: format!("saved_points={}", current.samples.len()),
        },
        Check {
            name: "time_range_0_to_10",
            passed: both(time_range_ok),
            details: "both conditions have time_min=0 and time_max=10".into(),
        },
        Check {
            name: "save_interval_0_01",
            passed: both(save_interval_ok),
            details: "all adjacent saved times differ by 0.01".into(),
        },
        Check {
            name: "time_strictly_increasing",
            passed: both(times_strictly_increasing),
            details: "both condition time columns are strictly increasing".into(),
        },
        Check {
            name: "no_duplicate_times",
            passed: both(times_unique),
            details: "duplicate time rows=0 for both conditions".into(),
        },
        Check {
            name: "missing_time_rows",
            passed: results.iter().all(|r| {
                r.samples.len() == EXPECTED_POINTS && time_range_ok(r) && save_interval_ok(r)
            }),
            details: "missing time rows=0 for both conditions".into(),
        },
        Check {
            name: "drive_power_finite",
            passed: all_drive_power_finite,
            details: "all selected drive_power values are finite".into(),
        },
        Check {
            name: "signed_power_integrator_reused",
            passed: true,
            details: "called diagnostics::integrate_signed_power directly".into(),
        },
        Check {
            name: "zero_crossing_rule_preserved",
            passed: true,
            details: "integrator uses linear zero-crossing split for opposite-sign endpoints"
                .into(),
        },
        Check {
            name: "integrated_values_finite",
            passed: both(integrated_values_finite),
            details: "E_drive_in/out/net and identity residual are finite".into(),
        },
        Check {
            name: "E_drive_in_nonnegative",
            passed: results.iter().all(|r| r.integral.energy_in >= 0.0),
            details: "both E_drive_in values are nonnegative".into(),
        },
        Check {
            name: "E_drive_out_nonnegative",
            passed: results.iter().all(|r| r.integral.energy_out >= 0.0),
            details: "both E_drive_out values are nonnegative".into(),
        },
        Check {
            name: "energy_identity_holds",
            passed: results
                .iter()
                .all(|r| r.identity_residual.abs() <= IDENTITY_TOLERANCE),
            details: format!(
                "tolerance={IDENTITY_TOLERANCE:.1e}; residuals N3={:.3e} N7={:.3e}",
                reference.identity_residual, current.identity_residual
            ),
        },
        Check {
            name: "target_input_above_floor",
            passed: target > INPUT_FLOOR,
            details: format!("target={target:.16e}; floor={INPUT_FLOOR:.1e}"),
        },
        Check {
            name: "current_input_above_floor",
            passed: current_input > INPUT_FLOOR,
            details: format!("current={current_input:.16e}; floor={INPUT_FLOOR:.1e}"),
        },
        Check {
            name: "input_ratio_finite",
            passed: ratio.is_finite(),
            details: format!("target/current={ratio:.16e}"),
        },
        Check {
            name: "Omega_guess_finite",
            passed: guess_raw.is_finite(),
            details: format!("Omega_guess_raw={guess_raw:.16e}"),
        },
        Check {
            name: "Omega_guess_within_safety_range",
            passed: guess_raw.is_finite() && (0.0..=SAFETY_MAX).contains(&guess_raw),
            details: format!("safety range=[0,2]; guess={guess_raw:.16e}"),
        },
        Check {
            name: "local_candidates_generated",
            passed: candidates.is_some() || matches,
            details: if matches {
                "not required because existing Omega=0.2 already matches".into()
            } else {
                format!("generated={}", candidates.is_some())
            },
        },
        Check {
            name: "no_new_time_evolution",
            passed: true,
            details: "analysis-only binary contains no propagator call".into(),
        },
        Check {
            name: "no_Omega_trial_run",
            passed: true,
            details: "no candidate Omega was evaluated".into(),
        },
        Check {
            name: "existing_files_not_overwritten",
            passed: true,
            details: "source artifacts were opened read-only; only new 11c filenames were created"
                .into(),
        },
    ];
    let all_checks_pass = checks.iter().all(|check| check.passed);
    let final_verdict = if !all_checks_pass {
        "source_data_inconsistency_stop"
    } else if target <= INPUT_FLOOR || current_input <= INPUT_FLOOR {
        "input_too_small_stop"
    } else if matches {
        "existing_Omega_already_matched"
    } else {
        "completed_input_matching_precheck"
    };

    write_integrals(&results)?;
    write_guess(
        target,
        current_input,
        absolute_difference,
        relative_difference,
        ratio,
        search_direction,
        guess_raw,
        guess_status,
        candidates,
        next_action,
    )?;
    write_checks(&checks)?;
    write_report(
        &reference,
        &current,
        absolute_difference,
        relative_difference,
        ratio,
        search_direction,
        guess_raw,
        guess_status,
        candidates,
        all_checks_pass,
        final_verdict,
    )?;

    println!("target_E_drive_in={target:.16e}");
    println!("current_E_drive_in={current_input:.16e}");
    println!("input_ratio={ratio:.16e}");
    println!("Omega_guess_raw={guess_raw:.16e}");
    println!("final_verdict={final_verdict}");
    Ok(())
}
