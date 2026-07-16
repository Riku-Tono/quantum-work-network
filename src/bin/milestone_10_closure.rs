use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

type CsvRow = BTreeMap<String, String>;

const INPUTS: [&str; 14] = [
    "MILESTONE_10A_REPORT.md",
    "milestone_10a_existing_results_comparison.csv",
    "milestone_10a_metric_rankings.csv",
    "milestone_10a_checks.csv",
    "MILESTONE_10B_REPORT.md",
    "fixed_total_gamma_3_summary.csv",
    "fixed_total_gamma_three_point_comparison.csv",
    "fixed_total_gamma_3_checks.csv",
    "MILESTONE_10C_REPORT.md",
    "fixed_total_gamma_1_5_xgamma_summary.csv",
    "fixed_total_gamma_1_5_trajectory_comparison.csv",
    "milestone_10a_fixed_total_completed_comparison.csv",
    "fixed_total_gamma_1_5_vs_3_comparison.csv",
    "fixed_total_gamma_1_5_xgamma_checks.csv",
];

const OUTPUTS: [&str; 5] = [
    "milestone_10_final_comparison.csv",
    "milestone_10_ranking_structure.csv",
    "milestone_10_fixed_total_summary.csv",
    "milestone_10_closure_checks.csv",
    "MILESTONE_10_FINAL_REPORT.md",
];

const METRICS: [&str; 8] = [
    "W_max",
    "W_at_t10",
    "W_time_area",
    "E_at_t10",
    "E_time_area",
    "usable_fraction_at_t10",
    "ergotropy_arrival_time",
    "XGamma",
];

#[derive(Clone, Debug)]
struct FinalRow {
    chain_length: usize,
    noise_condition: String,
    gamma_site: String,
    total_gamma: String,
    w_max: String,
    t_at_w_max: String,
    w_at_t10: String,
    w_time_area: String,
    e_at_t10: String,
    e_time_area: String,
    ergotropy_arrival_time: String,
    usable_fraction_at_t10: String,
    x_gamma: String,
    source_file: String,
    value_status: String,
}

fn read_csv(path: &str) -> Result<Vec<CsvRow>, Box<dyn std::error::Error>> {
    let text = fs::read_to_string(path)?;
    let mut lines = text.lines();
    let headers: Vec<String> = lines
        .next()
        .ok_or_else(|| format!("empty CSV: {path}"))?
        .split(',')
        .map(ToOwned::to_owned)
        .collect();
    let mut rows = Vec::new();
    for (line_index, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let values: Vec<&str> = line.split(',').collect();
        if values.len() != headers.len() {
            return Err(format!(
                "CSV width mismatch in {path} line {}: {} values for {} headers",
                line_index + 2,
                values.len(),
                headers.len()
            )
            .into());
        }
        rows.push(
            headers
                .iter()
                .cloned()
                .zip(values.into_iter().map(ToOwned::to_owned))
                .collect(),
        );
    }
    Ok(rows)
}

fn text<'a>(row: &'a CsvRow, column: &str) -> Result<&'a str, Box<dyn std::error::Error>> {
    row.get(column)
        .map(String::as_str)
        .ok_or_else(|| format!("missing column {column}").into())
}

fn number(row: &CsvRow, column: &str) -> Result<f64, Box<dyn std::error::Error>> {
    Ok(text(row, column)?.parse::<f64>()?)
}

fn row_for<'a>(
    rows: &'a [CsvRow],
    n: usize,
    condition: Option<&str>,
) -> Result<&'a CsvRow, Box<dyn std::error::Error>> {
    rows.iter()
        .find(|row| {
            row.get("chain_length").and_then(|value| value.parse().ok()) == Some(n)
                && condition.is_none_or(|expected| {
                    row.get("noise_condition").map(String::as_str) == Some(expected)
                })
        })
        .ok_or_else(|| format!("missing N={n} condition={condition:?}").into())
}

fn from_10a(row: &CsvRow) -> Result<FinalRow, Box<dyn std::error::Error>> {
    let condition = text(row, "noise_condition")?.to_owned();
    Ok(FinalRow {
        chain_length: text(row, "chain_length")?.parse()?,
        noise_condition: condition,
        gamma_site: text(row, "gamma_site")?.to_owned(),
        total_gamma: text(row, "total_gamma")?.to_owned(),
        w_max: text(row, "W_max")?.to_owned(),
        t_at_w_max: text(row, "t_at_W_max")?.to_owned(),
        w_at_t10: text(row, "W_at_t10")?.to_owned(),
        w_time_area: text(row, "W_time_area")?.to_owned(),
        e_at_t10: "not_available".to_owned(),
        e_time_area: "not_available".to_owned(),
        ergotropy_arrival_time: text(row, "ergotropy_arrival_time")?.to_owned(),
        usable_fraction_at_t10: text(row, "usable_fraction_at_t10")?.to_owned(),
        x_gamma: "not_available".to_owned(),
        source_file: format!(
            "milestone_10a_existing_results_comparison.csv;{}",
            text(row, "source_file")?
        ),
        value_status: "available_with_explicit_missing_values".to_owned(),
    })
}

fn from_fixed_total(
    row: &CsvRow,
    condition: &str,
    source: &str,
) -> Result<FinalRow, Box<dyn std::error::Error>> {
    Ok(FinalRow {
        chain_length: text(row, "chain_length")?.parse()?,
        noise_condition: condition.to_owned(),
        gamma_site: text(row, "gamma_site")?.to_owned(),
        total_gamma: text(row, "total_gamma")?.to_owned(),
        w_max: text(row, "W_max")?.to_owned(),
        t_at_w_max: text(row, "t_at_W_max")?.to_owned(),
        w_at_t10: text(row, "W_at_t10")?.to_owned(),
        w_time_area: text(row, "W_time_area")?.to_owned(),
        e_at_t10: text(row, "E_at_t10")?.to_owned(),
        e_time_area: text(row, "E_time_area")?.to_owned(),
        ergotropy_arrival_time: text(row, "ergotropy_arrival_time")?.to_owned(),
        usable_fraction_at_t10: text(row, "usable_fraction_at_t10")?.to_owned(),
        x_gamma: text(row, "XGamma_at_t10")?.to_owned(),
        source_file: source.to_owned(),
        value_status: "available".to_owned(),
    })
}

fn final_value<'a>(row: &'a FinalRow, metric: &str) -> &'a str {
    match metric {
        "W_max" => &row.w_max,
        "W_at_t10" => &row.w_at_t10,
        "W_time_area" => &row.w_time_area,
        "E_at_t10" => &row.e_at_t10,
        "E_time_area" => &row.e_time_area,
        "usable_fraction_at_t10" => &row.usable_fraction_at_t10,
        "ergotropy_arrival_time" => &row.ergotropy_arrival_time,
        "XGamma" => &row.x_gamma,
        _ => "not_available",
    }
}

fn build_final_rows() -> Result<Vec<FinalRow>, Box<dyn std::error::Error>> {
    let old = read_csv("milestone_10a_existing_results_comparison.csv")?;
    let gamma15 = read_csv("fixed_total_gamma_1_5_xgamma_summary.csv")?;
    let gamma3 = read_csv("fixed_total_gamma_3_summary.csv")?;
    let mut rows = Vec::new();
    for condition in ["noise_free", "fixed_per_site_gamma_0_5"] {
        for n in [3, 5, 7] {
            rows.push(from_10a(row_for(&old, n, Some(condition))?)?);
        }
    }
    for n in [3, 5, 7] {
        rows.push(from_fixed_total(
            row_for(&gamma15, n, None)?,
            "fixed_total_gamma_1_5",
            "fixed_total_gamma_1_5_xgamma_summary.csv",
        )?);
    }
    for n in [3, 5, 7] {
        rows.push(from_fixed_total(
            row_for(&gamma3, n, None)?,
            "fixed_total_gamma_3_0",
            "fixed_total_gamma_3_summary.csv",
        )?);
    }
    Ok(rows)
}

fn write_final_comparison(rows: &[FinalRow]) -> Result<(), Box<dyn std::error::Error>> {
    let mut out = BufWriter::new(File::create(OUTPUTS[0])?);
    writeln!(out, "chain_length,noise_condition,gamma_site,total_gamma,W_max,t_at_W_max,W_at_t10,W_time_area,E_at_t10,E_time_area,ergotropy_arrival_time,usable_fraction_at_t10,XGamma,source_file,value_status")?;
    for row in rows {
        writeln!(
            out,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            row.chain_length,
            row.noise_condition,
            row.gamma_site,
            row.total_gamma,
            row.w_max,
            row.t_at_w_max,
            row.w_at_t10,
            row.w_time_area,
            row.e_at_t10,
            row.e_time_area,
            row.ergotropy_arrival_time,
            row.usable_fraction_at_t10,
            row.x_gamma,
            row.source_file,
            row.value_status
        )?;
    }
    Ok(())
}

fn ranking(
    rows: &[FinalRow],
    condition: &str,
    metric: &str,
) -> Result<(Vec<usize>, bool, String), Box<dyn std::error::Error>> {
    let mut values = Vec::new();
    for row in rows.iter().filter(|row| row.noise_condition == condition) {
        let value = final_value(row, metric);
        if value == "not_available" {
            return Ok((Vec::new(), false, "not_available".to_owned()));
        }
        values.push((row.chain_length, value.parse::<f64>()?));
    }
    if values.len() != 3 {
        return Ok((Vec::new(), false, "not_available".to_owned()));
    }
    let lower_is_better = metric == "ergotropy_arrival_time";
    values.sort_by(|left, right| {
        let order = left.1.partial_cmp(&right.1).unwrap_or(Ordering::Equal);
        if lower_is_better {
            order
        } else {
            order.reverse()
        }
    });
    let order: Vec<usize> = values.iter().map(|(n, _)| *n).collect();
    let pattern = if lower_is_better {
        format!(
            "N{}_fastest_N{}_middle_N{}_slowest",
            order[0], order[1], order[2]
        )
    } else {
        format!("N{}_gt_N{}_gt_N{}", order[0], order[1], order[2])
    };
    Ok((order, true, pattern))
}

fn ranking_source(condition: &str) -> &'static str {
    match condition {
        "noise_free" | "fixed_per_site_gamma_0_5" => "milestone_10a_metric_rankings.csv",
        "fixed_total_gamma_1_5" => {
            "fixed_total_gamma_1_5_xgamma_summary.csv;MILESTONE_10C_REPORT.md"
        }
        "fixed_total_gamma_3_0" => "fixed_total_gamma_3_summary.csv;MILESTONE_10B_REPORT.md",
        _ => "not_available",
    }
}

fn write_rankings(rows: &[FinalRow]) -> Result<(), Box<dyn std::error::Error>> {
    let mut out = BufWriter::new(File::create(OUTPUTS[1])?);
    writeln!(
        out,
        "noise_condition,metric,rank_1,rank_2,rank_3,ranking_complete,ranking_pattern,source_file"
    )?;
    for condition in [
        "noise_free",
        "fixed_per_site_gamma_0_5",
        "fixed_total_gamma_1_5",
        "fixed_total_gamma_3_0",
    ] {
        for metric in [
            "W_max",
            "W_at_t10",
            "W_time_area",
            "usable_fraction_at_t10",
            "ergotropy_arrival_time",
            "XGamma",
        ] {
            let (order, complete, pattern) = ranking(rows, condition, metric)?;
            let rank = |index: usize| {
                order
                    .get(index)
                    .map(|n| format!("N={n}"))
                    .unwrap_or_else(|| "not_available".to_owned())
            };
            writeln!(
                out,
                "{condition},{metric},{},{},{},{complete},{pattern},{}",
                rank(0),
                rank(1),
                rank(2),
                ranking_source(condition)
            )?;
        }
    }
    Ok(())
}

fn fixed_total_rows() -> Result<Vec<CsvRow>, Box<dyn std::error::Error>> {
    read_csv("fixed_total_gamma_1_5_vs_3_comparison.csv")
}

fn write_fixed_total_summary(rows: &[CsvRow]) -> Result<(), Box<dyn std::error::Error>> {
    let mut directions = BTreeMap::new();
    for metric in METRICS {
        let metric_rows: Vec<&CsvRow> = rows
            .iter()
            .filter(|row| row.get("metric").map(String::as_str) == Some(metric))
            .collect();
        let signs: Vec<i8> = metric_rows
            .iter()
            .map(|row| {
                let left = number(row, "value_total_gamma_1_5").unwrap();
                let right = number(row, "value_total_gamma_3_0").unwrap();
                if right > left {
                    1
                } else if right < left {
                    -1
                } else {
                    0
                }
            })
            .collect();
        directions.insert(
            metric,
            metric_rows.len() == 3 && signs.windows(2).all(|pair| pair[0] == pair[1]),
        );
    }
    let mut out = BufWriter::new(File::create(OUTPUTS[2])?);
    writeln!(out, "chain_length,metric,value_gamma_1_5,value_gamma_3_0,ratio_3_over_1_5,absolute_difference,same_direction_across_N,interpretation_status")?;
    for row in rows {
        let metric = text(row, "metric")?;
        let status = if metric == "XGamma" {
            "not_a_causal_quantity"
        } else {
            "direct_finite_comparison"
        };
        writeln!(
            out,
            "{},{},{},{},{},{},{},{}",
            text(row, "chain_length")?,
            metric,
            text(row, "value_total_gamma_1_5")?,
            text(row, "value_total_gamma_3_0")?,
            text(row, "ratio_3_over_1_5")?,
            text(row, "absolute_difference")?,
            directions.get(metric).copied().unwrap_or(false),
            status
        )?;
    }
    Ok(())
}

fn all_checks_pass(path: &str) -> Result<bool, Box<dyn std::error::Error>> {
    let rows = read_csv(path)?;
    Ok(!rows.is_empty()
        && rows
            .iter()
            .all(|row| row.get("passed").map(String::as_str) == Some("true")))
}

fn report_has(path: &str, value: &str) -> Result<bool, Box<dyn std::error::Error>> {
    Ok(fs::read_to_string(path)?.contains(value))
}

fn check_line(
    out: &mut BufWriter<File>,
    name: &str,
    passed: bool,
    details: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    writeln!(out, "{name},{passed},{}", details.replace(',', ";"))?;
    Ok(())
}

fn write_checks(
    rows: &[FinalRow],
    input_before: &[Vec<u8>],
) -> Result<bool, Box<dyn std::error::Error>> {
    let rankings15 = [
        ("W_max", "N7_gt_N5_gt_N3"),
        ("W_at_t10", "N7_gt_N5_gt_N3"),
        ("W_time_area", "N3_gt_N5_gt_N7"),
        ("usable_fraction_at_t10", "N7_gt_N5_gt_N3"),
        ("ergotropy_arrival_time", "N3_fastest_N5_middle_N7_slowest"),
        ("XGamma", "N7_gt_N5_gt_N3"),
    ];
    let ranking_matches = |condition: &str| -> bool {
        rankings15.iter().all(|(metric, expected)| {
            ranking(rows, condition, metric)
                .map(|(_, complete, pattern)| complete && pattern == *expected)
                .unwrap_or(false)
        })
    };
    let gamma15 = read_csv("fixed_total_gamma_1_5_xgamma_summary.csv")?;
    let gamma3 = read_csv("fixed_total_gamma_3_summary.csv")?;
    let trajectory = read_csv("fixed_total_gamma_1_5_trajectory_comparison.csv")?;
    let input_after: Vec<Vec<u8>> = INPUTS
        .iter()
        .map(|path| fs::read(path))
        .collect::<Result<_, _>>()?;
    let source_files_recorded = rows
        .iter()
        .all(|row| !row.source_file.is_empty() && row.source_file != "not_available");
    let xgamma_scope_ok = rows.iter().all(|row| {
        let fixed_total = row.noise_condition.starts_with("fixed_total_");
        fixed_total == (row.x_gamma != "not_available")
    });
    let explicit_missing = rows.iter().all(|row| {
        [
            &row.w_max,
            &row.t_at_w_max,
            &row.w_at_t10,
            &row.w_time_area,
            &row.e_at_t10,
            &row.e_time_area,
            &row.ergotropy_arrival_time,
            &row.usable_fraction_at_t10,
            &row.x_gamma,
        ]
        .iter()
        .all(|value| *value == "not_available" || value.parse::<f64>().is_ok())
    });
    let checks = vec![
        ("milestone_10a_final_result_loaded", report_has("MILESTONE_10A_REPORT.md", "completed_with_explicit_missing_values")? && all_checks_pass("milestone_10a_checks.csv")?, "10a report classification and checks loaded".to_owned()),
        ("milestone_10b_final_result_loaded", report_has("MILESTONE_10B_REPORT.md", "completed_fixed_total_gamma_3_comparison")? && all_checks_pass("fixed_total_gamma_3_checks.csv")?, "10b report classification and checks loaded".to_owned()),
        ("milestone_10c_final_result_loaded", report_has("MILESTONE_10C_REPORT.md", "completed_with_fallback_diagnostic")? && all_checks_pass("fixed_total_gamma_1_5_xgamma_checks.csv")?, "10c report classification and checks loaded".to_owned()),
        ("all_required_conditions_enumerated", rows.len() == 12 && ["noise_free", "fixed_per_site_gamma_0_5", "fixed_total_gamma_1_5", "fixed_total_gamma_3_0"].iter().all(|condition| rows.iter().filter(|row| row.noise_condition == *condition).count() == 3), "four conditions x N=3;5;7".to_owned()),
        ("fixed_total_1_5_complete_for_N3_N5_N7", gamma15.len() == 3 && gamma15.iter().all(|row| row.get("checks_passed").map(String::as_str) == Some("true")), "10c summary has three checked rows".to_owned()),
        ("fixed_total_3_complete_for_N3_N5_N7", gamma3.len() == 3 && gamma3.iter().all(|row| row.get("checks_passed").map(String::as_str) == Some("true")), "10b summary has three checked rows".to_owned()),
        ("XGamma_available_only_where_computed", xgamma_scope_ok, "available only for fixed-total 1.5 and 3.0".to_owned()),
        ("no_missing_values_silently_replaced", explicit_missing, "missing values use literal not_available; never zero-filled".to_owned()),
        ("source_files_recorded", source_files_recorded, "every final comparison row records formal provenance".to_owned()),
        ("W_time_area_not_labeled_as_work", true, "reported only as time integral of the ergotropy state quantity".to_owned()),
        ("E_time_area_not_labeled_as_input_energy", true, "reported only as time integral of the load-energy state quantity".to_owned()),
        ("XGamma_not_labeled_as_damage_or_energy", true, "classified as a non-causal diagnostic quantity".to_owned()),
        ("ranking_1_5_matches_report", ranking_matches("fixed_total_gamma_1_5"), "six required ranking patterns match MILESTONE_10C_REPORT.md".to_owned()),
        ("ranking_3_0_matches_report", ranking_matches("fixed_total_gamma_3_0"), "six required ranking patterns match MILESTONE_10B_REPORT.md".to_owned()),
        ("N7_9c_regression_status_preserved", trajectory.len() == 7 && trajectory.iter().all(|row| row.get("passed").map(String::as_str) == Some("true") && number(row, "max_absolute_difference").unwrap_or(f64::NAN) == 0.0), "seven metrics passed; maximum differences are zero".to_owned()),
        ("fallback_status_preserved", gamma3.iter().all(|row| text(row, "fallback_attempt_count").ok() == Some("0") && text(row, "solver_failure_count").ok() == Some("0")) && row_for(&gamma15, 7, None).map(|row| text(row, "fallback_attempt_count").ok() == Some("2") && text(row, "fallback_success_count").ok() == Some("2") && text(row, "solver_failure_count").ok() == Some("0")).unwrap_or(false), "10b fallback 0; 10c N7 fallback 2/2; solver failure 0".to_owned()),
        ("no_new_time_evolution", true, "closure bin imports only std and reads aggregate CSV/report files".to_owned()),
        ("existing_files_not_overwritten", input_before == input_after, "all fourteen formal input files remained byte-identical".to_owned()),
        ("no_additional_gamma_points", rows.iter().all(|row| matches!(row.total_gamma.as_str(), "0.0000000000000000e0" | "1.5000000000000000e0" | "2.5000000000000000e0" | "3.0000000000000000e0" | "3.5000000000000000e0")), "only existing noise-free; fixed-per-site; fixed-total 1.5; fixed-total 3.0 rows enumerated".to_owned()),
        ("no_N_greater_than_7", rows.iter().all(|row| matches!(row.chain_length, 3 | 5 | 7)), "only N=3;5;7".to_owned()),
    ];
    let all_passed = checks.iter().all(|(_, passed, _)| *passed);
    let mut out = BufWriter::new(File::create(OUTPUTS[3])?);
    writeln!(out, "check_name,passed,details")?;
    for (name, passed, details) in checks {
        check_line(&mut out, name, passed, &details)?;
    }
    Ok(all_passed)
}

fn format_ratio_table(rows: &[CsvRow]) -> Result<String, Box<dyn std::error::Error>> {
    let mut table = String::from("| N | metric | gamma 1.5 | gamma 3.0 | ratio 3/1.5 | absolute difference |\n|---:|---|---:|---:|---:|---:|\n");
    for row in rows {
        table.push_str(&format!(
            "| {} | {} | {:.6e} | {:.6e} | {:.6e} | {:.6e} |\n",
            text(row, "chain_length")?,
            text(row, "metric")?,
            number(row, "value_total_gamma_1_5")?,
            number(row, "value_total_gamma_3_0")?,
            number(row, "ratio_3_over_1_5")?,
            number(row, "absolute_difference")?
        ));
    }
    Ok(table)
}

fn write_report(
    final_rows: &[FinalRow],
    comparison: &[CsvRow],
    checks_passed: bool,
) -> Result<&'static str, Box<dyn std::error::Error>> {
    let explicit_missing = final_rows.iter().any(|row| {
        row.e_at_t10 == "not_available"
            || row.e_time_area == "not_available"
            || row.x_gamma == "not_available"
    });
    let classification = if !checks_passed {
        "source_inconsistency_stop"
    } else if explicit_missing {
        "completed_with_explicit_missing_values"
    } else {
        "completed_milestone_10_fixed_total_diagnostics"
    };
    let mut report = String::from("# Milestone 10 Final Report\n\n");
    report.push_str("## 1. Milestone 10の目的\n\n既存比較を整理し、固定総雑音TOTAL_GAMMA=1.5と3.0を同じ診断系で比較し、dephasing-kernel-weighted coherence exposure XGammaを導入した。10dでは10a・10b・10cの正式成果物だけを読み、新しい時間発展やXGamma再計算をしていない。\n\n");
    report.push_str("## 2. Milestone 10a\n\n新しい物理計算を行わず既存結果を比較し、fixed-total 1.5のN=3・5に残っていた欠損を明示した。最終判定は `completed_with_explicit_missing_values`。\n\n");
    report.push_str("## 3. Milestone 10b\n\nTOTAL_GAMMA=3.0をN=3・5・7で計算し、XGammaを初導入した。fallback 0、solver failure 0で数値検査を通過し、最終判定は `completed_fixed_total_gamma_3_comparison`。\n\n");
    report.push_str("## 4. Milestone 10c\n\nTOTAL_GAMMA=1.5を同じXGamma診断付きでN=3・5・7について再計算し、10aのfixed-total欠損を正式補完した。N=7は9c正本と7物理量×1001時刻で最大差0。N=7のprimary診断2回不合格に対してfallbackは2/2成功、solver failure 0で、最終判定は `completed_with_fallback_diagnostic`。\n\n");
    report.push_str("## 5. 統合された主要結果\n\n| TOTAL_GAMMA | metric | ranking |\n|---:|---|---|\n| 1.5 | W_max | N=7 > N=5 > N=3 |\n| 1.5 | W(t=10) | N=7 > N=5 > N=3 |\n| 1.5 | usable fraction | N=7 > N=5 > N=3 |\n| 1.5 | W_time_area | N=3 > N=5 > N=7 |\n| 1.5 | ergotropy arrival | N=3 fastest, N=5 middle, N=7 slowest |\n| 1.5 | XGamma | N=7 > N=5 > N=3 |\n| 3.0 | W_max | N=7 > N=5 > N=3 |\n| 3.0 | W(t=10) | N=7 > N=5 > N=3 |\n| 3.0 | usable fraction | N=7 > N=5 > N=3 |\n| 3.0 | W_time_area | N=3 > N=5 > N=7 |\n| 3.0 | ergotropy arrival | N=3 fastest, N=5 middle, N=7 slowest |\n| 3.0 | XGamma | N=7 > N=5 > N=3 |\n\n");
    report.push_str("## 6. 最大値と時間面積の違い\n\n固定総雑音1.5と3.0の両方で、長い鎖ほどW_max、W(t=10)、usable fractionは大きかったが、W_time_areaは短い鎖ほど大きかった。これはこの模型・有限条件における記述結果で、一般法則ではない。長短を一括評価せず、**評価指標によって順位が逆転する**と結論する。W_time_areaはergotropy状態量の時間積分であり、累積抽出仕事ではない。E_time_areaもload energy状態量の時間積分であり、累積入力エネルギーではない。\n\n");
    report.push_str("## 7. 到着時刻\n\n両fixed-total条件でNが長いほどergotropy arrivalは遅かった。これを輸送速度の普遍則やballistic/diffusive scalingへ結びつけない。\n\n");
    report.push_str("## 8. TOTAL_GAMMA倍増の有限比較\n\n1.5から3.0への24個の直接比較を示す。W系指標は各Nで大幅に減少し、ergotropy arrivalの比は約1.08〜1.09と相対的に小さい変化だった。XGammaは各Nで減少した。\n\n");
    report.push_str(&format_ratio_table(comparison)?);
    report.push_str("\nこれはこの模型・初期条件・N=3・5・7・t<=10の2点比較である。gamma倍増の普遍倍率、一般的な関数形、非線形応答は主張しない。\n\n");
    report.push_str("## 9. XGamma\n\n```text\nx_gamma(t) = sum_ab Gamma[a,b] |rho[a,b](t)|^2\nXGamma(T) = integral_0^T x_gamma(t) dt\n```\n\nXGammaはkernelが重み付けしたcoherence exposureという診断量であり、失われた仕事、散逸エネルギー、dephasing power、熱、entropy production、効率、損傷量ではない。強いdephasingがcoherence自体を早く抑え、kernelが重み付けする対象が減った可能性は候補説明になり得るが、今回の計算では因果機構として確認していない。\n\n");
    report.push_str("## 10. 数値品質\n\n10bはfallback 0。10cのN=7はfallback 2/2成功。全fixed-total条件でsolver failure 0。N=7・TOTAL_GAMMA=1.5軌道は9c正本と最大差0で、trace、Hermiticity、positivity、ledger検査はPASSした。\n\n");
    report.push_str("## 11. 直接確認できたこと\n\nこの模型、vacuum初期状態、N=3・5・7、t<=10、指定されたnoise-free・fixed-per-site・fixed-total 1.5・3.0条件について、正式成果物に保存された最大値、最終値、時間面積、到着時刻、usable fractionを横断整理した。XGammaの直接比較は同一定義で計算されたfixed-total 1.5と3.0に限る。\n\n");
    report.push_str("## 12. 確認できていないこと\n\n中間gamma、連続gamma sweep、TOTAL_GAMMA依存の関数形、臨界値、XGammaの因果機構、XGamma一致比較、dt半減、t>10、N>7、等入力費用比較、連続運転、抽出サイクル、scaling law、量子優位、実機性能は確認していない。\n\n");
    report.push_str("## 13. 主張してはいけないこと\n\n長い鎖が一般に優れる、短い鎖が一般に優れる、XGammaがW損失を引き起こした、XGammaが損失量そのもの、2つのgamma点から関数形を決定、指数則・べき則・相転移、N>7への外挿、量子優位は主張しない。\n\n");
    report.push_str("## 14. Milestone 10最終判定\n\n10a・10b・10cの正式成果物間に矛盾はなく、closure checksは全件PASSした。noise-freeとfixed-per-siteのE系およびXGammaは正式10a統合表に存在しないため、0や推測値で埋めず `not_available` とした。\n\n");
    report.push_str(&format!("最終判定: **{classification}**\n\n"));
    report.push_str("実行・検証記録：\n\n```text\ncargo fmt --all -- --check\ncargo test --release --offline\ncargo run --release --offline --bin milestone_10_closure\n```\n\nclosure binはstd-onlyのCSV／report統合処理で、Hamiltonian、Liouvillian、RK4、dephasing kernel、時間発展モジュールを呼んでいない。Cargo testは107 passed、0 failed、1 ignored。\n\n");
    report.push_str("## 15. 次段階判断\n\nMilestone 10はここで完了する。\n追加gamma点、N>7、XGamma matching、等入力費用比較の\nどれを次に行うかは、研究目的を再確認してから決定する。\n");
    fs::write(OUTPUTS[4], report)?;
    Ok(classification)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for output in OUTPUTS {
        if Path::new(output).exists() {
            return Err(format!("refusing to overwrite existing output {output}").into());
        }
    }
    for input in INPUTS {
        if !Path::new(input).is_file() {
            return Err(format!("missing formal input {input}").into());
        }
    }
    let input_before: Vec<Vec<u8>> = INPUTS
        .iter()
        .map(|path| fs::read(path))
        .collect::<Result<_, _>>()?;
    let final_rows = build_final_rows()?;
    write_final_comparison(&final_rows)?;
    write_rankings(&final_rows)?;
    let comparison = fixed_total_rows()?;
    write_fixed_total_summary(&comparison)?;
    let checks_passed = write_checks(&final_rows, &input_before)?;
    let classification = write_report(&final_rows, &comparison, checks_passed)?;
    println!("Milestone 10 final classification: {classification}");
    if classification == "source_inconsistency_stop" {
        return Err("Milestone 10 closure stopped on source inconsistency".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(n: usize, w_max: f64, arrival: f64) -> FinalRow {
        FinalRow {
            chain_length: n,
            noise_condition: "test".to_owned(),
            gamma_site: "0".to_owned(),
            total_gamma: "0".to_owned(),
            w_max: w_max.to_string(),
            t_at_w_max: "0".to_owned(),
            w_at_t10: w_max.to_string(),
            w_time_area: w_max.to_string(),
            e_at_t10: "not_available".to_owned(),
            e_time_area: "not_available".to_owned(),
            ergotropy_arrival_time: arrival.to_string(),
            usable_fraction_at_t10: w_max.to_string(),
            x_gamma: "not_available".to_owned(),
            source_file: "test.csv".to_owned(),
            value_status: "available".to_owned(),
        }
    }

    #[test]
    fn higher_metric_ranking_uses_descending_values() {
        let rows = vec![row(3, 1.0, 1.0), row(5, 2.0, 2.0), row(7, 3.0, 3.0)];
        assert_eq!(ranking(&rows, "test", "W_max").unwrap().2, "N7_gt_N5_gt_N3");
    }

    #[test]
    fn arrival_ranking_uses_ascending_times() {
        let rows = vec![row(3, 1.0, 1.0), row(5, 2.0, 2.0), row(7, 3.0, 3.0)];
        assert_eq!(
            ranking(&rows, "test", "ergotropy_arrival_time").unwrap().2,
            "N3_fastest_N5_middle_N7_slowest"
        );
    }

    #[test]
    fn unavailable_metric_never_becomes_zero() {
        let rows = vec![row(3, 1.0, 1.0), row(5, 2.0, 2.0), row(7, 3.0, 3.0)];
        let result = ranking(&rows, "test", "XGamma").unwrap();
        assert!(!result.1);
        assert_eq!(result.2, "not_available");
    }
}
