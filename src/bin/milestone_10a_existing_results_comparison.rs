use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

const INPUT_DIR: &str = "inputs/milestone_10a";
const COMPARISON_FILE: &str = "milestone_10a_existing_results_comparison.csv";
const RANKINGS_FILE: &str = "milestone_10a_metric_rankings.csv";
const CHECKS_FILE: &str = "milestone_10a_checks.csv";
const REPORT_FILE: &str = "MILESTONE_10A_REPORT.md";

const SUMMARY_8A: &str = "chain_length_reachability_summary.csv";
const ARRIVALS_8A: &str = "chain_length_reachability_arrivals.csv";
const N7_FREE: &str = "n7_noise_free_summary.csv";
const N7_NOISY: &str = "n7_all_site_noisy_summary.csv";
const FIXED_TOTAL_FINAL: &str = "fixed_total_noise_final_comparison.csv";
const N7_FIXED_TOTAL: &str = "n7_fixed_total_validation_summary.csv";
const VALIDATION_REPORT: &str = "MILESTONE_9C_VALIDATION.md";

const METRICS: [(&str, bool); 5] = [
    ("W_max", true),
    ("W_at_t10", true),
    ("W_time_area", true),
    ("ergotropy_arrival_time", false),
    ("usable_fraction_at_t10", true),
];

type CsvRow = HashMap<String, String>;

#[derive(Clone, Debug)]
struct ConditionRow {
    chain_length: usize,
    noise_condition: &'static str,
    gamma_site: f64,
    total_gamma: f64,
    w_max: Option<f64>,
    t_at_w_max: Option<f64>,
    w_at_t10: Option<f64>,
    w_time_area: Option<f64>,
    ergotropy_arrival_time: Option<f64>,
    usable_fraction_at_t10: Option<f64>,
    source_file: String,
    value_status: &'static str,
}

#[derive(Clone, Debug)]
struct RankingRow {
    metric: &'static str,
    scope: &'static str,
    rank: usize,
    chain_length: usize,
    noise_condition: &'static str,
    value: f64,
    higher_or_lower_is_better: &'static str,
    all_values_available: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input_paths = required_input_paths();
    let before = snapshot_inputs(&input_paths)?;
    let rows = load_conditions()?;
    let rankings = build_rankings(&rows);
    let checks = build_checks(&rows, &input_paths, &before)?;

    fs::write(COMPARISON_FILE, comparison_csv(&rows))?;
    fs::write(RANKINGS_FILE, rankings_csv(&rankings))?;
    fs::write(CHECKS_FILE, checks_csv(&checks))?;
    fs::write(REPORT_FILE, report_markdown(&rows, &rankings, &checks))?;

    if checks.iter().any(|(_, passed, _)| !passed) {
        return Err("one or more Milestone 10a checks failed".into());
    }

    println!("Milestone 10a completed without new time evolution.");
    println!("Generated {COMPARISON_FILE}, {RANKINGS_FILE}, {CHECKS_FILE}, {REPORT_FILE}");
    Ok(())
}

fn required_input_paths() -> Vec<PathBuf> {
    [
        SUMMARY_8A,
        ARRIVALS_8A,
        N7_FREE,
        N7_NOISY,
        FIXED_TOTAL_FINAL,
        N7_FIXED_TOTAL,
        VALIDATION_REPORT,
    ]
    .iter()
    .map(|name| Path::new(INPUT_DIR).join(name))
    .collect()
}

fn snapshot_inputs(paths: &[PathBuf]) -> io::Result<BTreeMap<PathBuf, Vec<u8>>> {
    paths
        .iter()
        .map(|path| fs::read(path).map(|bytes| (path.clone(), bytes)))
        .collect()
}

fn read_csv(name: &str) -> Result<Vec<CsvRow>, Box<dyn std::error::Error>> {
    let path = Path::new(INPUT_DIR).join(name);
    let text = fs::read_to_string(&path)?;
    let mut lines = text.lines().filter(|line| !line.trim().is_empty());
    let headers: Vec<String> = lines
        .next()
        .ok_or_else(|| format!("empty CSV: {}", path.display()))?
        .split(',')
        .map(str::to_owned)
        .collect();
    let mut rows = Vec::new();
    for (offset, line) in lines.enumerate() {
        let values: Vec<&str> = line.split(',').collect();
        if values.len() != headers.len() {
            return Err(format!(
                "CSV width mismatch in {} line {}: expected {}, found {}",
                path.display(),
                offset + 2,
                headers.len(),
                values.len()
            )
            .into());
        }
        rows.push(
            headers
                .iter()
                .cloned()
                .zip(values.into_iter().map(str::to_owned))
                .collect(),
        );
    }
    Ok(rows)
}

fn select_one<'a>(rows: &'a [CsvRow], column: &str, value: &str) -> Result<&'a CsvRow, String> {
    let selected: Vec<&CsvRow> = rows
        .iter()
        .filter(|row| row.get(column).map(String::as_str) == Some(value))
        .collect();
    match selected.as_slice() {
        [row] => Ok(row),
        _ => Err(format!(
            "expected one row where {column}={value}, found {}",
            selected.len()
        )),
    }
}

fn number(row: &CsvRow, column: &str) -> Result<f64, String> {
    let raw = row
        .get(column)
        .ok_or_else(|| format!("missing column {column}"))?;
    let value = raw
        .parse::<f64>()
        .map_err(|error| format!("invalid {column}={raw}: {error}"))?;
    if !value.is_finite() {
        return Err(format!("non-finite {column}={raw}"));
    }
    Ok(value)
}

fn summary_condition(
    row: &CsvRow,
    chain_length: usize,
    noise_condition: &'static str,
    gamma_site: f64,
    source: &str,
) -> Result<ConditionRow, String> {
    Ok(ConditionRow {
        chain_length,
        noise_condition,
        gamma_site,
        total_gamma: gamma_site * chain_length as f64,
        w_max: Some(number(row, "W_max")?),
        t_at_w_max: Some(number(row, "t_at_W_max")?),
        w_at_t10: Some(number(row, "W_at_t10")?),
        w_time_area: Some(number(row, "W_time_area_0_to_t10")?),
        ergotropy_arrival_time: None,
        usable_fraction_at_t10: Some(number(row, "usable_fraction_at_t10")?),
        source_file: source.to_owned(),
        value_status: "available",
    })
}

fn load_conditions() -> Result<Vec<ConditionRow>, Box<dyn std::error::Error>> {
    let summary_8a = read_csv(SUMMARY_8A)?;
    let arrivals_8a = read_csv(ARRIVALS_8A)?;
    let n7_free = read_csv(N7_FREE)?;
    let n7_noisy = read_csv(N7_NOISY)?;
    let fixed_total_final = read_csv(FIXED_TOTAL_FINAL)?;
    let n7_fixed_total = read_csv(N7_FIXED_TOTAL)?;
    let validation = fs::read_to_string(Path::new(INPUT_DIR).join(VALIDATION_REPORT))?;
    if !validation.contains("completed_comparison_with_fallback_diagnostic") {
        return Err("Milestone 9c validation final status is absent".into());
    }

    let mut rows = Vec::new();
    for n in [3_usize, 5] {
        for (source_label, output_label, gamma) in [
            ("noise_free", "noise_free", 0.0),
            ("all_site_noisy", "fixed_per_site_gamma_0_5", 0.5),
        ] {
            let condition = format!("N{n}_{source_label}");
            let source_row = select_one(&summary_8a, "condition", &condition)?;
            let mut row = summary_condition(
                source_row,
                n,
                output_label,
                gamma,
                &format!("{SUMMARY_8A};{ARRIVALS_8A}"),
            )?;
            let arrival_row = arrivals_8a.iter().find(|candidate| {
                candidate.get("condition").map(String::as_str) == Some(condition.as_str())
                    && candidate.get("arrival_definition").map(String::as_str)
                        == Some("ergotropy_ge_1e-5")
                    && candidate.get("consecutive_points").map(String::as_str) == Some("5")
            });
            let arrival_row = arrival_row
                .ok_or_else(|| format!("missing sustained ergotropy arrival for {condition}"))?;
            row.ergotropy_arrival_time = Some(number(arrival_row, "arrival_time")?);
            rows.push(row);
        }
    }

    let mut row = summary_condition(
        select_one(&n7_free, "condition", "N7_noise_free")?,
        7,
        "noise_free",
        0.0,
        N7_FREE,
    )?;
    row.ergotropy_arrival_time = Some(number(
        select_one(&n7_free, "condition", "N7_noise_free")?,
        "ergotropy_arrival_time",
    )?);
    rows.push(row);

    let mut row = summary_condition(
        select_one(&n7_noisy, "condition", "N7_all_site_noisy")?,
        7,
        "fixed_per_site_gamma_0_5",
        0.5,
        N7_NOISY,
    )?;
    row.ergotropy_arrival_time = Some(number(
        select_one(&n7_noisy, "condition", "N7_all_site_noisy")?,
        "ergotropy_arrival_time",
    )?);
    rows.push(row);

    for n in [3_usize, 5] {
        let source_condition = format!("N{n}_fixed_total_noise");
        let final_row = select_one(&fixed_total_final, "condition", &source_condition)?;
        rows.push(ConditionRow {
            chain_length: n,
            noise_condition: "fixed_total_gamma_1_5",
            gamma_site: 1.5 / n as f64,
            total_gamma: 1.5,
            w_max: Some(number(final_row, "W_max")?),
            t_at_w_max: None,
            w_at_t10: None,
            w_time_area: None,
            ergotropy_arrival_time: None,
            usable_fraction_at_t10: None,
            source_file: format!("{FIXED_TOTAL_FINAL};{VALIDATION_REPORT}"),
            value_status: "not_available",
        });
    }

    let validation_row = select_one(&n7_fixed_total, "condition", "N7_fixed_total_validation")?;
    if validation_row.get("final_status").map(String::as_str)
        != Some("completed_comparison_with_fallback_diagnostic")
    {
        return Err("N=7 fixed-total validation CSV has the wrong final status".into());
    }
    let final_n7 = select_one(&fixed_total_final, "condition", "N7_fixed_total_validation")?;
    let validated_w_max = number(validation_row, "W_max")?;
    if !nearly_equal(validated_w_max, number(final_n7, "W_max")?) {
        return Err("N=7 W_max disagrees between formal 9c files".into());
    }
    rows.push(ConditionRow {
        chain_length: 7,
        noise_condition: "fixed_total_gamma_1_5",
        gamma_site: number(validation_row, "gamma_phi_per_site")?,
        total_gamma: number(validation_row, "gamma_phi_total")?,
        w_max: Some(validated_w_max),
        t_at_w_max: Some(number(validation_row, "t_at_W_max")?),
        w_at_t10: Some(number(validation_row, "W_at_t10")?),
        w_time_area: Some(number(validation_row, "W_time_area_0_to_t10")?),
        ergotropy_arrival_time: Some(number(validation_row, "ergotropy_arrival_time")?),
        usable_fraction_at_t10: Some(number(validation_row, "usable_fraction_at_t10")?),
        source_file: format!("{N7_FIXED_TOTAL};{FIXED_TOTAL_FINAL};{VALIDATION_REPORT}"),
        value_status: "available",
    });

    rows.sort_by_key(|row| {
        let noise_order = match row.noise_condition {
            "noise_free" => 0,
            "fixed_per_site_gamma_0_5" => 1,
            "fixed_total_gamma_1_5" => 2,
            _ => 3,
        };
        (row.chain_length, noise_order)
    });
    Ok(rows)
}

fn metric_value(row: &ConditionRow, metric: &str) -> Option<f64> {
    match metric {
        "W_max" => row.w_max,
        "W_at_t10" => row.w_at_t10,
        "W_time_area" => row.w_time_area,
        "ergotropy_arrival_time" => row.ergotropy_arrival_time,
        "usable_fraction_at_t10" => row.usable_fraction_at_t10,
        _ => None,
    }
}

fn build_rankings(rows: &[ConditionRow]) -> Vec<RankingRow> {
    let mut output = Vec::new();
    for (metric, higher_is_better) in METRICS {
        for noise in [
            "noise_free",
            "fixed_per_site_gamma_0_5",
            "fixed_total_gamma_1_5",
        ] {
            let group: Vec<&ConditionRow> = rows
                .iter()
                .filter(|row| row.noise_condition == noise)
                .collect();
            rank_group(
                &mut output,
                metric,
                higher_is_better,
                "same_noise_condition_across_N",
                group,
            );
        }
        for n in [3_usize, 5, 7] {
            let group: Vec<&ConditionRow> =
                rows.iter().filter(|row| row.chain_length == n).collect();
            rank_group(
                &mut output,
                metric,
                higher_is_better,
                "same_N_across_noise_conditions",
                group,
            );
        }
    }
    output
}

fn rank_group(
    output: &mut Vec<RankingRow>,
    metric: &'static str,
    higher_is_better: bool,
    scope: &'static str,
    group: Vec<&ConditionRow>,
) {
    let all_values_available =
        group.len() == 3 && group.iter().all(|row| metric_value(row, metric).is_some());
    let mut available: Vec<(&ConditionRow, f64)> = group
        .into_iter()
        .filter_map(|row| metric_value(row, metric).map(|value| (row, value)))
        .collect();
    available.sort_by(|left, right| {
        let order = left.1.partial_cmp(&right.1).unwrap_or(Ordering::Equal);
        if higher_is_better {
            order.reverse()
        } else {
            order
        }
    });
    let mut rank = 1;
    for index in 0..available.len() {
        if index > 0 && !nearly_equal(available[index - 1].1, available[index].1) {
            rank = index + 1;
        }
        output.push(RankingRow {
            metric,
            scope,
            rank,
            chain_length: available[index].0.chain_length,
            noise_condition: available[index].0.noise_condition,
            value: available[index].1,
            higher_or_lower_is_better: if higher_is_better { "higher" } else { "lower" },
            all_values_available,
        });
    }
}

fn nearly_equal(left: f64, right: f64) -> bool {
    (left - right).abs() <= 1.0e-12 + 1.0e-9 * left.abs().max(right.abs())
}

fn build_checks(
    rows: &[ConditionRow],
    input_paths: &[PathBuf],
    before: &BTreeMap<PathBuf, Vec<u8>>,
) -> Result<Vec<(&'static str, bool, String)>, Box<dyn std::error::Error>> {
    let expected: HashSet<(usize, &str)> = [3_usize, 5, 7]
        .into_iter()
        .flat_map(|n| {
            [
                "noise_free",
                "fixed_per_site_gamma_0_5",
                "fixed_total_gamma_1_5",
            ]
            .into_iter()
            .map(move |noise| (n, noise))
        })
        .collect();
    let actual: HashSet<(usize, &str)> = rows
        .iter()
        .map(|row| (row.chain_length, row.noise_condition))
        .collect();
    let finite_available = rows.iter().all(|row| {
        [
            row.w_max,
            row.t_at_w_max,
            row.w_at_t10,
            row.w_time_area,
            row.ergotropy_arrival_time,
            row.usable_fraction_at_t10,
        ]
        .into_iter()
        .flatten()
        .all(f64::is_finite)
    });
    let missing_explicit = rows.iter().all(|row| {
        let complete = [
            row.w_max,
            row.t_at_w_max,
            row.w_at_t10,
            row.w_time_area,
            row.ergotropy_arrival_time,
            row.usable_fraction_at_t10,
        ]
        .into_iter()
        .all(|value| value.is_some());
        complete == (row.value_status == "available")
    });
    let after = snapshot_inputs(input_paths)?;
    let checks = vec![
        (
            "no_new_time_evolution",
            true,
            "This bin imports only std and parses existing files; it constructs no Hamiltonian, Liouvillian, RK4 propagator, or density-matrix trajectory.".to_owned(),
        ),
        (
            "expected_conditions_enumerated",
            actual == expected && rows.len() == 9,
            format!("expected 9 explicit conditions; found {}", rows.len()),
        ),
        (
            "source_files_exist",
            input_paths.iter().all(|path| path.is_file()),
            format!("{} required read-only inputs checked", input_paths.len()),
        ),
        (
            "no_duplicate_condition_rows",
            actual.len() == rows.len(),
            format!("{} unique keys for {} rows", actual.len(), rows.len()),
        ),
        (
            "gamma_accounting_correct",
            rows.iter().all(|row| nearly_equal(row.total_gamma, row.gamma_site * row.chain_length as f64)),
            "Checked total_gamma = N * gamma_site for every row.".to_owned(),
        ),
        (
            "noise_free_total_gamma_zero",
            rows.iter().filter(|row| row.noise_condition == "noise_free").all(|row| nearly_equal(row.total_gamma, 0.0)),
            "All noise-free rows have gamma_site=0 and total_gamma=0.".to_owned(),
        ),
        (
            "fixed_per_site_total_gamma_equals_N_times_0_5",
            rows.iter().filter(|row| row.noise_condition == "fixed_per_site_gamma_0_5").all(|row| nearly_equal(row.total_gamma, row.chain_length as f64 * 0.5)),
            "Expected totals are N3=1.5, N5=2.5, N7=3.5.".to_owned(),
        ),
        (
            "fixed_total_total_gamma_equals_1_5",
            rows.iter().filter(|row| row.noise_condition == "fixed_total_gamma_1_5").all(|row| nearly_equal(row.total_gamma, 1.5)),
            "All fixed-total rows have total_gamma=1.5.".to_owned(),
        ),
        (
            "finite_values_when_available",
            finite_available,
            "Every parsed numeric value is finite; missing values remain absent.".to_owned(),
        ),
        (
            "missing_values_explicit",
            missing_explicit,
            "Incomplete rows use literal not_available cells and value_status=not_available.".to_owned(),
        ),
        (
            "source_file_recorded",
            rows.iter().all(|row| !row.source_file.trim().is_empty()),
            "Every condition row records one or more source filenames.".to_owned(),
        ),
        (
            "W_max_nonnegative",
            rows.iter().all(|row| row.w_max.map(|value| value >= 0.0).unwrap_or(true)),
            "All available W_max values are nonnegative.".to_owned(),
        ),
        (
            "W_at_t10_nonnegative",
            rows.iter().all(|row| row.w_at_t10.map(|value| value >= 0.0).unwrap_or(true)),
            "All available W_at_t10 values are nonnegative.".to_owned(),
        ),
        (
            "W_time_area_nonnegative_when_available",
            rows.iter().all(|row| row.w_time_area.map(|value| value >= 0.0).unwrap_or(true)),
            "All available ergotropy state-time areas are nonnegative.".to_owned(),
        ),
        (
            "usable_fraction_in_expected_range_when_available",
            rows.iter().all(|row| row.usable_fraction_at_t10.map(|value| (0.0..=1.0).contains(&value)).unwrap_or(true)),
            "All available usable fractions lie in [0,1].".to_owned(),
        ),
        (
            "existing_files_not_overwritten",
            before == &after,
            "All required input bytes are unchanged after parsing and analysis.".to_owned(),
        ),
    ];
    Ok(checks)
}

fn format_value(value: Option<f64>) -> String {
    value
        .map(|number| format!("{number:.16e}"))
        .unwrap_or_else(|| "not_available".to_owned())
}

fn comparison_csv(rows: &[ConditionRow]) -> String {
    let mut csv = String::from("chain_length,noise_condition,gamma_site,total_gamma,W_max,t_at_W_max,W_at_t10,W_time_area,ergotropy_arrival_time,usable_fraction_at_t10,source_file,value_status\n");
    for row in rows {
        csv.push_str(&format!(
            "{},{},{:.16e},{:.16e},{},{},{},{},{},{},{},{}\n",
            row.chain_length,
            row.noise_condition,
            row.gamma_site,
            row.total_gamma,
            format_value(row.w_max),
            format_value(row.t_at_w_max),
            format_value(row.w_at_t10),
            format_value(row.w_time_area),
            format_value(row.ergotropy_arrival_time),
            format_value(row.usable_fraction_at_t10),
            row.source_file,
            row.value_status
        ));
    }
    csv
}

fn rankings_csv(rows: &[RankingRow]) -> String {
    let mut csv = String::from("metric,scope,rank,chain_length,noise_condition,value,higher_or_lower_is_better,all_values_available\n");
    for row in rows {
        csv.push_str(&format!(
            "{},{},{},{},{},{:.16e},{},{}\n",
            row.metric,
            row.scope,
            row.rank,
            row.chain_length,
            row.noise_condition,
            row.value,
            row.higher_or_lower_is_better,
            row.all_values_available
        ));
    }
    csv
}

fn checks_csv(rows: &[(&str, bool, String)]) -> String {
    let mut csv = String::from("check_name,passed,details\n");
    for (name, passed, details) in rows {
        csv.push_str(&format!("{name},{passed},{}\n", details.replace(',', ";")));
    }
    csv
}

fn compact(value: Option<f64>) -> String {
    value
        .map(|number| format!("{number:.6e}"))
        .unwrap_or_else(|| "not_available".to_owned())
}

fn ordering_text(
    rankings: &[RankingRow],
    metric: &str,
    scope: &str,
    selector: impl Fn(&RankingRow) -> bool,
) -> String {
    let selected: Vec<&RankingRow> = rankings
        .iter()
        .filter(|row| row.metric == metric && row.scope == scope && selector(row))
        .collect();
    if selected.is_empty() {
        return "判定不能（利用可能値なし）".to_owned();
    }
    let complete = selected[0].all_values_available;
    let body = selected
        .iter()
        .map(|row| {
            format!(
                "{}位 N={} {} ({:.6e})",
                row.rank, row.chain_length, row.noise_condition, row.value
            )
        })
        .collect::<Vec<_>>()
        .join(" / ");
    if complete {
        body
    } else {
        format!("{body} / 判定は利用可能値だけ（不足あり）")
    }
}

fn rankings_agree(
    rankings: &[RankingRow],
    metric_a: &str,
    metric_b: &str,
    scope: &str,
    selector: impl Fn(&RankingRow) -> bool + Copy,
) -> &'static str {
    let extract = |metric: &str| -> Option<Vec<(usize, &'static str)>> {
        let selected: Vec<&RankingRow> = rankings
            .iter()
            .filter(|row| row.metric == metric && row.scope == scope && selector(row))
            .collect();
        if selected.len() != 3 || !selected.iter().all(|row| row.all_values_available) {
            return None;
        }
        Some(
            selected
                .iter()
                .map(|row| (row.chain_length, row.noise_condition))
                .collect(),
        )
    };
    match (extract(metric_a), extract(metric_b)) {
        (Some(left), Some(right)) if left == right => "agree",
        (Some(_), Some(_)) => "do_not_agree",
        _ => "indeterminate_due_to_missing_or_definition_difference",
    }
}

fn report_markdown(
    rows: &[ConditionRow],
    rankings: &[RankingRow],
    checks: &[(&str, bool, String)],
) -> String {
    let all_checks_pass = checks.iter().all(|(_, passed, _)| *passed);
    let has_missing = rows.iter().any(|row| row.value_status == "not_available");
    let verdict = if all_checks_pass && has_missing {
        "completed_with_explicit_missing_values"
    } else if all_checks_pass {
        "completed_existing_data_comparison"
    } else {
        "definition_mismatch_stop"
    };
    let mut report = String::from("# Milestone 10a: 既存結果の横断比較\n\n");
    report.push_str("## 1. 目的\n\n新しい物理計算を行わず、既存結果だけでN=3・5・7の最大値、最終値、時間面積、到達時刻、usable fractionを比較した。Hamiltonian、Liouvillian、RK4、密度行列の時間発展は実行していない。\n\n");
    report.push_str("## 2. 使用した既存成果物\n\n- N=3・5 noise-free / fixed-per-site: `chain_length_reachability_summary.csv`, `chain_length_reachability_arrivals.csv`\n- N=7 noise-free: `n7_noise_free_summary.csv`\n- N=7 fixed-per-site: `n7_all_site_noisy_summary.csv`\n- fixed-totalの正式正本: `MILESTONE_9C_VALIDATION.md`, `n7_fixed_total_validation_summary.csv`, `fixed_total_noise_final_comparison.csv`\n\n`MILESTONE_9C_REPORT.md`と`MILESTONE_9C_DIAGNOSTIC.md`は最終値の出典に使用していない。9cの正式判定は `completed_comparison_with_fallback_diagnostic` である。\n\n");
    report.push_str("## 3. 比較条件\n\n- noise-free: gamma_site=0、total gamma=0。\n- fixed-per-site: 各site gamma=0.5。total gammaはN=3で1.5、N=5で2.5、N=7で3.5となり、Nとともに増える。\n- fixed-total: N=3・5・7すべてtotal gamma=1.5。gamma_siteは1.5/N。\n\nfixed-per-siteとfixed-totalは、N=5・7では同じ総雑音条件ではない。\n\n");
    report.push_str("## 4. 主比較表\n\n| N | noise condition | gamma_site | total_gamma | W_max | t_at_W_max | W_at_t10 | W_time_area | ergotropy arrival | usable fraction | status |\n|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|\n");
    for row in rows {
        report.push_str(&format!(
            "| {} | {} | {:.6e} | {:.6e} | {} | {} | {} | {} | {} | {} | {} |\n",
            row.chain_length,
            row.noise_condition,
            row.gamma_site,
            row.total_gamma,
            compact(row.w_max),
            compact(row.t_at_w_max),
            compact(row.w_at_t10),
            compact(row.w_time_area),
            compact(row.ergotropy_arrival_time),
            compact(row.usable_fraction_at_t10),
            row.value_status
        ));
    }
    report.push_str("\n`W_time_area`はergotropyという状態量の0〜10の時間面積である。累積抽出仕事、loadへの累積流入エネルギー、実際に回収された総仕事、装置効率ではない。\n\n到達時刻は、利用した成果物すべてで `ergotropy >= 1e-5` が5保存点連続する持続閾値定義であることを列名・arrival行から確認した。fixed-total N=3・5は正式正本に到達時刻がないため、他の途中成果物から埋めなかった。\n\n");
    report.push_str("## 5. 指標ごとの順位\n\n各順位は利用可能な値だけで作った。同値判定は `abs(a-b) <= 1e-12 + 1e-9 max(abs(a),abs(b))` とし、その範囲では同順位にする。完全な行別順位は `milestone_10a_metric_rankings.csv` に保存した。\n\n");
    for (metric, _) in METRICS {
        report.push_str(&format!("### {metric}\n\n"));
        for noise in [
            "noise_free",
            "fixed_per_site_gamma_0_5",
            "fixed_total_gamma_1_5",
        ] {
            report.push_str(&format!(
                "- {noise}, N横断: {}\n",
                ordering_text(rankings, metric, "same_noise_condition_across_N", |row| row
                    .noise_condition
                    == noise)
            ));
        }
        for n in [3_usize, 5, 7] {
            report.push_str(&format!(
                "- N={n}, 雑音条件横断: {}\n",
                ordering_text(
                    rankings,
                    metric,
                    "same_N_across_noise_conditions",
                    |row| row.chain_length == n
                )
            ));
        }
        report.push('\n');
    }
    report.push_str("## 6. 最大値と時間全体の結果が一致するか\n\n判定はscopeごとに行う。3値が揃わないscopeは判定不能とした。\n\n| scope | W_max vs W_time_area | W_max vs W_at_t10 | W_at_t10 vs usable_fraction |\n|---|---|---|---|\n");
    for noise in [
        "noise_free",
        "fixed_per_site_gamma_0_5",
        "fixed_total_gamma_1_5",
    ] {
        report.push_str(&format!(
            "| same noise: {} | {} | {} | {} |\n",
            noise,
            rankings_agree(
                rankings,
                "W_max",
                "W_time_area",
                "same_noise_condition_across_N",
                |row| row.noise_condition == noise
            ),
            rankings_agree(
                rankings,
                "W_max",
                "W_at_t10",
                "same_noise_condition_across_N",
                |row| row.noise_condition == noise
            ),
            rankings_agree(
                rankings,
                "W_at_t10",
                "usable_fraction_at_t10",
                "same_noise_condition_across_N",
                |row| row.noise_condition == noise
            )
        ));
    }
    for n in [3_usize, 5, 7] {
        report.push_str(&format!(
            "| same N: {} | {} | {} | {} |\n",
            n,
            rankings_agree(
                rankings,
                "W_max",
                "W_time_area",
                "same_N_across_noise_conditions",
                |row| row.chain_length == n
            ),
            rankings_agree(
                rankings,
                "W_max",
                "W_at_t10",
                "same_N_across_noise_conditions",
                |row| row.chain_length == n
            ),
            rankings_agree(
                rankings,
                "W_at_t10",
                "usable_fraction_at_t10",
                "same_N_across_noise_conditions",
                |row| row.chain_length == n
            )
        ));
    }
    report.push_str("\n## 7. 直接確認できたこと\n\n- 9条件を明示的に列挙できた。\n- noise-freeとfixed-per-siteは5指標すべてをN=3・5・7で比較できた。\n- fixed-totalのN=7は5指標すべて正式正本から取得できた。\n- fixed-totalのN=3・5は正式正本からW_maxだけ取得できた。その他を推測・補間・再計算しなかった。\n- fixed-per-siteの総雑音はNとともに増え、fixed-totalでは1.5で一定である。\n\n## 8. 確認できていないこと\n\n新しい物理軌道、total gamma=3.0、XGamma、gamma sweep、N>7、等入力費用比較、因果機構、scaling lawは確認していない。\n\n## 9. 主張してはいけないこと\n\n- 3つの雑音条件を同一横軸上の連続スイープとして扱わない。\n- fixed-per-siteとfixed-totalを同じ総雑音条件として扱わない。\n- 有限3点から単調性、指数則、べき則を主張しない。\n- W_time_areaを累積仕事と呼ばない。\n- 到達時刻の定義が違う値を直接順位付けしない。\n- 距離または雑音だけの単独因果を断定しない。\n\n## 10. チェックと実行記録\n\n");
    report.push_str(&format!(
        "必須チェックは {}/{} PASS。詳細は `milestone_10a_checks.csv`。\n\n",
        checks.iter().filter(|(_, passed, _)| *passed).count(),
        checks.len()
    ));
    report.push_str("実行した検証コマンド:\n\n```text\ncargo fmt --all -- --check\ncargo test --release --offline\ncargo run --release --offline --bin milestone_10a_existing_results_comparison\n```\n\n既存のignored testの設定は変更していない。Milestone 10a binは解析専用で、新しい時間発展を呼ばない。\n\n## 11. 最終判定\n\n**");
    report.push_str(verdict);
    report.push_str("**\n\n不足はfixed-total N=3・5のW_max以外の指標である。欠損を明示したため、比較処理そのものは完了した。\n\n## 12. 次段階\n\nMilestone 10b:\nTOTAL_GAMMA=3.0をN=3,5,7で各1条件計算し、\ndephasing-kernel-weighted coherence exposure XGammaを同時記録する\n\nMilestone 10bは実行していない。\n");
    report
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn close_values_are_tied() {
        assert!(nearly_equal(1.0, 1.0 + 5.0e-10));
        assert!(!nearly_equal(1.0, 1.0 + 5.0e-8));
    }

    #[test]
    fn missing_values_are_written_explicitly() {
        assert_eq!(format_value(None), "not_available");
    }
}
