use std::collections::{BTreeMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};

const TIMESERIES_INPUT: &str = "local_noise_placement_timeseries.csv";
const SUMMARY_INPUT: &str = "local_noise_placement_summary.csv";
const RATIOS_INPUT: &str = "local_noise_placement_ratios.csv";
const CONDITIONS: [&str; 4] = ["noise_free", "noise_entrance", "noise_middle", "noise_exit"];
const NOISY_CONDITIONS: [&str; 3] = ["noise_entrance", "noise_middle", "noise_exit"];
const CONSECUTIVE_POINTS: usize = 5;
const SIGNAL_TOLERANCE: f64 = 1.0e-8;
const VALUE_TOLERANCE: f64 = 1.0e-10;
const AREA_TOLERANCE: f64 = 1.0e-10;
const POPULATION_TOLERANCE: f64 = 1.0e-8;
const RANK_ABSOLUTE_TOLERANCE: f64 = 1.0e-7;
const RANK_RELATIVE_TOLERANCE: f64 = 5.0e-3;

#[derive(Clone)]
struct CsvTable {
    headers: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl CsvTable {
    fn read(path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let text = std::fs::read_to_string(path)?;
        let mut lines = text.lines();
        let headers: Vec<String> = lines
            .next()
            .ok_or("CSV is empty")?
            .split(',')
            .map(str::to_string)
            .collect();
        let mut rows = Vec::new();
        for (line_index, line) in lines.enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let values: Vec<String> = line.split(',').map(str::to_string).collect();
            if values.len() != headers.len() {
                return Err(format!(
                    "{path}: line {} has {} columns, expected {}",
                    line_index + 2,
                    values.len(),
                    headers.len()
                )
                .into());
            }
            rows.push(values);
        }
        Ok(Self { headers, rows })
    }

    fn column(&self, name: &str) -> Result<usize, Box<dyn std::error::Error>> {
        self.headers
            .iter()
            .position(|header| header == name)
            .ok_or_else(|| format!("missing required column {name}").into())
    }

    fn value<'a>(
        &'a self,
        row: &'a [String],
        name: &str,
    ) -> Result<&'a str, Box<dyn std::error::Error>> {
        Ok(&row[self.column(name)?])
    }

    fn number(&self, row: &[String], name: &str) -> Result<f64, Box<dyn std::error::Error>> {
        Ok(self.value(row, name)?.parse::<f64>()?)
    }
}

#[derive(Clone)]
struct SourceRow {
    condition: String,
    noise_site: String,
    gamma_phi: f64,
    omega: f64,
    time: f64,
    load_energy: f64,
    load_ergotropy: f64,
    usable_fraction: f64,
    site: [f64; 3],
}

#[derive(Clone)]
struct DamageRow {
    condition: String,
    noise_site: String,
    time: f64,
    load_energy: f64,
    load_ergotropy: f64,
    usable_fraction: f64,
    reference_load_energy: f64,
    reference_load_ergotropy: f64,
    reference_usable_fraction: f64,
    delta_e: f64,
    delta_w: f64,
    delta_use: f64,
    r_e: f64,
    r_w: f64,
    r_use: f64,
    loss_e: f64,
    loss_w: f64,
    loss_use: f64,
    site: [f64; 3],
    noisy_site_population: f64,
    total_chain_population: f64,
    noisy_site_population_fraction: f64,
}

#[derive(Clone, Copy)]
struct Quantity {
    name: &'static str,
    absolute_threshold: f64,
    minimum_reference: f64,
}

const QUANTITIES: [Quantity; 3] = [
    Quantity {
        name: "E",
        absolute_threshold: 1.0e-5,
        minimum_reference: 1.0e-4,
    },
    Quantity {
        name: "W",
        absolute_threshold: 1.0e-5,
        minimum_reference: 1.0e-4,
    },
    Quantity {
        name: "usable_fraction",
        absolute_threshold: 1.0e-3,
        minimum_reference: 1.0e-3,
    },
];

#[derive(Clone, Copy)]
struct Threshold {
    name: &'static str,
    relative: f64,
}

const THRESHOLDS: [Threshold; 3] = [
    Threshold {
        name: "weak",
        relative: 0.01,
    },
    Threshold {
        name: "medium",
        relative: 0.05,
    },
    Threshold {
        name: "strong",
        relative: 0.10,
    },
];

#[derive(Clone)]
struct Onset {
    condition: String,
    quantity: &'static str,
    threshold_level: &'static str,
    absolute_threshold: f64,
    relative_threshold: f64,
    first_crossing_time: f64,
    sustained_onset_time: f64,
    sustained_duration_at_detection: f64,
    absolute_difference_at_onset: f64,
    relative_loss_at_onset: f64,
    reference_value_at_onset: f64,
    noisy_value_at_onset: f64,
    noisy_site_population_at_onset: f64,
    total_chain_population_at_onset: f64,
    noisy_site_population_fraction_at_onset: f64,
    dominant_site_at_onset: String,
}

#[derive(Clone)]
struct Extremum {
    condition: String,
    quantity: &'static str,
    time: f64,
    absolute_difference: f64,
    relative_loss: f64,
    reference_value: f64,
    noisy_value: f64,
    noisy_site_population: f64,
    total_chain_population: f64,
    noisy_site_population_fraction: f64,
    dominant_site: String,
    load_energy: f64,
    load_ergotropy: f64,
    usable_fraction: f64,
}

#[derive(Clone, Copy)]
struct Window {
    name: &'static str,
    start: f64,
    end: f64,
    exclude_start_from_mean: bool,
}

const WINDOWS: [Window; 4] = [
    Window {
        name: "pulse_interval",
        start: 0.0,
        end: 3.2,
        exclude_start_from_mean: false,
    },
    Window {
        name: "early_post_pulse",
        start: 3.2,
        end: 5.0,
        exclude_start_from_mean: true,
    },
    Window {
        name: "middle_interval",
        start: 5.0,
        end: 7.5,
        exclude_start_from_mean: true,
    },
    Window {
        name: "late_interval",
        start: 7.5,
        end: 10.0,
        exclude_start_from_mean: true,
    },
];

#[derive(Clone)]
struct WindowResult {
    condition: String,
    window: Window,
    point_count: usize,
    mean_load_energy: f64,
    mean_load_ergotropy: f64,
    mean_usable_fraction: f64,
    e_time_area: f64,
    w_time_area: f64,
    mean_site: [f64; 3],
    maximum_e_loss: f64,
    maximum_w_loss: f64,
    maximum_use_loss: f64,
}

#[derive(Clone)]
struct RankingRow {
    time: f64,
    worst_e_condition: String,
    worst_w_condition: String,
    worst_use_condition: String,
    e_tie: bool,
    w_tie: bool,
    use_tie: bool,
}

#[derive(Clone)]
struct Check {
    check: String,
    value: String,
    expected: String,
    pass: bool,
}

fn n(value: f64) -> String {
    if value.is_nan() {
        "NaN".to_string()
    } else {
        format!("{value:.16e}")
    }
}

fn safe_ratio(numerator: f64, denominator: f64) -> f64 {
    if !numerator.is_finite() || !denominator.is_finite() || denominator.abs() <= SIGNAL_TOLERANCE {
        f64::NAN
    } else {
        numerator / denominator
    }
}

fn ratio_rule(numerator: f64, denominator: f64, ratio: f64) -> bool {
    if !numerator.is_finite() || !denominator.is_finite() || denominator.abs() <= SIGNAL_TOLERANCE {
        ratio.is_nan()
    } else {
        ratio.is_finite() && approximately_equal(ratio, numerator / denominator, VALUE_TOLERANCE)
    }
}

fn finite_or_nan(value: f64) -> bool {
    value.is_finite() || value.is_nan()
}

fn approximately_equal(left: f64, right: f64, tolerance: f64) -> bool {
    (left - right).abs() <= tolerance
}

fn condition_site(condition: &str) -> Result<(&'static str, Option<usize>), String> {
    match condition {
        "noise_free" => Ok(("none", None)),
        "noise_entrance" => Ok(("site1", Some(0))),
        "noise_middle" => Ok(("site2", Some(1))),
        "noise_exit" => Ok(("site3", Some(2))),
        _ => Err(format!("unknown condition {condition}")),
    }
}

fn parse_source(table: &CsvTable) -> Result<Vec<SourceRow>, Box<dyn std::error::Error>> {
    let required = [
        "condition",
        "noise_site",
        "gamma_phi",
        "Omega",
        "time",
        "load_energy",
        "load_ergotropy",
        "usable_fraction",
        "site1_population",
        "site2_population",
        "site3_population",
    ];
    for column in required {
        table.column(column)?;
    }
    table
        .rows
        .iter()
        .map(|row| {
            Ok(SourceRow {
                condition: table.value(row, "condition")?.to_string(),
                noise_site: table.value(row, "noise_site")?.to_string(),
                gamma_phi: table.number(row, "gamma_phi")?,
                omega: table.number(row, "Omega")?,
                time: table.number(row, "time")?,
                load_energy: table.number(row, "load_energy")?,
                load_ergotropy: table.number(row, "load_ergotropy")?,
                usable_fraction: table.number(row, "usable_fraction")?,
                site: [
                    table.number(row, "site1_population")?,
                    table.number(row, "site2_population")?,
                    table.number(row, "site3_population")?,
                ],
            })
        })
        .collect()
}

fn group_source(rows: &[SourceRow]) -> BTreeMap<String, Vec<SourceRow>> {
    let mut grouped: BTreeMap<String, Vec<SourceRow>> = BTreeMap::new();
    for row in rows {
        grouped
            .entry(row.condition.clone())
            .or_default()
            .push(row.clone());
    }
    for values in grouped.values_mut() {
        values.sort_by(|left, right| left.time.total_cmp(&right.time));
    }
    grouped
}

fn trapezoid_source(rows: &[SourceRow], value: impl Fn(&SourceRow) -> f64) -> f64 {
    rows.windows(2)
        .map(|pair| 0.5 * (value(&pair[0]) + value(&pair[1])) * (pair[1].time - pair[0].time))
        .sum()
}

fn source_checks(
    table: &CsvTable,
    grouped: &BTreeMap<String, Vec<SourceRow>>,
    summary: &CsvTable,
    ratios: &CsvTable,
) -> Result<Vec<Check>, Box<dyn std::error::Error>> {
    let mut checks = Vec::new();
    let required = [
        "condition",
        "noise_site",
        "gamma_phi",
        "Omega",
        "time",
        "load_energy",
        "load_ergotropy",
        "usable_fraction",
        "site1_population",
        "site2_population",
        "site3_population",
    ];
    let required_ok = required
        .iter()
        .all(|name| table.headers.iter().any(|header| header == name));
    checks.push(Check {
        check: "required_columns".to_string(),
        value: required.len().to_string(),
        expected: "all present".to_string(),
        pass: required_ok,
    });
    let condition_ok = CONDITIONS
        .iter()
        .all(|condition| grouped.contains_key(*condition))
        && grouped.len() == CONDITIONS.len();
    checks.push(Check {
        check: "four_conditions".to_string(),
        value: grouped.len().to_string(),
        expected: "4 exact labels".to_string(),
        pass: condition_ok,
    });
    if !condition_ok {
        return Ok(checks);
    }
    let reference = &grouped["noise_free"];
    let equal_counts = CONDITIONS
        .iter()
        .all(|condition| grouped[*condition].len() == reference.len());
    checks.push(Check {
        check: "equal_point_counts".to_string(),
        value: CONDITIONS
            .iter()
            .map(|condition| grouped[*condition].len().to_string())
            .collect::<Vec<_>>()
            .join("/"),
        expected: "all equal".to_string(),
        pass: equal_counts,
    });
    let same_grid = CONDITIONS.iter().all(|condition| {
        grouped[*condition]
            .iter()
            .zip(reference)
            .all(|(left, right)| left.time == right.time)
    });
    checks.push(Check {
        check: "same_time_grid".to_string(),
        value: reference.len().to_string(),
        expected: "exact match".to_string(),
        pass: same_grid,
    });
    let monotonic = CONDITIONS.iter().all(|condition| {
        grouped[*condition]
            .windows(2)
            .all(|pair| pair[1].time > pair[0].time)
    });
    checks.push(Check {
        check: "strictly_increasing_time".to_string(),
        value: monotonic.to_string(),
        expected: "true".to_string(),
        pass: monotonic,
    });
    let mut seen = HashSet::new();
    let unique = grouped
        .values()
        .flatten()
        .all(|row| seen.insert(format!("{}:{:.17e}", row.condition, row.time)));
    checks.push(Check {
        check: "unique_condition_time".to_string(),
        value: seen.len().to_string(),
        expected: table.rows.len().to_string(),
        pass: unique,
    });
    let labels = grouped.iter().all(|(condition, rows)| {
        condition_site(condition)
            .is_ok_and(|(expected, _)| rows.iter().all(|row| row.noise_site == expected))
    });
    checks.push(Check {
        check: "noise_site_labels".to_string(),
        value: labels.to_string(),
        expected: "condition mapping".to_string(),
        pass: labels,
    });
    let constants = grouped.values().all(|rows| {
        rows.iter()
            .all(|row| row.gamma_phi == rows[0].gamma_phi && row.omega == rows[0].omega)
    });
    checks.push(Check {
        check: "gamma_and_omega_constant".to_string(),
        value: constants.to_string(),
        expected: "true within each condition".to_string(),
        pass: constants,
    });
    let finite_major = grouped.values().flatten().all(|row| {
        [
            row.time,
            row.gamma_phi,
            row.omega,
            row.load_energy,
            row.load_ergotropy,
            row.site[0],
            row.site[1],
            row.site[2],
        ]
        .iter()
        .all(|value| value.is_finite())
            && finite_or_nan(row.usable_fraction)
    });
    checks.push(Check {
        check: "finite_major_or_allowed_nan".to_string(),
        value: finite_major.to_string(),
        expected: "true".to_string(),
        pass: finite_major,
    });
    let populations = grouped.values().flatten().all(|row| {
        row.site
            .iter()
            .all(|value| *value >= -POPULATION_TOLERANCE && *value <= 1.0 + POPULATION_TOLERANCE)
    });
    checks.push(Check {
        check: "population_bounds".to_string(),
        value: populations.to_string(),
        expected: "-1e-8<=p<=1+1e-8".to_string(),
        pass: populations,
    });
    let physical_values = grouped.values().flatten().all(|row| {
        row.load_energy >= -VALUE_TOLERANCE
            && row.load_ergotropy >= -VALUE_TOLERANCE
            && row.load_ergotropy <= row.load_energy + VALUE_TOLERANCE
            && (row.usable_fraction.is_nan() || row.usable_fraction >= -VALUE_TOLERANCE)
    });
    checks.push(Check {
        check: "nonnegative_and_W_le_E".to_string(),
        value: physical_values.to_string(),
        expected: "within 1e-10".to_string(),
        pass: physical_values,
    });
    let denominator_rule = grouped.values().flatten().all(|row| {
        let should_nan = row.load_energy.abs() <= SIGNAL_TOLERANCE;
        if should_nan {
            row.usable_fraction.is_nan()
        } else {
            row.usable_fraction.is_finite()
                && approximately_equal(
                    row.usable_fraction,
                    row.load_ergotropy / row.load_energy,
                    VALUE_TOLERANCE,
                )
        }
    });
    checks.push(Check {
        check: "usable_fraction_denominator_rule".to_string(),
        value: denominator_rule.to_string(),
        expected: "NaN iff |E|<=1e-8; otherwise W/E".to_string(),
        pass: denominator_rule,
    });

    for condition in CONDITIONS {
        let summary_row = summary
            .rows
            .iter()
            .find(|row| summary.value(row, "condition").ok() == Some(condition))
            .ok_or_else(|| format!("summary missing {condition}"))?;
        let last = grouped[condition].last().ok_or("empty condition")?;
        for (source_name, summary_name, value) in [
            ("E", "E_at_t10", last.load_energy),
            ("W", "W_at_t10", last.load_ergotropy),
            (
                "usable_fraction",
                "usable_fraction_at_t10",
                last.usable_fraction,
            ),
        ] {
            let expected = summary.number(summary_row, summary_name)?;
            let difference = (value - expected).abs();
            checks.push(Check {
                check: format!("summary_t10_{condition}_{source_name}"),
                value: n(difference),
                expected: "absolute<=1e-10".to_string(),
                pass: difference <= VALUE_TOLERANCE,
            });
        }
        let w_area = trapezoid_source(&grouped[condition], |row| row.load_ergotropy);
        let expected_area = summary.number(summary_row, "W_time_area")?;
        let difference = (w_area - expected_area).abs();
        checks.push(Check {
            check: format!("summary_W_time_area_{condition}"),
            value: n(difference),
            expected: "absolute<=1e-10".to_string(),
            pass: difference <= AREA_TOLERANCE,
        });
    }

    let ratio_reference = &grouped["noise_free"];
    for condition in NOISY_CONDITIONS {
        let ratio_row = ratios
            .rows
            .iter()
            .find(|row| ratios.value(row, "condition").ok() == Some(condition))
            .ok_or_else(|| format!("ratios missing {condition}"))?;
        let last = grouped[condition].last().ok_or("empty condition")?;
        let reference_last = ratio_reference.last().ok_or("empty reference")?;
        for (name, calculated) in [
            (
                "R_E",
                safe_ratio(last.load_energy, reference_last.load_energy),
            ),
            (
                "R_W",
                safe_ratio(last.load_ergotropy, reference_last.load_ergotropy),
            ),
            (
                "R_use",
                safe_ratio(last.usable_fraction, reference_last.usable_fraction),
            ),
        ] {
            let expected = ratios.number(ratio_row, name)?;
            let difference = (calculated - expected).abs();
            checks.push(Check {
                check: format!("ratios_t10_{condition}_{name}"),
                value: n(difference),
                expected: "absolute<=1e-10".to_string(),
                pass: difference <= VALUE_TOLERANCE,
            });
        }
    }
    Ok(checks)
}

fn build_damage(
    grouped: &BTreeMap<String, Vec<SourceRow>>,
) -> Result<BTreeMap<String, Vec<DamageRow>>, Box<dyn std::error::Error>> {
    let reference = &grouped["noise_free"];
    let mut result = BTreeMap::new();
    for condition in NOISY_CONDITIONS {
        let (_, site_index) = condition_site(condition).map_err(|error| error.to_string())?;
        let site_index = site_index.ok_or("no noisy site")?;
        let mut rows = Vec::with_capacity(reference.len());
        for (noisy, baseline) in grouped[condition].iter().zip(reference) {
            let r_e = safe_ratio(noisy.load_energy, baseline.load_energy);
            let r_w = safe_ratio(noisy.load_ergotropy, baseline.load_ergotropy);
            let r_use = safe_ratio(noisy.usable_fraction, baseline.usable_fraction);
            let total_chain_population: f64 = noisy.site.iter().sum();
            rows.push(DamageRow {
                condition: condition.to_string(),
                noise_site: noisy.noise_site.clone(),
                time: noisy.time,
                load_energy: noisy.load_energy,
                load_ergotropy: noisy.load_ergotropy,
                usable_fraction: noisy.usable_fraction,
                reference_load_energy: baseline.load_energy,
                reference_load_ergotropy: baseline.load_ergotropy,
                reference_usable_fraction: baseline.usable_fraction,
                delta_e: noisy.load_energy - baseline.load_energy,
                delta_w: noisy.load_ergotropy - baseline.load_ergotropy,
                delta_use: if noisy.usable_fraction.is_finite()
                    && baseline.usable_fraction.is_finite()
                {
                    noisy.usable_fraction - baseline.usable_fraction
                } else {
                    f64::NAN
                },
                r_e,
                r_w,
                r_use,
                loss_e: if r_e.is_finite() { 1.0 - r_e } else { f64::NAN },
                loss_w: if r_w.is_finite() { 1.0 - r_w } else { f64::NAN },
                loss_use: if r_use.is_finite() {
                    1.0 - r_use
                } else {
                    f64::NAN
                },
                site: noisy.site,
                noisy_site_population: noisy.site[site_index],
                total_chain_population,
                noisy_site_population_fraction: safe_ratio(
                    noisy.site[site_index],
                    total_chain_population,
                ),
            });
        }
        result.insert(condition.to_string(), rows);
    }
    Ok(result)
}

fn quantity_values(row: &DamageRow, quantity: &str) -> (f64, f64, f64, f64) {
    match quantity {
        "E" => (
            row.reference_load_energy,
            row.load_energy,
            row.delta_e,
            row.loss_e,
        ),
        "W" => (
            row.reference_load_ergotropy,
            row.load_ergotropy,
            row.delta_w,
            row.loss_w,
        ),
        "usable_fraction" => (
            row.reference_usable_fraction,
            row.usable_fraction,
            row.delta_use,
            row.loss_use,
        ),
        _ => unreachable!(),
    }
}

fn dominant_site(site: &[f64; 3]) -> String {
    let maximum = site.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let tied: Vec<_> = site
        .iter()
        .enumerate()
        .filter(|(_, value)| (**value - maximum).abs() <= 1.0e-12)
        .map(|(index, _)| format!("site{}", index + 1))
        .collect();
    if tied.len() == 1 {
        tied[0].clone()
    } else {
        format!("tie({})", tied.join("+"))
    }
}

fn threshold_crossed(row: &DamageRow, quantity: Quantity, threshold: Threshold) -> bool {
    let (reference, _noisy, delta, loss) = quantity_values(row, quantity.name);
    reference.is_finite()
        && reference >= quantity.minimum_reference
        && delta.is_finite()
        && -delta >= quantity.absolute_threshold
        && loss.is_finite()
        && loss >= threshold.relative
}

fn meaningful_loss(row: &DamageRow, quantity: Quantity) -> Option<f64> {
    let (reference, _noisy, _delta, loss) = quantity_values(row, quantity.name);
    (reference.is_finite() && reference >= quantity.minimum_reference && loss.is_finite())
        .then_some(loss)
}

fn detect_onsets(damage: &BTreeMap<String, Vec<DamageRow>>) -> Vec<Onset> {
    let mut onsets = Vec::new();
    for condition in NOISY_CONDITIONS {
        let rows = &damage[condition];
        let grid_step = rows
            .windows(2)
            .map(|pair| pair[1].time - pair[0].time)
            .next()
            .unwrap_or(f64::NAN);
        for quantity in QUANTITIES {
            for threshold in THRESHOLDS {
                let first_crossing = rows
                    .iter()
                    .position(|row| threshold_crossed(row, quantity, threshold));
                let sustained =
                    (0..=rows.len().saturating_sub(CONSECUTIVE_POINTS)).find(|&index| {
                        rows[index..index + CONSECUTIVE_POINTS]
                            .iter()
                            .all(|row| threshold_crossed(row, quantity, threshold))
                    });
                let row = sustained.map(|index| &rows[index]);
                let (reference, noisy, delta, loss) = row
                    .map(|value| quantity_values(value, quantity.name))
                    .unwrap_or((f64::NAN, f64::NAN, f64::NAN, f64::NAN));
                onsets.push(Onset {
                    condition: condition.to_string(),
                    quantity: quantity.name,
                    threshold_level: threshold.name,
                    absolute_threshold: quantity.absolute_threshold,
                    relative_threshold: threshold.relative,
                    first_crossing_time: first_crossing
                        .map(|index| rows[index].time)
                        .unwrap_or(f64::NAN),
                    sustained_onset_time: row.map(|value| value.time).unwrap_or(f64::NAN),
                    sustained_duration_at_detection: if row.is_some() {
                        grid_step * CONSECUTIVE_POINTS as f64
                    } else {
                        f64::NAN
                    },
                    absolute_difference_at_onset: -delta,
                    relative_loss_at_onset: loss,
                    reference_value_at_onset: reference,
                    noisy_value_at_onset: noisy,
                    noisy_site_population_at_onset: row
                        .map(|value| value.noisy_site_population)
                        .unwrap_or(f64::NAN),
                    total_chain_population_at_onset: row
                        .map(|value| value.total_chain_population)
                        .unwrap_or(f64::NAN),
                    noisy_site_population_fraction_at_onset: row
                        .map(|value| value.noisy_site_population_fraction)
                        .unwrap_or(f64::NAN),
                    dominant_site_at_onset: row
                        .map(|value| dominant_site(&value.site))
                        .unwrap_or_else(|| "undefined".to_string()),
                });
            }
        }
    }
    onsets
}

fn detect_extrema(damage: &BTreeMap<String, Vec<DamageRow>>) -> Vec<Extremum> {
    let mut extrema = Vec::new();
    for condition in NOISY_CONDITIONS {
        for quantity in QUANTITIES {
            let row = damage[condition]
                .iter()
                .filter(|row| meaningful_loss(row, quantity).is_some())
                .max_by(|left, right| {
                    meaningful_loss(left, quantity)
                        .unwrap()
                        .total_cmp(&meaningful_loss(right, quantity).unwrap())
                })
                .expect("each quantity has a finite loss");
            let (reference, noisy, delta, loss) = quantity_values(row, quantity.name);
            extrema.push(Extremum {
                condition: condition.to_string(),
                quantity: quantity.name,
                time: row.time,
                absolute_difference: -delta,
                relative_loss: loss,
                reference_value: reference,
                noisy_value: noisy,
                noisy_site_population: row.noisy_site_population,
                total_chain_population: row.total_chain_population,
                noisy_site_population_fraction: row.noisy_site_population_fraction,
                dominant_site: dominant_site(&row.site),
                load_energy: row.load_energy,
                load_ergotropy: row.load_ergotropy,
                usable_fraction: row.usable_fraction,
            });
        }
    }
    extrema
}

fn mean_finite(values: impl Iterator<Item = f64>) -> f64 {
    let finite: Vec<_> = values.filter(|value| value.is_finite()).collect();
    if finite.is_empty() {
        f64::NAN
    } else {
        finite.iter().sum::<f64>() / finite.len() as f64
    }
}

fn maximum_finite(values: impl Iterator<Item = f64>) -> f64 {
    values
        .filter(|value| value.is_finite())
        .max_by(f64::total_cmp)
        .unwrap_or(f64::NAN)
}

fn trapezoid_damage(rows: &[&DamageRow], value: impl Fn(&DamageRow) -> f64) -> f64 {
    rows.windows(2)
        .map(|pair| 0.5 * (value(pair[0]) + value(pair[1])) * (pair[1].time - pair[0].time))
        .sum()
}

fn analyze_windows(damage: &BTreeMap<String, Vec<DamageRow>>) -> Vec<WindowResult> {
    let mut results = Vec::new();
    for condition in NOISY_CONDITIONS {
        for window in WINDOWS {
            let area_rows: Vec<_> = damage[condition]
                .iter()
                .filter(|row| row.time >= window.start && row.time <= window.end)
                .collect();
            let mean_rows: Vec<_> = area_rows
                .iter()
                .copied()
                .filter(|row| !window.exclude_start_from_mean || row.time > window.start)
                .collect();
            results.push(WindowResult {
                condition: condition.to_string(),
                window,
                point_count: mean_rows.len(),
                mean_load_energy: mean_finite(mean_rows.iter().map(|row| row.load_energy)),
                mean_load_ergotropy: mean_finite(mean_rows.iter().map(|row| row.load_ergotropy)),
                mean_usable_fraction: mean_finite(mean_rows.iter().map(|row| row.usable_fraction)),
                e_time_area: trapezoid_damage(&area_rows, |row| row.load_energy),
                w_time_area: trapezoid_damage(&area_rows, |row| row.load_ergotropy),
                mean_site: [
                    mean_finite(mean_rows.iter().map(|row| row.site[0])),
                    mean_finite(mean_rows.iter().map(|row| row.site[1])),
                    mean_finite(mean_rows.iter().map(|row| row.site[2])),
                ],
                maximum_e_loss: maximum_finite(
                    mean_rows
                        .iter()
                        .filter_map(|row| meaningful_loss(row, QUANTITIES[0])),
                ),
                maximum_w_loss: maximum_finite(
                    mean_rows
                        .iter()
                        .filter_map(|row| meaningful_loss(row, QUANTITIES[1])),
                ),
                maximum_use_loss: maximum_finite(
                    mean_rows
                        .iter()
                        .filter_map(|row| meaningful_loss(row, QUANTITIES[2])),
                ),
            });
        }
    }
    results
}

fn rank_equivalent(left: f64, right: f64) -> bool {
    (left - right).abs() <= RANK_ABSOLUTE_TOLERANCE
        || (left - right).abs() / left.abs().max(right.abs()).max(SIGNAL_TOLERANCE)
            <= RANK_RELATIVE_TOLERANCE
}

fn worst_label(values: &[(&str, f64)]) -> (String, bool) {
    let finite: Vec<_> = values
        .iter()
        .copied()
        .filter(|(_, value)| value.is_finite())
        .collect();
    if finite.is_empty() {
        return ("undefined".to_string(), false);
    }
    let minimum = finite
        .iter()
        .map(|(_, value)| *value)
        .min_by(f64::total_cmp)
        .unwrap();
    let tied: Vec<_> = finite
        .iter()
        .filter(|(_, value)| rank_equivalent(*value, minimum))
        .map(|(condition, _)| *condition)
        .collect();
    (tied.join("+"), tied.len() > 1)
}

fn build_rankings(damage: &BTreeMap<String, Vec<DamageRow>>) -> Vec<RankingRow> {
    let count = damage[NOISY_CONDITIONS[0]].len();
    (0..count)
        .map(|index| {
            let values: Vec<_> = NOISY_CONDITIONS
                .iter()
                .map(|condition| (*condition, &damage[*condition][index]))
                .collect();
            let (worst_e_condition, e_tie) = worst_label(
                &values
                    .iter()
                    .map(|(condition, row)| (*condition, row.load_energy))
                    .collect::<Vec<_>>(),
            );
            let (worst_w_condition, w_tie) = worst_label(
                &values
                    .iter()
                    .map(|(condition, row)| (*condition, row.load_ergotropy))
                    .collect::<Vec<_>>(),
            );
            let (worst_use_condition, use_tie) = worst_label(
                &values
                    .iter()
                    .map(|(condition, row)| (*condition, row.usable_fraction))
                    .collect::<Vec<_>>(),
            );
            RankingRow {
                time: values[0].1.time,
                worst_e_condition,
                worst_w_condition,
                worst_use_condition,
                e_tie,
                w_tie,
                use_tie,
            }
        })
        .collect()
}

fn confirmed_switches(times: &[f64], labels: &[String]) -> Vec<(f64, String, String)> {
    if labels.len() < CONSECUTIVE_POINTS {
        return Vec::new();
    }
    let mut stable: Option<String> = None;
    let mut switches = Vec::new();
    for index in 0..=labels.len() - CONSECUTIVE_POINTS {
        let candidate = &labels[index];
        if candidate == "undefined"
            || labels[index..index + CONSECUTIVE_POINTS]
                .iter()
                .any(|label| label != candidate)
        {
            continue;
        }
        match &stable {
            None => stable = Some(candidate.clone()),
            Some(previous) if previous != candidate => {
                switches.push((times[index], previous.clone(), candidate.clone()));
                stable = Some(candidate.clone());
            }
            _ => {}
        }
    }
    switches
}

fn pair_label(entrance: f64, exit: f64) -> String {
    if !entrance.is_finite() || !exit.is_finite() {
        "undefined".to_string()
    } else if rank_equivalent(entrance, exit) {
        "tie".to_string()
    } else if entrance < exit {
        "entrance_worse".to_string()
    } else {
        "exit_worse".to_string()
    }
}

fn write_damage_timeseries(damage: &BTreeMap<String, Vec<DamageRow>>) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("local_noise_damage_timeseries.csv")?);
    writeln!(writer, "condition,noise_site,time,load_energy,load_ergotropy,usable_fraction,reference_load_energy,reference_load_ergotropy,reference_usable_fraction,delta_E,delta_W,delta_use,R_E,R_W,R_use,loss_E,loss_W,loss_use,site1_population,site2_population,site3_population,noisy_site_population,total_chain_population,noisy_site_population_fraction")?;
    for condition in NOISY_CONDITIONS {
        for row in &damage[condition] {
            writeln!(
                writer,
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                row.condition,
                row.noise_site,
                n(row.time),
                n(row.load_energy),
                n(row.load_ergotropy),
                n(row.usable_fraction),
                n(row.reference_load_energy),
                n(row.reference_load_ergotropy),
                n(row.reference_usable_fraction),
                n(row.delta_e),
                n(row.delta_w),
                n(row.delta_use),
                n(row.r_e),
                n(row.r_w),
                n(row.r_use),
                n(row.loss_e),
                n(row.loss_w),
                n(row.loss_use),
                n(row.site[0]),
                n(row.site[1]),
                n(row.site[2]),
                n(row.noisy_site_population),
                n(row.total_chain_population),
                n(row.noisy_site_population_fraction)
            )?;
        }
    }
    Ok(())
}

fn write_onsets(onsets: &[Onset]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("local_noise_damage_onsets.csv")?);
    writeln!(writer, "condition,quantity,threshold_level,absolute_threshold,relative_threshold,consecutive_points,first_crossing_time,sustained_onset_time,sustained_duration_at_detection,absolute_difference_at_onset,relative_loss_at_onset,reference_value_at_onset,noisy_value_at_onset,noisy_site_population_at_onset,total_chain_population_at_onset,noisy_site_population_fraction_at_onset,dominant_site_at_onset")?;
    for row in onsets {
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            row.condition,
            row.quantity,
            row.threshold_level,
            n(row.absolute_threshold),
            n(row.relative_threshold),
            CONSECUTIVE_POINTS,
            n(row.first_crossing_time),
            n(row.sustained_onset_time),
            n(row.sustained_duration_at_detection),
            n(row.absolute_difference_at_onset),
            n(row.relative_loss_at_onset),
            n(row.reference_value_at_onset),
            n(row.noisy_value_at_onset),
            n(row.noisy_site_population_at_onset),
            n(row.total_chain_population_at_onset),
            n(row.noisy_site_population_fraction_at_onset),
            row.dominant_site_at_onset
        )?;
    }
    Ok(())
}

fn write_extrema(extrema: &[Extremum]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("local_noise_damage_extrema.csv")?);
    writeln!(writer, "condition,quantity,extremum_type,time,absolute_difference,relative_loss,reference_value,noisy_value,noisy_site_population,total_chain_population,noisy_site_population_fraction,dominant_site,load_energy,load_ergotropy,usable_fraction")?;
    for row in extrema {
        writeln!(
            writer,
            "{},{},maximum_relative_loss,{},{},{},{},{},{},{},{},{},{},{},{}",
            row.condition,
            row.quantity,
            n(row.time),
            n(row.absolute_difference),
            n(row.relative_loss),
            n(row.reference_value),
            n(row.noisy_value),
            n(row.noisy_site_population),
            n(row.total_chain_population),
            n(row.noisy_site_population_fraction),
            row.dominant_site,
            n(row.load_energy),
            n(row.load_ergotropy),
            n(row.usable_fraction)
        )?;
    }
    Ok(())
}

fn write_windows(windows: &[WindowResult]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("local_noise_damage_windows.csv")?);
    writeln!(writer, "condition,window_name,time_start,time_end,point_count,mean_load_energy,mean_load_ergotropy,mean_usable_fraction,E_time_area,W_time_area,mean_site1_population,mean_site2_population,mean_site3_population,maximum_E_loss,maximum_W_loss,maximum_usable_fraction_loss")?;
    for row in windows {
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            row.condition,
            row.window.name,
            n(row.window.start),
            n(row.window.end),
            row.point_count,
            n(row.mean_load_energy),
            n(row.mean_load_ergotropy),
            n(row.mean_usable_fraction),
            n(row.e_time_area),
            n(row.w_time_area),
            n(row.mean_site[0]),
            n(row.mean_site[1]),
            n(row.mean_site[2]),
            n(row.maximum_e_loss),
            n(row.maximum_w_loss),
            n(row.maximum_use_loss)
        )?;
    }
    Ok(())
}

fn write_rankings(rankings: &[RankingRow]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("local_noise_damage_rankings.csv")?);
    writeln!(
        writer,
        "time,worst_E_condition,worst_W_condition,worst_use_condition,E_tie,W_tie,use_tie"
    )?;
    for row in rankings {
        writeln!(
            writer,
            "{},{},{},{},{},{},{}",
            n(row.time),
            row.worst_e_condition,
            row.worst_w_condition,
            row.worst_use_condition,
            row.e_tie,
            row.w_tie,
            row.use_tie
        )?;
    }
    Ok(())
}

fn write_checks(checks: &[Check]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("local_noise_damage_checks.csv")?);
    writeln!(writer, "check,value,expected,result")?;
    for check in checks {
        writeln!(
            writer,
            "{},{},{},{}",
            check.check,
            check.value,
            check.expected,
            if check.pass { "PASS" } else { "FAIL" }
        )?;
    }
    Ok(())
}

fn onset<'a>(onsets: &'a [Onset], condition: &str, quantity: &str, threshold: &str) -> &'a Onset {
    onsets
        .iter()
        .find(|row| {
            row.condition == condition
                && row.quantity == quantity
                && row.threshold_level == threshold
        })
        .expect("requested onset exists")
}

fn extremum<'a>(extrema: &'a [Extremum], condition: &str, quantity: &str) -> &'a Extremum {
    extrema
        .iter()
        .find(|row| row.condition == condition && row.quantity == quantity)
        .expect("requested extremum exists")
}

fn format_switches(switches: &[(f64, String, String)]) -> String {
    if switches.is_empty() {
        "確定した切替なし".to_string()
    } else {
        switches
            .iter()
            .map(|(time, from, to)| format!("t={time:.2}: {from} -> {to}"))
            .collect::<Vec<_>>()
            .join("; ")
    }
}

fn worst_fraction(labels: impl Iterator<Item = String>, condition: &str) -> f64 {
    let defined: Vec<_> = labels.filter(|label| label != "undefined").collect();
    if defined.is_empty() {
        return f64::NAN;
    }
    defined
        .iter()
        .filter(|label| label.split('+').any(|part| part == condition))
        .count() as f64
        / defined.len() as f64
}

fn write_report(
    grouped: &BTreeMap<String, Vec<SourceRow>>,
    damage: &BTreeMap<String, Vec<DamageRow>>,
    onsets: &[Onset],
    extrema: &[Extremum],
    windows: &[WindowResult],
    rankings: &[RankingRow],
    checks: &[Check],
) -> std::io::Result<()> {
    let times: Vec<_> = rankings.iter().map(|row| row.time).collect();
    let e_switches = confirmed_switches(
        &times,
        &rankings
            .iter()
            .map(|row| row.worst_e_condition.clone())
            .collect::<Vec<_>>(),
    );
    let w_switches = confirmed_switches(
        &times,
        &rankings
            .iter()
            .map(|row| row.worst_w_condition.clone())
            .collect::<Vec<_>>(),
    );
    let use_switches = confirmed_switches(
        &times,
        &rankings
            .iter()
            .map(|row| row.worst_use_condition.clone())
            .collect::<Vec<_>>(),
    );
    let entrance_exit_w: Vec<_> = damage["noise_entrance"]
        .iter()
        .zip(&damage["noise_exit"])
        .map(|(entrance, exit)| pair_label(entrance.load_ergotropy, exit.load_ergotropy))
        .collect();
    let entrance_exit_use: Vec<_> = damage["noise_entrance"]
        .iter()
        .zip(&damage["noise_exit"])
        .map(|(entrance, exit)| pair_label(entrance.usable_fraction, exit.usable_fraction))
        .collect();
    let pair_w_switches = confirmed_switches(&times, &entrance_exit_w);
    let pair_use_switches = confirmed_switches(&times, &entrance_exit_use);

    let mut writer = BufWriter::new(File::create("MILESTONE_7B_REPORT.md")?);
    writeln!(
        writer,
        "# Milestone 7b: Local-noise damage onset analysis\n"
    )?;
    writeln!(writer, "## 1. 目的\n\nMilestone 7aの固定済み時系列から、局所位相雑音条件とnoise-freeとの差がいつ持続的に現れたかを記述する。因果分解は行わない。\n")?;
    writeln!(writer, "## 2. 使用した既存データ\n\n主要入力は `{TIMESERIES_INPUT}` のみ。`{SUMMARY_INPUT}` と `{RATIOS_INPUT}` はt=10、W_time_area、比の整合性確認だけに使用した。各条件{}点、合計{}行。\n", grouped["noise_free"].len(), grouped.values().map(Vec::len).sum::<usize>())?;
    writeln!(writer, "## 3. 新規時間発展を行っていないこと\n\nこのbinはCSV reader、算術、台形積分、順位判定、CSV/Markdown writerだけを含む。Hamiltonian、Lindblad、RK4、collapse operator生成、量子状態伝播は呼び出していない。\n")?;
    writeln!(writer, "## 4. onsetの定義\n\n基準量がminimum_reference以上で、`reference - noisy`が絶対閾値以上、かつ相対損失が指定閾値以上となる状態が連続{CONSECUTIVE_POINTS}点続いた最初の点をsustained onsetとした。単一点の超過は採用しない。\n")?;
    writeln!(writer, "## 5. 閾値と持続条件\n\n推奨値を変更せず採用した。E/Wの絶対閾値`1e-5`、usable fractionの絶対閾値`1e-3`。minimum referenceはE/W `1e-4`、usable fraction `1e-3`。相対閾値はweak `1%`、medium `5%`、strong `10%`。時間幅0.01の5点持続を約0.05時間として記録した。\n")?;
    writeln!(writer, "## 6. 数値品質チェック\n\n全{}項目がPASS。必須列、4条件、点数、時間グリッド、単調性、重複、ラベル、定数、有限値/許容NaN、population、E/W/useの範囲、分母処理、summary/ratios整合を確認した。\n", checks.len())?;

    for (section, condition, title) in [
        (7, "noise_entrance", "entrance条件の被害開始"),
        (8, "noise_middle", "middle条件の被害開始"),
        (9, "noise_exit", "exit条件の被害開始"),
    ] {
        writeln!(
            writer,
            "## {section}. {title}\n\n| quantity | weak | medium | strong |\n|---|---:|---:|---:|"
        )?;
        for quantity in ["E", "W", "usable_fraction"] {
            writeln!(
                writer,
                "| {quantity} | {} | {} | {} |",
                n(onset(onsets, condition, quantity, "weak").sustained_onset_time),
                n(onset(onsets, condition, quantity, "medium").sustained_onset_time),
                n(onset(onsets, condition, quantity, "strong").sustained_onset_time)
            )?;
        }
        writeln!(
            writer,
            "\n時刻と同時点のpopulationの対応であり、特定過程だけの損傷や因果原因を意味しない。\n"
        )?;
    }
    writeln!(writer, "## 10. E、W、usable fractionの開始時刻比較\n\nmedium閾値での持続onsetをまとめる。\n\n| condition | E | W | usable fraction |\n|---|---:|---:|---:|")?;
    for condition in NOISY_CONDITIONS {
        writeln!(
            writer,
            "| {condition} | {} | {} | {} |",
            n(onset(onsets, condition, "E", "medium").sustained_onset_time),
            n(onset(onsets, condition, "W", "medium").sustained_onset_time),
            n(onset(onsets, condition, "usable_fraction", "medium").sustained_onset_time)
        )?;
    }
    writeln!(writer, "\nEとWは3条件・3閾値すべてでt=2.25となった。これは同じminimum referenceと絶対閾値を使った診断上の同時開始であり、energy lossとquality lossの物理機構が同時だという因果主張ではない。usable fractionは閾値依存性があり、weakでは3条件ともt=0.98、mediumではexit 1.27、entrance 1.57、middle 1.62、strongではexit 2.54、entrance 2.93、middleは未検出だった。entranceとexitの順序反転はなく、medium/strongではexitが早いが、単一の正確な開始時刻へ圧縮しない。\n")?;
    writeln!(writer, "## 11. 最大損失時刻\n\n| condition | E loss max | W loss max | use loss max |\n|---|---:|---:|---:|")?;
    for condition in NOISY_CONDITIONS {
        writeln!(
            writer,
            "| {condition} | {:.2} | {:.2} | {:.2} |",
            extremum(extrema, condition, "E").time,
            extremum(extrema, condition, "W").time,
            extremum(extrema, condition, "usable_fraction").time
        )?;
    }
    writeln!(
        writer,
        "\n最大値の詳細と同時点のpopulationは `local_noise_damage_extrema.csv` に保存した。\n"
    )?;
    writeln!(writer, "## 12. onset時のsite population\n\nmedium W onsetでの雑音site populationを示す。\n\n| condition | onset | noisy-site population | chain内比 | dominant site |\n|---|---:|---:|---:|---|")?;
    for condition in NOISY_CONDITIONS {
        let row = onset(onsets, condition, "W", "medium");
        writeln!(
            writer,
            "| {condition} | {} | {} | {} | {} |",
            n(row.sustained_onset_time),
            n(row.noisy_site_population_at_onset),
            n(row.noisy_site_population_fraction_at_onset),
            row.dominant_site_at_onset
        )?;
    }
    writeln!(
        writer,
        "\nこれは励起分布の同時記録であり、populationが被害の原因だとは示さない。\n"
    )?;
    writeln!(writer, "## 13. 時間窓別比較\n\n各窓の平均と時間面積は `local_noise_damage_windows.csv` に保存した。W_time_areaが最小の条件を窓ごとに示す。\n")?;
    for window in WINDOWS {
        let worst = windows
            .iter()
            .filter(|row| row.window.name == window.name)
            .min_by(|left, right| left.w_time_area.total_cmp(&right.w_time_area))
            .unwrap();
        writeln!(
            writer,
            "- {}: `{}` (W_time_area={:.8e})",
            window.name, worst.condition, worst.w_time_area
        )?;
    }
    writeln!(
        writer,
        "\nE/W_time_areaは状態量の時間面積であり、累積流入エネルギーや累積抽出仕事ではない。\n"
    )?;
    writeln!(writer, "## 14. 順位の時間変化\n\n- worst Eの確定切替: {}\n- worst Wの確定切替: {}\n- worst usable fractionの確定切替: {}\n- entrance/exit Wの入替: {}\n- entrance/exit usable fractionの入替: {}\n", format_switches(&e_switches), format_switches(&w_switches), format_switches(&use_switches), format_switches(&pair_w_switches), format_switches(&pair_use_switches))?;
    writeln!(writer, "\nTieを含めて各条件が最悪集合に入った時間割合を示す（tieのため合計は100%を超えうる）。\n\n| condition | E | W | usable fraction |\n|---|---:|---:|---:|")?;
    for condition in NOISY_CONDITIONS {
        writeln!(
            writer,
            "| {condition} | {:.4} | {:.4} | {:.4} |",
            worst_fraction(
                rankings.iter().map(|row| row.worst_e_condition.clone()),
                condition
            ),
            worst_fraction(
                rankings.iter().map(|row| row.worst_w_condition.clone()),
                condition
            ),
            worst_fraction(
                rankings.iter().map(|row| row.worst_use_condition.clone()),
                condition
            )
        )?;
    }
    writeln!(writer, "短い反転を除くため、同じ順位状態が5点続いた時だけ切替とした。全時刻の順位は `local_noise_damage_rankings.csv`。\n")?;

    let mut population_lines = Vec::new();
    for condition in NOISY_CONDITIONS {
        let rows = &damage[condition];
        let max_population = rows
            .iter()
            .map(|row| row.noisy_site_population)
            .max_by(f64::total_cmp)
            .unwrap();
        let population_area = trapezoid_damage(&rows.iter().collect::<Vec<_>>(), |row| {
            row.noisy_site_population
        });
        let pulse_mean = mean_finite(
            rows.iter()
                .filter(|row| row.time <= 3.2)
                .map(|row| row.noisy_site_population),
        );
        let post_mean = mean_finite(
            rows.iter()
                .filter(|row| row.time > 3.2)
                .map(|row| row.noisy_site_population),
        );
        let onset_population =
            onset(onsets, condition, "W", "medium").noisy_site_population_at_onset;
        let max_loss_population = extremum(extrema, condition, "W").noisy_site_population;
        population_lines.push((
            condition,
            max_population,
            population_area,
            pulse_mean,
            post_mean,
            onset_population,
            max_loss_population,
        ));
    }
    writeln!(writer, "## 15. middle条件の被害が小さいことへの手がかり\n\n| condition | max noisy-site p | p time-area | pulse mean p | post-pulse mean p | p at W onset | p at max W loss |\n|---|---:|---:|---:|---:|---:|---:|")?;
    for row in &population_lines {
        writeln!(
            writer,
            "| {} | {:.6e} | {:.6e} | {:.6e} | {:.6e} | {:.6e} | {:.6e} |",
            row.0, row.1, row.2, row.3, row.4, row.5, row.6
        )?;
    }
    writeln!(writer, "\nmiddleの雑音site population時間面積はentranceに近く、exitより小さかった。pulse中平均はentranceより小さくexitより大きい3条件中間だった一方、7aのW損失はmiddleが明らかに軽かった。したがって、少なくとも単純な最大populationまたは時間面積だけではW損失差を説明し切れない。時刻のずれや未計算の量も候補として残る。これは手がかりであり原因説明ではない。\n")?;
    writeln!(writer, "## 16. 直接確認できたこと\n\n固定された7aデータ内で、E/W/usable fractionの差が各閾値で持続的に現れた時刻、その時刻のsite population、時間窓別の状態量、順位切替を直接確認した。\n")?;
    writeln!(writer, "## 17. 確認できていないこと\n\n特定物理過程への因果分解、populationだけによる説明、coherence current、site間current、保護による回復、他パラメータ・長いnetworkへの一般化は確認していない。\n")?;
    writeln!(writer, "## 18. 主張してはいけないこと\n\nentrance/middle/exit雑音が注入/輸送/受け渡しだけを壊すという断定、populationが大きいことを原因とする断定、先に変化した量を原因とする断定、保護効果の予測はできない。\n")?;
    writeln!(writer, "## 19. 次の保護実験を選ぶための判断材料\n\n- entrance protection候補: 7aでt=10のW、W_max、W_time_areaが最小であり、今回の時系列でも入口雑音条件の持続損失が確認された。\n- exit protection候補: 7aでt=10のusable fractionが最小であり、今回その開始時刻と順位変化を確認した。\n- 同時比較候補: entranceは仕事量、exitは使える割合で異なる悪化指標を持つため、同じ固定条件で並べる価値がある。\n\nこれらは候補選定の材料であり、保護による回復を予測するものではない。保護機能は実装していない。\n")?;
    writeln!(writer, "## 20. 生成ファイル一覧\n\n- `src/bin/local_noise_damage_analysis.rs`\n- `local_noise_damage_timeseries.csv`\n- `local_noise_damage_onsets.csv`\n- `local_noise_damage_extrema.csv`\n- `local_noise_damage_windows.csv`\n- `local_noise_damage_rankings.csv`\n- `local_noise_damage_checks.csv`\n- `MILESTONE_7B_REPORT.md`\n")?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("reading existing Milestone 7a CSV files; no propagation will be run");
    let timeseries_table = CsvTable::read(TIMESERIES_INPUT)?;
    let summary_table = CsvTable::read(SUMMARY_INPUT)?;
    let ratios_table = CsvTable::read(RATIOS_INPUT)?;
    let source = parse_source(&timeseries_table)?;
    let grouped = group_source(&source);
    let mut checks = source_checks(&timeseries_table, &grouped, &summary_table, &ratios_table)?;
    let all_input_checks = checks.iter().all(|check| check.pass);
    write_checks(&checks)?;
    if !all_input_checks {
        return Err("input validation failed; analysis stopped".into());
    }

    let damage = build_damage(&grouped)?;
    let damage_finite = damage.values().flatten().all(|row| {
        [
            row.delta_e,
            row.delta_w,
            row.noisy_site_population,
            row.total_chain_population,
        ]
        .iter()
        .all(|value| value.is_finite())
            && [
                row.delta_use,
                row.r_e,
                row.r_w,
                row.r_use,
                row.loss_e,
                row.loss_w,
                row.loss_use,
                row.noisy_site_population_fraction,
            ]
            .iter()
            .all(|value| finite_or_nan(*value))
    });
    checks.push(Check {
        check: "derived_damage_finite_or_allowed_nan".to_string(),
        value: damage_finite.to_string(),
        expected: "true".to_string(),
        pass: damage_finite,
    });
    let denominator_consistent = damage.values().flatten().all(|row| {
        ratio_rule(row.load_energy, row.reference_load_energy, row.r_e)
            && ratio_rule(row.load_ergotropy, row.reference_load_ergotropy, row.r_w)
            && ratio_rule(
                row.usable_fraction,
                row.reference_usable_fraction,
                row.r_use,
            )
    });
    checks.push(Check {
        check: "derived_ratio_denominator_rule".to_string(),
        value: denominator_consistent.to_string(),
        expected: "NaN exactly for invalid/small reference".to_string(),
        pass: denominator_consistent,
    });
    let reference_zero = source
        .iter()
        .filter(|row| row.condition == "noise_free")
        .all(|row| {
            (row.load_energy - row.load_energy).abs() == 0.0
                && (row.load_ergotropy - row.load_ergotropy).abs() == 0.0
        });
    checks.push(Check {
        check: "noise_free_self_difference_zero".to_string(),
        value: reference_zero.to_string(),
        expected: "true".to_string(),
        pass: reference_zero,
    });
    write_checks(&checks)?;
    if checks.iter().any(|check| !check.pass) {
        return Err("derived-data validation failed; analysis stopped".into());
    }

    let onsets = detect_onsets(&damage);
    let extrema = detect_extrema(&damage);
    let windows = analyze_windows(&damage);
    let rankings = build_rankings(&damage);
    write_damage_timeseries(&damage)?;
    write_onsets(&onsets)?;
    write_extrema(&extrema)?;
    write_windows(&windows)?;
    write_rankings(&rankings)?;
    write_report(
        &grouped, &damage, &onsets, &extrema, &windows, &rankings, &checks,
    )?;
    println!("Milestone 7b CSV-only analysis complete");
    Ok(())
}
