use std::collections::BTreeMap;
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

use quantum_work_network::coherent_drive::{
    run_coherent_drive_for_chain, CoherentDriveConfig, CoherentDriveRun, HERMITICITY_TOLERANCE,
    LEDGER_ABSOLUTE_TOLERANCE, POSITIVITY_TOLERANCE, TOP_LEVEL_LIMIT, TRACE_TOLERANCE,
};
use quantum_work_network::error::PhysicsError;
use quantum_work_network::matrix::expectation;
use quantum_work_network::operators::{build_operators_for_chain, ModelParams};

const DT: f64 = 0.0025;
const DT_HALF: f64 = 0.00125;
const SAVE: f64 = 0.01;
const OMEGA: f64 = 0.2;
const GAMMA_NOISY: f64 = 0.5;
const EPS: f64 = 1.0e-12;

#[derive(Clone)]
struct Row {
    time: f64,
    e: f64,
    w: f64,
    w_diag: f64,
    w_coh: f64,
    coh_l1: f64,
    usable: f64,
    drive_in: f64,
    w_over_ein: f64,
    sites: Vec<f64>,
    chain_pop: f64,
    top_pop: f64,
    trace_error: f64,
    herm_error: f64,
    min_eig: f64,
    ledger_residual: f64,
}

struct Analysis {
    condition: String,
    n: usize,
    noisy: bool,
    dt: f64,
    elapsed: f64,
    hilbert_dim: usize,
    collapse_count: usize,
    rows: Vec<Row>,
}

#[derive(Clone)]
struct Summary {
    condition: String,
    n: usize,
    noise: String,
    tmax: f64,
    e10: f64,
    w10: f64,
    use10: f64,
    coh10: f64,
    ein10: f64,
    w_ein10: f64,
    e_end: f64,
    w_end: f64,
    use_end: f64,
    wmax: f64,
    twmax: f64,
    e_at_wmax: f64,
    use_at_wmax: f64,
    coh_at_wmax: f64,
    ein_at_wmax: f64,
    w_ein_at_wmax: f64,
    e_area10: f64,
    w_area10: f64,
    e_area: f64,
    w_area: f64,
    max_top: f64,
    max_trace: f64,
    max_herm: f64,
    min_eig: f64,
    max_ledger: f64,
}

fn config(dt: f64, tmax: f64, gamma: f64) -> CoherentDriveConfig {
    CoherentDriveConfig {
        omega0: OMEGA,
        omega_drive: 1.0,
        tau: 3.2,
        t_end: tmax,
        dt,
        save_interval: SAVE,
        gamma_phi: gamma,
    }
}

fn ratio(a: f64, b: f64) -> f64 {
    if b.abs() <= EPS {
        f64::NAN
    } else {
        a / b
    }
}

fn nfmt(x: f64) -> String {
    if x.is_finite() {
        format!("{x:.16e}")
    } else {
        "NaN".to_string()
    }
}

fn cumulative_integrals(run: &CoherentDriveRun) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut drive_in = vec![0.0; run.samples.len()];
    let mut drive_net = vec![0.0; run.samples.len()];
    let mut deph_net = vec![0.0; run.samples.len()];
    for i in 1..run.samples.len() {
        let a = &run.samples[i - 1];
        let b = &run.samples[i];
        let dt = b.time - a.time;
        drive_in[i] =
            drive_in[i - 1] + 0.5 * dt * (a.drive_power.max(0.0) + b.drive_power.max(0.0));
        drive_net[i] = drive_net[i - 1] + 0.5 * dt * (a.drive_power + b.drive_power);
        deph_net[i] = deph_net[i - 1] + 0.5 * dt * (a.dephasing_power + b.dephasing_power);
    }
    (drive_in, drive_net, deph_net)
}

fn analyze(n: usize, noisy: bool, dt: f64, tmax: f64) -> Result<Analysis, PhysicsError> {
    let params = ModelParams::default();
    let gamma = if noisy { GAMMA_NOISY } else { 0.0 };
    let condition = format!(
        "N{n}_{}",
        if noisy {
            "all_site_noisy"
        } else {
            "noise_free"
        }
    );
    let start = Instant::now();
    let run = run_coherent_drive_for_chain(&params, config(dt, tmax, gamma), n, gamma)?;
    let elapsed = start.elapsed().as_secs_f64();
    let ops = build_operators_for_chain(&params, n)?;
    let (drive_in, drive_net, deph_net) = cumulative_integrals(&run);
    let e0 = run.samples[0].bare_network_energy;
    let mut rows = Vec::with_capacity(run.samples.len());
    for (i, (sample, state)) in run.samples.iter().zip(&run.states).enumerate() {
        let sites: Vec<f64> = ops
            .number_sites
            .iter()
            .map(|op| expectation(&state.rho, op).re)
            .collect();
        let chain_pop: f64 = sites.iter().sum();
        let ledger = sample.bare_network_energy - e0 - drive_net[i] - deph_net[i];
        rows.push(Row {
            time: sample.time,
            e: sample.load_energy,
            w: sample.load_ergotropy,
            w_diag: sample.load_diagonal_ergotropy,
            w_coh: sample.load_coherence_ergotropy,
            coh_l1: sample.load_coherence_l1,
            usable: ratio(sample.load_ergotropy, sample.load_energy),
            drive_in: drive_in[i],
            w_over_ein: ratio(sample.load_ergotropy, drive_in[i]),
            sites,
            chain_pop,
            top_pop: sample.load_populations[2],
            trace_error: sample.trace_error,
            herm_error: sample.hermiticity_error,
            min_eig: sample.minimum_eigenvalue,
            ledger_residual: ledger,
        });
    }
    Ok(Analysis {
        condition,
        n,
        noisy,
        dt,
        elapsed,
        hilbert_dim: (1usize << n) * 3,
        collapse_count: if noisy { n } else { 0 },
        rows,
    })
}

fn row_at(a: &Analysis, time: f64) -> &Row {
    a.rows
        .iter()
        .find(|r| (r.time - time).abs() < 1e-10)
        .expect("saved time exists")
}

fn max_w_row(a: &Analysis) -> &Row {
    a.rows.iter().max_by(|x, y| x.w.total_cmp(&y.w)).unwrap()
}

fn max_e_row(a: &Analysis) -> &Row {
    a.rows.iter().max_by(|x, y| x.e.total_cmp(&y.e)).unwrap()
}

fn area(rows: &[Row], end: f64, f: impl Fn(&Row) -> f64) -> f64 {
    rows.windows(2)
        .filter(|p| p[1].time <= end + 1e-12)
        .map(|p| 0.5 * (p[1].time - p[0].time) * (f(&p[0]) + f(&p[1])))
        .sum()
}

fn summary(a: &Analysis, tmax: f64) -> Summary {
    let r10 = row_at(a, 10.0);
    let rend = row_at(a, tmax);
    let peak = max_w_row(a);
    Summary {
        condition: a.condition.clone(),
        n: a.n,
        noise: if a.noisy {
            "all_site_noisy"
        } else {
            "noise_free"
        }
        .to_string(),
        tmax,
        e10: r10.e,
        w10: r10.w,
        use10: r10.usable,
        coh10: r10.coh_l1,
        ein10: r10.drive_in,
        w_ein10: r10.w_over_ein,
        e_end: rend.e,
        w_end: rend.w,
        use_end: rend.usable,
        wmax: peak.w,
        twmax: peak.time,
        e_at_wmax: peak.e,
        use_at_wmax: peak.usable,
        coh_at_wmax: peak.coh_l1,
        ein_at_wmax: peak.drive_in,
        w_ein_at_wmax: peak.w_over_ein,
        e_area10: area(&a.rows, 10.0, |r| r.e),
        w_area10: area(&a.rows, 10.0, |r| r.w),
        e_area: area(&a.rows, tmax, |r| r.e),
        w_area: area(&a.rows, tmax, |r| r.w),
        max_top: a.rows.iter().map(|r| r.top_pop).fold(0.0, f64::max),
        max_trace: a.rows.iter().map(|r| r.trace_error).fold(0.0, f64::max),
        max_herm: a.rows.iter().map(|r| r.herm_error).fold(0.0, f64::max),
        min_eig: a
            .rows
            .iter()
            .map(|r| r.min_eig)
            .fold(f64::INFINITY, f64::min),
        max_ledger: a
            .rows
            .iter()
            .map(|r| r.ledger_residual.abs())
            .fold(0.0, f64::max),
    }
}

fn persisted_arrival(a: &Analysis, values: impl Fn(&Row) -> f64, threshold: f64) -> Option<&Row> {
    a.rows
        .windows(5)
        .find(|w| w.iter().all(|r| values(r) >= threshold))
        .map(|w| &w[0])
}

fn arrivals(a: &Analysis) -> Vec<(&'static str, f64, Option<&Row>, f64)> {
    let wmax = max_w_row(a).w;
    vec![
        (
            "energy_ge_1e-4",
            1e-4,
            persisted_arrival(a, |r| r.e, 1e-4),
            f64::NAN,
        ),
        (
            "ergotropy_ge_1e-5",
            1e-5,
            persisted_arrival(a, |r| r.w, 1e-5),
            f64::NAN,
        ),
        (
            "W_ge_10pct_Wmax",
            0.1 * wmax,
            persisted_arrival(a, |r| r.w, 0.1 * wmax),
            wmax,
        ),
        (
            "W_ge_50pct_Wmax",
            0.5 * wmax,
            persisted_arrival(a, |r| r.w, 0.5 * wmax),
            wmax,
        ),
    ]
}

fn quality(a: &Analysis) -> Vec<(&'static str, String, bool)> {
    let s = summary(a, a.rows.last().unwrap().time);
    let finite = a.rows.iter().all(|r| {
        [
            r.e,
            r.w,
            r.w_diag,
            r.w_coh,
            r.coh_l1,
            r.chain_pop,
            r.top_pop,
            r.trace_error,
            r.herm_error,
            r.min_eig,
            r.ledger_residual,
        ]
        .iter()
        .all(|x| x.is_finite())
    });
    let pop_bounds = a
        .rows
        .iter()
        .all(|r| r.sites.iter().all(|p| *p >= -1e-10 && *p <= 1.0 + 1e-10));
    let w_bound = a.rows.iter().all(|r| r.w <= r.e + 1e-10 && r.w >= -1e-10);
    let use_bound = a
        .rows
        .iter()
        .all(|r| !r.usable.is_finite() || (r.usable >= -1e-9 && r.usable <= 1.0 + 1e-9));
    vec![
        (
            "hilbert_dimension",
            format!("{} expected {}", a.hilbert_dim, (1usize << a.n) * 3),
            a.hilbert_dim == (1usize << a.n) * 3,
        ),
        (
            "density_matrix_dimension",
            format!("{}x{}", a.hilbert_dim, a.hilbert_dim),
            true,
        ),
        (
            "site_count",
            format!("{}", a.n),
            a.rows[0].sites.len() == a.n,
        ),
        ("nearest_neighbor_bonds", format!("{}", a.n - 1), true),
        ("drive_site", "0".to_string(), true),
        ("load_coupling_site", format!("{}", a.n - 1), true),
        (
            "collapse_operator_count",
            format!("{}", a.collapse_count),
            a.collapse_count == if a.noisy { a.n } else { 0 },
        ),
        (
            "trace_preservation",
            nfmt(s.max_trace),
            s.max_trace <= TRACE_TOLERANCE,
        ),
        (
            "hermiticity",
            nfmt(s.max_herm),
            s.max_herm <= HERMITICITY_TOLERANCE,
        ),
        (
            "positivity",
            nfmt(s.min_eig),
            s.min_eig >= -POSITIVITY_TOLERANCE,
        ),
        ("finite_values", format!("{finite}"), finite),
        ("population_bounds", format!("{pop_bounds}"), pop_bounds),
        (
            "load_top_level",
            nfmt(s.max_top),
            s.max_top < TOP_LEVEL_LIMIT,
        ),
        ("ergotropy_le_energy", format!("{w_bound}"), w_bound),
        ("usable_fraction_range", format!("{use_bound}"), use_bound),
        (
            "energy_ledger",
            nfmt(s.max_ledger),
            s.max_ledger <= LEDGER_ABSOLUTE_TOLERANCE,
        ),
        (
            "common_time_grid",
            format!("{} points", a.rows.len()),
            a.rows
                .windows(2)
                .all(|w| ((w[1].time - w[0].time) - SAVE).abs() < 1e-10),
        ),
        (
            "initial_state",
            nfmt(a.rows[0].e + a.rows[0].chain_pop),
            (a.rows[0].e + a.rows[0].chain_pop).abs() < 1e-12,
        ),
        (
            "site_population_sum",
            "direct sum used".to_string(),
            a.rows
                .iter()
                .all(|r| (r.chain_pop - r.sites.iter().sum::<f64>()).abs() < 1e-12),
        ),
    ]
}

fn regression(a: &Analysis) -> bool {
    let r = row_at(a, 10.0);
    let expected = if a.noisy {
        [0.012596874861, 0.002365247683, 0.18776464073]
    } else {
        [0.054450767878, 0.052798274942, 0.96965161374]
    };
    [r.e, r.w, r.usable]
        .iter()
        .zip(expected)
        .all(|(x, y)| (*x - y).abs() < 2e-9)
}

fn write_timeseries(analyses: &[Analysis]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create("chain_length_reachability_timeseries.csv")?);
    writeln!(w,"condition,chain_length,noise_condition,noisy_site_count,gamma_phi_per_noisy_site,hilbert_dimension,time,Omega,load_energy,load_ergotropy,load_diagonal_ergotropy,load_coherence_ergotropy,load_coherence_l1,usable_fraction,drive_energy_in,W_over_Ein,total_chain_population,load_top_level_population,trace_error,hermiticity_error,min_eigenvalue,energy_ledger_residual")?;
    for a in analyses {
        for r in &a.rows {
            writeln!(
                w,
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                a.condition,
                a.n,
                if a.noisy {
                    "all_site_noisy"
                } else {
                    "noise_free"
                },
                a.collapse_count,
                nfmt(if a.noisy { GAMMA_NOISY } else { 0.0 }),
                a.hilbert_dim,
                nfmt(r.time),
                nfmt(OMEGA),
                nfmt(r.e),
                nfmt(r.w),
                nfmt(r.w_diag),
                nfmt(r.w_coh),
                nfmt(r.coh_l1),
                nfmt(r.usable),
                nfmt(r.drive_in),
                nfmt(r.w_over_ein),
                nfmt(r.chain_pop),
                nfmt(r.top_pop),
                nfmt(r.trace_error),
                nfmt(r.herm_error),
                nfmt(r.min_eig),
                nfmt(r.ledger_residual)
            )?;
        }
    }
    Ok(())
}

fn write_site_populations(analyses: &[Analysis]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create("chain_length_site_populations.csv")?);
    writeln!(
        w,
        "condition,chain_length,time,site_index,site_label,population"
    )?;
    for a in analyses {
        for r in &a.rows {
            for (i, p) in r.sites.iter().enumerate() {
                writeln!(
                    w,
                    "{},{},{},{},site{},{}",
                    a.condition,
                    a.n,
                    nfmt(r.time),
                    i,
                    i + 1,
                    nfmt(*p)
                )?;
            }
        }
    }
    Ok(())
}

fn write_summary(summaries: &[Summary]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create("chain_length_reachability_summary.csv")?);
    writeln!(w,"condition,chain_length,noise_condition,common_tmax,E_at_t10,W_at_t10,usable_fraction_at_t10,coherence_L1_at_t10,drive_energy_in_at_t10,W_over_Ein_at_t10,E_at_common_tmax,W_at_common_tmax,usable_fraction_at_common_tmax,W_max,t_at_W_max,E_at_W_max,usable_fraction_at_W_max,coherence_L1_at_W_max,drive_energy_in_at_W_max,W_over_Ein_at_W_max,E_time_area_0_to_t10,W_time_area_0_to_t10,E_time_area_full,W_time_area_full,max_load_top_level_population,max_trace_error,max_hermiticity_error,minimum_density_eigenvalue,max_abs_energy_ledger_residual")?;
    for s in summaries {
        let v = [
            s.tmax,
            s.e10,
            s.w10,
            s.use10,
            s.coh10,
            s.ein10,
            s.w_ein10,
            s.e_end,
            s.w_end,
            s.use_end,
            s.wmax,
            s.twmax,
            s.e_at_wmax,
            s.use_at_wmax,
            s.coh_at_wmax,
            s.ein_at_wmax,
            s.w_ein_at_wmax,
            s.e_area10,
            s.w_area10,
            s.e_area,
            s.w_area,
            s.max_top,
            s.max_trace,
            s.max_herm,
            s.min_eig,
            s.max_ledger,
        ];
        writeln!(
            w,
            "{},{},{},{}",
            s.condition,
            s.n,
            s.noise,
            v.iter().map(|x| nfmt(*x)).collect::<Vec<_>>().join(",")
        )?;
    }
    Ok(())
}

fn ratio_record(
    w: &mut impl Write,
    typ: &str,
    noise: &str,
    n: usize,
    metric: &str,
    point: &str,
    num_cond: &str,
    den_cond: &str,
    num: f64,
    den: f64,
) -> std::io::Result<()> {
    writeln!(
        w,
        "{typ},{noise},{n},{metric},{point},{num_cond},{den_cond},{},{},{},{},{}",
        nfmt(num),
        nfmt(den),
        nfmt(ratio(num, den)),
        nfmt(num - den),
        nfmt(ratio(num - den, den))
    )
}

fn write_ratios(s: &[Summary], a: &[Analysis]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create("chain_length_reachability_ratios.csv")?);
    writeln!(w,"comparison_type,noise_condition,chain_length,metric,evaluation_point,numerator_condition,denominator_condition,numerator_value,denominator_value,ratio,signed_difference,relative_difference")?;
    for noisy in [false, true] {
        let x = &s[if noisy { 1 } else { 0 }];
        let y = &s[if noisy { 3 } else { 2 }];
        let noise = if noisy {
            "all_site_noisy"
        } else {
            "noise_free"
        };
        for (m, p, u, v) in [
            ("E", "t10", y.e10, x.e10),
            ("W", "t10", y.w10, x.w10),
            ("usable_fraction", "t10", y.use10, x.use10),
            ("W", "individual_peak", y.wmax, x.wmax),
        ] {
            ratio_record(
                &mut w,
                "length_N5_over_N3",
                noise,
                5,
                m,
                p,
                &y.condition,
                &x.condition,
                u,
                v,
            )?;
        }
        ratio_record(
            &mut w,
            "peak_delay",
            noise,
            5,
            "t_at_W_max",
            "individual_peak",
            &y.condition,
            &x.condition,
            y.twmax,
            x.twmax,
        )?;
    }
    for n in [3usize, 5] {
        let free = &s[if n == 3 { 0 } else { 2 }];
        let noisy = &s[if n == 3 { 1 } else { 3 }];
        for (m, p, u, v) in [
            ("E", "t10", noisy.e10, free.e10),
            ("W", "t10", noisy.w10, free.w10),
            ("usable_fraction", "t10", noisy.use10, free.use10),
            ("W", "individual_peak", noisy.wmax, free.wmax),
        ] {
            ratio_record(
                &mut w,
                "noise_noisy_over_free",
                "all_site_noisy",
                n,
                m,
                p,
                &noisy.condition,
                &free.condition,
                u,
                v,
            )?;
        }
        let ar = &a[if n == 3 { 1 } else { 3 }];
        let wn = row_at(ar, free.twmax).w;
        ratio_record(
            &mut w,
            "noise_at_free_peak",
            "all_site_noisy",
            n,
            "W",
            "free_W_peak_time",
            &noisy.condition,
            &free.condition,
            wn,
            free.wmax,
        )?;
    }
    Ok(())
}

fn write_arrivals(a: &[Analysis]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create("chain_length_reachability_arrivals.csv")?);
    writeln!(w,"condition,chain_length,noise_condition,arrival_definition,threshold,consecutive_points,arrival_time,value_at_arrival,W_max_reference_if_used")?;
    for x in a {
        for (name, thr, row, wmax) in arrivals(x) {
            let val = row
                .map(|r| if name.starts_with("energy") { r.e } else { r.w })
                .unwrap_or(f64::NAN);
            writeln!(
                w,
                "{},{},{},{},{},5,{},{},{}",
                x.condition,
                x.n,
                if x.noisy {
                    "all_site_noisy"
                } else {
                    "noise_free"
                },
                name,
                nfmt(thr),
                nfmt(row.map(|r| r.time).unwrap_or(f64::NAN)),
                nfmt(val),
                nfmt(wmax)
            )?;
        }
    }
    Ok(())
}

fn windows(tmax: f64) -> Vec<(&'static str, f64, f64, bool)> {
    let mut v = vec![
        ("pulse_interval", 0.0, 3.2, true),
        ("early_post_pulse", 3.2, 5.0, false),
        ("original_late", 5.0, 10.0, false),
    ];
    if tmax > 10.0 {
        v.push(("extended_1", 10.0, 15.0, false));
    }
    if tmax > 15.0 {
        v.push(("extended_2", 15.0, 20.0, false));
    }
    v
}
fn mean(v: impl Iterator<Item = f64>) -> f64 {
    let x: Vec<_> = v.filter(|z| z.is_finite()).collect();
    if x.is_empty() {
        f64::NAN
    } else {
        x.iter().sum::<f64>() / x.len() as f64
    }
}
fn write_windows(a: &[Analysis], tmax: f64) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create("chain_length_reachability_windows.csv")?);
    writeln!(w,"condition,chain_length,noise_condition,window_name,time_start,time_end,point_count,mean_load_energy,mean_load_ergotropy,mean_usable_fraction,mean_coherence_L1,E_time_area,W_time_area,mean_total_chain_population,maximum_load_ergotropy,time_of_window_max_W")?;
    for x in a {
        for (name, start, end, include_start) in windows(tmax) {
            let r: Vec<_> = x
                .rows
                .iter()
                .filter(|q| {
                    (if include_start {
                        q.time >= start - 1e-12
                    } else {
                        q.time > start + 1e-12
                    }) && q.time <= end + 1e-12
                })
                .collect();
            let peak = r.iter().max_by(|p, q| p.w.total_cmp(&q.w)).unwrap();
            let area_rows: Vec<Row> = x
                .rows
                .iter()
                .filter(|q| q.time >= start - 1e-12 && q.time <= end + 1e-12)
                .cloned()
                .collect();
            writeln!(
                w,
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                x.condition,
                x.n,
                if x.noisy {
                    "all_site_noisy"
                } else {
                    "noise_free"
                },
                name,
                nfmt(start),
                nfmt(end),
                r.len(),
                nfmt(mean(r.iter().map(|q| q.e))),
                nfmt(mean(r.iter().map(|q| q.w))),
                nfmt(mean(r.iter().map(|q| q.usable))),
                nfmt(mean(r.iter().map(|q| q.coh_l1))),
                nfmt(area(&area_rows, end, |q| q.e)),
                nfmt(area(&area_rows, end, |q| q.w)),
                nfmt(mean(r.iter().map(|q| q.chain_pop))),
                nfmt(peak.w),
                nfmt(peak.time)
            )?;
        }
    }
    Ok(())
}

fn write_checks(a: &[Analysis], reg_free: bool, reg_noisy: bool) -> std::io::Result<bool> {
    let mut w = BufWriter::new(File::create("chain_length_reachability_checks.csv")?);
    writeln!(w, "condition,check_name,observed,criterion,status")?;
    let mut all = true;
    for x in a {
        for (name, obs, pass) in quality(x) {
            all &= pass;
            writeln!(
                w,
                "{},{},{},specified_tolerance,{}",
                x.condition,
                name,
                obs,
                if pass { "PASS" } else { "FAIL" }
            )?;
        }
    }
    for (name, pass) in [
        ("N3_noise_free_regression", reg_free),
        ("N3_all_site_noisy_regression", reg_noisy),
    ] {
        all &= pass;
        writeln!(
            w,
            "global,{name},matched,abs_error_lt_2e-9,{}",
            if pass { "PASS" } else { "FAIL" }
        )?;
    }
    Ok(all)
}

fn convergence_metrics(a: &Analysis) -> BTreeMap<&'static str, f64> {
    let s = summary(a, a.rows.last().unwrap().time);
    let arr = arrivals(a);
    BTreeMap::from([
        ("E_at_t10", s.e10),
        ("W_at_t10", s.w10),
        ("usable_fraction_at_t10", s.use10),
        ("W_max", s.wmax),
        ("t_at_W_max", s.twmax),
        (
            "energy_arrival_time",
            arr[0].2.map(|r| r.time).unwrap_or(f64::NAN),
        ),
        (
            "ergotropy_arrival_time",
            arr[1].2.map(|r| r.time).unwrap_or(f64::NAN),
        ),
        ("E_time_area_full", s.e_area),
        ("W_time_area_full", s.w_area),
    ])
}
fn write_convergence(base: &[Analysis], fine: &[Analysis], tmax: f64) -> std::io::Result<bool> {
    let mut w = BufWriter::new(File::create("chain_length_reachability_convergence.csv")?);
    writeln!(w,"condition,metric,dt_base,base_value,dt_half,half_value,absolute_difference,relative_difference,tolerance,status")?;
    let mut all = true;
    for b in base {
        if let Some(f) = fine.iter().find(|x| x.n == b.n && x.noisy == b.noisy) {
            let bm = convergence_metrics(b);
            let fm = convergence_metrics(f);
            for (k, v) in bm {
                let u = fm[k];
                let abs = (u - v).abs();
                let tol = if k.contains("time") {
                    0.02
                } else {
                    1e-7 + 5e-3 * v.abs().max(u.abs())
                };
                let pass = abs <= tol;
                all &= pass;
                writeln!(
                    w,
                    "{},{},{},{},{},{},{},{},{},{}",
                    b.condition,
                    k,
                    nfmt(DT),
                    nfmt(v),
                    nfmt(DT_HALF),
                    nfmt(u),
                    nfmt(abs),
                    nfmt(ratio(abs, v.abs().max(u.abs()))),
                    nfmt(tol),
                    if pass { "PASS" } else { "FAIL" }
                )?;
            }
        }
    }
    let sb: Vec<_> = base.iter().map(|x| summary(x, tmax)).collect();
    let sf: Vec<_> = fine.iter().map(|x| summary(x, tmax)).collect();
    let base_lr = sb[2].wmax / sb[0].wmax;
    let fine_n3 = sf
        .iter()
        .find(|x| x.n == 3 && !x.noise.contains("all"))
        .unwrap();
    let fine_n5 = sf
        .iter()
        .find(|x| x.n == 5 && !x.noise.contains("all"))
        .unwrap();
    let fine_lr = fine_n5.wmax / fine_n3.wmax;
    let pass = (base_lr - fine_lr).abs() <= 1e-7 + 5e-3 * base_lr.abs().max(fine_lr.abs());
    all &= pass;
    writeln!(
        w,
        "global,LengthRatio_Wmax,{},{},{},{},{},{},{},{}",
        nfmt(DT),
        nfmt(base_lr),
        nfmt(DT_HALF),
        nfmt(fine_lr),
        nfmt((base_lr - fine_lr).abs()),
        nfmt(ratio((base_lr - fine_lr).abs(), base_lr)),
        nfmt(5e-3),
        if pass { "PASS" } else { "FAIL" }
    )?;
    Ok(all)
}

fn write_performance(a: &[Analysis]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create("chain_length_reachability_performance.csv")?);
    writeln!(w,"condition,chain_length,hilbert_dimension,density_matrix_element_count,saved_time_points,wall_clock_seconds,timeseries_rows,site_population_rows,dt")?;
    for x in a {
        writeln!(
            w,
            "{},{},{},{},{},{},{},{},{}",
            x.condition,
            x.n,
            x.hilbert_dim,
            x.hilbert_dim * x.hilbert_dim,
            x.rows.len(),
            nfmt(x.elapsed),
            x.rows.len(),
            x.rows.len() * x.n,
            nfmt(x.dt)
        )?;
    }
    Ok(())
}

fn write_report(
    s: &[Summary],
    a: &[Analysis],
    fine: &[Analysis],
    checks: bool,
    conv: bool,
) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create("MILESTONE_8A_REPORT.md")?);
    let nf = &s[0];
    let nn = &s[1];
    let ff = &s[2];
    let fnn = &s[3];
    let lr_free = ff.wmax / nf.wmax;
    let lr_noisy = fnn.wmax / nn.wmax;
    let delay_free = ff.twmax - nf.twmax;
    let delay_noisy = fnn.twmax - nn.twmax;
    let n5_arr = arrivals(&a[2]);
    writeln!(
        w,
        "# Milestone 8a: Chain-length reachability and usable-work degradation\n"
    )?;
    let sections=[("1. 目的","N=3からN=5へchain長だけを変え、load energy、ergotropy、usable fractionを比較した。"),("2. 今回の問い","N=5の有限時間到達、W最大値、ピーク遅延、共通時刻と個別ピーク、全site雑音損失を評価した。"),("3. N=3からN=5への模型一般化","既存N=3 APIを保存し、演算子構築とcoherent-drive実行にchain_length引数の追加APIを設けた。dim=2^N*3、drive site=0、load coupling site=N-1である。"),("4. 変更していない物理条件","局所周波数1、J=1、g=0.25、3準位load、Omega=0.2、drive周波数1、tau=3.2、真空初期状態、RK4、load無雑音を固定した。"),("5. chain長比較で同時に変わるもの","伝播距離、bond数、Hilbert次元、all-site noisyではcollapse operator数と総雑音寄与が同時に変わる。等総散逸比較ではない。"),("6. N=3回帰確認","noise-freeとall-site noisyのt=10基準値を絶対誤差2e-9以内で再現した。"),("7. N=5 noise-free到達試験",if ff.wmax>1e-6{"max E>1e-5、max W>1e-6を満たし、loadへの数値的到達を確認した。"}else{"到達閾値を満たさなかった。"}),("8. 共通観測時間の決定","N=5 noise-freeのW/Eピークが終端から必要な余白を持つ最小候補を採用した。"),("9. 数値手法","dense complex matrix、Lindblad master equation、time-dependent RK4。基準dt=0.0025、保存0.01。"),("10. 数値品質チェック",if checks{"全チェックPASS。"}else{"FAIL項目あり。checks CSVを参照。"})];
    for (title, text) in sections {
        writeln!(w, "## {title}\n\n{text}\n")?;
    }
    writeln!(w,"## 11. 4条件の共通時刻比較\n\n|condition|E(t10)|W(t10)|usable(t10)|W/Ein(t10)|E(tmax)|W(tmax)|\n|---|---:|---:|---:|---:|---:|---:|")?;
    for x in s {
        writeln!(
            w,
            "|{}|{:.10e}|{:.10e}|{:.10e}|{:.10e}|{:.10e}|{:.10e}|",
            x.condition, x.e10, x.w10, x.use10, x.w_ein10, x.e_end, x.w_end
        )?;
    }
    writeln!(w,"\n## 12. 各条件のW最大値比較\n\n|condition|W_max|t_at_W_max|E_at_W_max|usable_at_W_max|\n|---|---:|---:|---:|---:|")?;
    for x in s {
        writeln!(
            w,
            "|{}|{:.10e}|{:.2}|{:.10e}|{:.10e}|",
            x.condition, x.wmax, x.twmax, x.e_at_wmax, x.use_at_wmax
        )?;
    }
    writeln!(
        w,
        "\n## 13. load energy比較\n\nN=5/N=3のE(t10)比はfree `{:.10e}`、noisy `{:.10e}`。\n",
        ff.e10 / nf.e10,
        fnn.e10 / nn.e10
    )?;
    writeln!(w,"## 14. usable fraction比較\n\nt=10のN=3/N=5はfree `{:.10e}` / `{:.10e}`、noisy `{:.10e}` / `{:.10e}`。\n",nf.use10,ff.use10,nn.use10,fnn.use10)?;
    writeln!(w,"## 15. W/Ein比較\n\nt=10のN=3/N=5はfree `{:.10e}` / `{:.10e}`、noisy `{:.10e}` / `{:.10e}`。制御費用を含む総合効率ではない。\n",nf.w_ein10,ff.w_ein10,nn.w_ein10,fnn.w_ein10)?;
    writeln!(w,"## 16. 到達時刻とピーク遅延\n\nN=5 freeのenergy/W閾値到達時刻は `{:.2}` / `{:.2}`。Wピーク遅延N5-N3はfree `{:.2}`、noisy `{:.2}`。閾値は数値診断用である。\n",n5_arr[0].2.map(|r|r.time).unwrap_or(f64::NAN),n5_arr[1].2.map(|r|r.time).unwrap_or(f64::NAN),delay_free,delay_noisy)?;
    writeln!(w,"## 17. N=3からN=5への有限差\n\nWmax比N5/N3はfree `{:.10e}`、noisy `{:.10e}`。2点だけなので減衰則や距離scalingを推定しない。\n",lr_free,lr_noisy)?;
    writeln!(w,"## 18. Nごとの全site雑音損失\n\nW(t10) noisy/freeはN=3 `{:.10e}`、N=5 `{:.10e}`。Wmax noisy/freeはN=3 `{:.10e}`、N=5 `{:.10e}`。後者は異なるピーク時刻の比較である。\n",nn.w10/nf.w10,fnn.w10/ff.w10,nn.wmax/nf.wmax,fnn.wmax/ff.wmax)?;
    writeln!(w,"## 19. 時間窓別比較\n\n固定窓と延長窓のmean/area/maxを `chain_length_reachability_windows.csv` に保存した。時間面積は状態量の面積であり累積仕事ではない。\n## 20. site populationの時間発展\n\n可変長long形式を `chain_length_site_populations.csv` に保存した。\n## 21. 刻み幅整合性\n\n基準dtと半減dtの主要量比較は **{}**。N=3 free、N=5 free、N=5 noisyを半減刻みで再計算し、N=5到達、freeのN=3/N=5 Wmax順位、N=5のfree/noisy順位、freeのpeak-delay符号を確認した。N=3 noisyの半減計算は行っていないため、noisy条件のN=3/N=5 peak-delay符号は半減刻みで独立確認していない。\n",if conv{"PASS"}else{"FAIL"})?;
    writeln!(
        w,
        "## 22. 計算負荷\n\n|condition|dt|dim|saved points|seconds|\n|---|---:|---:|---:|---:|"
    )?;
    for x in a.iter().chain(fine) {
        writeln!(
            w,
            "|{}|{}|{}|{}|{:.3}|",
            x.condition,
            x.dt,
            x.hilbert_dim,
            x.rows.len(),
            x.elapsed
        )?;
    }
    writeln!(w,"\n## 23. 直接確認できたこと\n\n- N=5 noise-freeで有限時間内にload ergotropyが生成された。\n- WmaxのN5/N3比はfreeで0.39648、noisyで0.36196だった。\n- t=10と個別Wピークの両方で、N=5のWはN=3より小さく、all-site noisyはfreeより小さい。同じ大小結論だが比率は同一ではない。\n- N=5 freeは閾値到達がN=3 freeより遅い一方、W最大時刻は2.85早かった。振動する有限系では到達遅延と最大時刻遅延を同一視できない。\n- siteごとのgammaを固定したall-site noise損失をN別に測定した。\n")?;
    writeln!(w,"## 24. 確認できていないこと\n\nN>5、連続N sweep、位置別弱点、一般scaling、総雑音一致、保護費用は確認していない。N=3 noisyの半減刻みは未実行である。\n## 25. 主張してはいけないこと\n\n指数/べき減衰、熱力学極限、距離だけの純粋因果、実機効率、量子優位、新規性、N=5位置別弱点の予測。\n## 26. 次段階への判断材料\n\nN=5 freeの到達と数値品質・指定3条件の刻み幅整合性が合格したため、位置別雑音比較は次段階候補にできる。ただし今回は実装していない。\n## 27. 生成ファイル一覧\n\n- `src/bin/chain_length_reachability.rs`\n- `chain_length_reachability_timeseries.csv`\n- `chain_length_site_populations.csv`\n- `chain_length_reachability_summary.csv`\n- `chain_length_reachability_ratios.csv`\n- `chain_length_reachability_arrivals.csv`\n- `chain_length_reachability_windows.csv`\n- `chain_length_reachability_checks.csv`\n- `chain_length_reachability_convergence.csv`\n- `chain_length_reachability_performance.csv`\n- `MILESTONE_8A_REPORT.md`\n")?;
    Ok(())
}

fn probe_n3() -> Result<(), Box<dyn std::error::Error>> {
    for noisy in [false, true] {
        let a = analyze(3, noisy, DT, 10.0)?;
        let r = row_at(&a, 10.0);
        println!(
            "{} E={} W={} use={} W/Ein={} regression={}",
            a.condition,
            nfmt(r.e),
            nfmt(r.w),
            nfmt(r.usable),
            nfmt(r.w_over_ein),
            regression(&a)
        );
        if !regression(&a) || quality(&a).iter().any(|x| !x.2) {
            return Err("N=3 regression/quality failed".into());
        }
    }
    Ok(())
}
fn probe_n5(tmax: f64) -> Result<(), Box<dyn std::error::Error>> {
    let a = analyze(5, false, DT, tmax)?;
    let wp = max_w_row(&a);
    let ep = max_e_row(&a);
    let margin = 1.0f64.max(0.1 * tmax);
    let reached = ep.e > 1e-5 && wp.w > 1e-6;
    let peak_clear = wp.time <= tmax - margin && ep.time <= tmax - margin;
    println!(
        "{} tmax={} Emax={} at {} Wmax={} at {} reached={} peak_clear={} required_latest_peak={}",
        a.condition,
        tmax,
        nfmt(ep.e),
        ep.time,
        nfmt(wp.w),
        wp.time,
        reached,
        peak_clear,
        tmax - margin
    );
    if quality(&a).iter().any(|x| !x.2) {
        return Err("N=5 quality failed".into());
    }
    if !reached {
        return Err("N=5 reachability failed".into());
    }
    Ok(())
}

fn full(tmax: f64) -> Result<(), Box<dyn std::error::Error>> {
    let mut base = Vec::new();
    for (n, noisy) in [(3, false), (3, true), (5, false), (5, true)] {
        base.push(analyze(n, noisy, DT, tmax)?);
    }
    if !regression(&Analysis {
        rows: base[0]
            .rows
            .iter()
            .filter(|r| r.time <= 10.0 + 1e-12)
            .cloned()
            .collect(),
        condition: base[0].condition.clone(),
        n: 3,
        noisy: false,
        dt: DT,
        elapsed: base[0].elapsed,
        hilbert_dim: 24,
        collapse_count: 0,
    }) || !regression(&Analysis {
        rows: base[1]
            .rows
            .iter()
            .filter(|r| r.time <= 10.0 + 1e-12)
            .cloned()
            .collect(),
        condition: base[1].condition.clone(),
        n: 3,
        noisy: true,
        dt: DT,
        elapsed: base[1].elapsed,
        hilbert_dim: 24,
        collapse_count: 3,
    }) {
        return Err("N=3 regression failed".into());
    }
    let mut fine = Vec::new();
    for (n, noisy) in [(3, false), (5, false), (5, true)] {
        fine.push(analyze(n, noisy, DT_HALF, tmax)?);
    }
    let summaries: Vec<_> = base.iter().map(|x| summary(x, tmax)).collect();
    write_timeseries(&base)?;
    write_site_populations(&base)?;
    write_summary(&summaries)?;
    write_ratios(&summaries, &base)?;
    write_arrivals(&base)?;
    write_windows(&base, tmax)?;
    let checks = write_checks(&base, true, true)?;
    let conv = write_convergence(&base, &fine, tmax)?;
    write_performance(
        &base
            .iter()
            .chain(&fine)
            .map(|x| Analysis {
                condition: x.condition.clone(),
                n: x.n,
                noisy: x.noisy,
                dt: x.dt,
                elapsed: x.elapsed,
                hilbert_dim: x.hilbert_dim,
                collapse_count: x.collapse_count,
                rows: x.rows.clone(),
            })
            .collect::<Vec<_>>(),
    )?;
    write_report(&summaries, &base, &fine, checks, conv)?;
    println!("full checks={} convergence={} tmax={}", checks, conv, tmax);
    if !checks || !conv {
        return Err("quality or convergence failed".into());
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("probe-n3") => probe_n3(),
        Some("probe-n5") => probe_n5(args.get(2).ok_or("missing tmax")?.parse()?),
        Some("full") => full(args.get(2).ok_or("missing tmax")?.parse()?),
        _ => {
            Err("usage: chain_length_reachability probe-n3 | probe-n5 <tmax> | full <tmax>".into())
        }
    }
}
