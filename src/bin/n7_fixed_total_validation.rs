mod validation {
    include!("fixed_total_noise_comparison.rs");

    use nalgebra::linalg::Schur;

    const V_OUTPUTS: [&str; 9] = [
        "robust_eigen_diagnostic_unit_checks.csv",
        "n7_fixed_total_validation_timeseries.csv",
        "n7_fixed_total_validation_eigen_diagnostics.csv",
        "n7_fixed_total_validation_summary.csv",
        "n7_fixed_total_validation_trajectory_comparison.csv",
        "n7_fixed_total_validation_checks.csv",
        "n7_fixed_total_validation_performance.csv",
        "fixed_total_noise_final_comparison.csv",
        "MILESTONE_9C_VALIDATION.md",
    ];
    const SUM_TRACE_TOL: f64 = 1.0e-10;
    const EIG_IMAG_TOL: f64 = 1.0e-10;
    const INPUT_HERM_TOL: f64 = 1.0e-10;
    const TRAJECTORY_TOL: f64 = 1.0e-12;
    const N: usize = 7;
    const LOAD_DIM: usize = 3;
    const DIM: usize = 384;
    const GAMMA_SITE: f64 = TOTAL_GAMMA / N as f64;

    #[derive(Clone, Debug)]
    struct SolverAttempt {
        attempted: bool,
        all_finite: bool,
        minimum: f64,
        max_imag: f64,
        sum_trace_difference: f64,
        pass: bool,
    }

    impl SolverAttempt {
        fn not_attempted() -> Self {
            Self {
                attempted: false,
                all_finite: false,
                minimum: f64::NAN,
                max_imag: f64::NAN,
                sum_trace_difference: f64::NAN,
                pass: false,
            }
        }
    }

    #[derive(Clone, Debug)]
    struct PositivityDiagnostic {
        time: f64,
        input_finite: bool,
        trace_error: f64,
        hermiticity_error: f64,
        correction_norm: f64,
        primary: SolverAttempt,
        fallback: SolverAttempt,
        selected_solver: &'static str,
        selected_minimum: f64,
        selected_sum_trace_difference: f64,
        fallback_used: bool,
        positivity_pass: bool,
        solver_failure: bool,
    }

    #[derive(Clone)]
    struct ValidationRun {
        rows: Vec<Row>,
        eigen: Vec<PositivityDiagnostic>,
        construction_seconds: f64,
        kernel_construction_seconds: f64,
        propagation_seconds: f64,
        diagnostics_seconds: f64,
        total_seconds: f64,
        working_set_before: u64,
        working_set_after: u64,
        peak_working_set: u64,
    }

    #[derive(Clone)]
    struct ComparisonRow {
        metric: String,
        compared_rows: usize,
        max_difference: f64,
        differing_rows: usize,
        pass: bool,
    }

    #[derive(Clone)]
    struct UnitCheck {
        name: &'static str,
        observed: String,
        expected: &'static str,
        pass: bool,
    }

    fn vfmt(x: f64) -> String {
        if x.is_nan() {
            "NaN".into()
        } else if x == f64::INFINITY {
            "+Inf".into()
        } else if x == f64::NEG_INFINITY {
            "-Inf".into()
        } else {
            format!("{x:.16e}")
        }
    }

    fn all_elements_finite(rho: &ComplexMatrix) -> bool {
        rho.iter().all(|z| z.re.is_finite() && z.im.is_finite())
    }

    fn eigenvalue_sum_trace_pass(difference: f64, tolerance: f64) -> bool {
        difference.is_finite() && difference <= tolerance
    }

    fn nonfinite_minimum(values: impl Iterator<Item = f64>) -> f64 {
        let values: Vec<f64> = values.collect();
        if values.iter().any(|x| x.is_nan()) {
            f64::NAN
        } else if values.iter().any(|x| *x == f64::NEG_INFINITY) {
            f64::NEG_INFINITY
        } else if values.iter().any(|x| *x == f64::INFINITY) {
            f64::INFINITY
        } else {
            values.into_iter().fold(f64::INFINITY, f64::min)
        }
    }

    fn try_symmetric_eigen(rho_h: &ComplexMatrix) -> SolverAttempt {
        let values = SymmetricEigen::new(rho_h.clone()).eigenvalues;
        let all_finite = values.iter().all(|x| x.is_finite());
        let minimum = if all_finite {
            values.iter().copied().fold(f64::INFINITY, f64::min)
        } else {
            nonfinite_minimum(values.iter().copied())
        };
        let sum: f64 = values.iter().sum();
        let difference = (C64::new(sum, 0.0) - rho_h.trace()).norm();
        let pass = all_finite
            && minimum.is_finite()
            && sum.is_finite()
            && eigenvalue_sum_trace_pass(difference, SUM_TRACE_TOL);
        SolverAttempt {
            attempted: true,
            all_finite,
            minimum,
            max_imag: 0.0,
            sum_trace_difference: difference,
            pass,
        }
    }

    fn try_complex_schur(rho_h: &ComplexMatrix) -> SolverAttempt {
        let (_, triangular) = Schur::new(rho_h.clone()).unpack();
        let values: Vec<C64> = (0..triangular.nrows())
            .map(|i| triangular[(i, i)])
            .collect();
        let all_finite = values.iter().all(|z| z.re.is_finite() && z.im.is_finite());
        let minimum = if all_finite {
            values.iter().map(|z| z.re).fold(f64::INFINITY, f64::min)
        } else {
            nonfinite_minimum(values.iter().map(|z| z.re))
        };
        let max_imag = values.iter().map(|z| z.im.abs()).fold(0.0, f64::max);
        let sum: C64 = values.iter().copied().sum();
        let difference = (sum - rho_h.trace()).norm();
        let pass = all_finite
            && minimum.is_finite()
            && max_imag.is_finite()
            && max_imag <= EIG_IMAG_TOL
            && eigenvalue_sum_trace_pass(difference, SUM_TRACE_TOL);
        SolverAttempt {
            attempted: true,
            all_finite,
            minimum,
            max_imag,
            sum_trace_difference: difference,
            pass,
        }
    }

    fn select_attempts(
        primary: SolverAttempt,
        fallback: SolverAttempt,
    ) -> (&'static str, f64, f64, bool, bool) {
        if primary.pass {
            (
                "symmetric_eigen",
                primary.minimum,
                primary.sum_trace_difference,
                false,
                false,
            )
        } else if fallback.pass {
            (
                "complex_schur_fallback",
                fallback.minimum,
                fallback.sum_trace_difference,
                true,
                false,
            )
        } else {
            ("none", f64::NAN, f64::NAN, fallback.attempted, true)
        }
    }

    fn evaluate_positivity(rho: &ComplexMatrix, time: f64) -> PositivityDiagnostic {
        let input_finite =
            all_elements_finite(rho) && rho.trace().re.is_finite() && rho.trace().im.is_finite();
        let trace_error = (rho.trace() - C64::new(1.0, 0.0)).norm();
        let herm_error = hermiticity_error(rho);
        let rho_h = (rho + rho.adjoint()) * C64::new(0.5, 0.0);
        let correction = frobenius_norm(&(&rho_h - rho));
        if !input_finite || herm_error > INPUT_HERM_TOL {
            return PositivityDiagnostic {
                time,
                input_finite,
                trace_error,
                hermiticity_error: herm_error,
                correction_norm: correction,
                primary: SolverAttempt::not_attempted(),
                fallback: SolverAttempt::not_attempted(),
                selected_solver: "state_input_invalid",
                selected_minimum: f64::NAN,
                selected_sum_trace_difference: f64::NAN,
                fallback_used: false,
                positivity_pass: false,
                solver_failure: true,
            };
        }
        let primary = try_symmetric_eigen(&rho_h);
        let fallback = if primary.pass {
            SolverAttempt::not_attempted()
        } else {
            try_complex_schur(&rho_h)
        };
        let (selected, minimum, difference, fallback_used, solver_failure) =
            select_attempts(primary.clone(), fallback.clone());
        PositivityDiagnostic {
            time,
            input_finite,
            trace_error,
            hermiticity_error: herm_error,
            correction_norm: correction,
            primary,
            fallback,
            selected_solver: selected,
            selected_minimum: minimum,
            selected_sum_trace_difference: difference,
            fallback_used,
            positivity_pass: !solver_failure && minimum >= -POS_TOL,
            solver_failure,
        }
    }

    fn diagnose_validation(
        spec: Spec,
        rho: &ComplexMatrix,
        time: f64,
        ops: &Operators,
        params: &ModelParams,
        drive_in: f64,
        drive_net: f64,
        dephasing_net: f64,
        bare0: f64,
        drive_power: C64,
        dephasing_power: C64,
        positivity: &PositivityDiagnostic,
    ) -> Result<Row, Box<dyn std::error::Error>> {
        let load = partial_trace(rho, &ops.dims, &[spec.n])?;
        let h_load = ComplexMatrix::from_diagonal(&nalgebra::DVector::from_iterator(
            params.load_dim,
            (0..params.load_dim).map(|level| C64::new(level as f64 * params.omega_load, 0.0)),
        ));
        let work = ergotropy(&load, &h_load, 1.0e-9)?;
        let mut diagonal = ComplexMatrix::zeros(params.load_dim, params.load_dim);
        for level in 0..params.load_dim {
            diagonal[(level, level)] = load[(level, level)];
        }
        let diagonal_work = ergotropy(&diagonal, &h_load, 1.0e-9)?.ergotropy;
        let coherence_l1: f64 = (0..params.load_dim)
            .flat_map(|i| (0..params.load_dim).map(move |j| (i, j)))
            .filter(|(i, j)| i != j)
            .map(|(i, j)| load[(i, j)].norm())
            .sum();
        let sites: Vec<f64> = ops
            .number_sites
            .iter()
            .map(|number| expectation(rho, number).re)
            .collect();
        let chain_population = sites.iter().sum();
        let bare_energy = expectation(rho, &ops.h_total).re;
        let load_populations = [load[(0, 0)].re, load[(1, 1)].re, load[(2, 2)].re];
        let values = [
            work.energy,
            work.ergotropy,
            diagonal_work,
            coherence_l1,
            drive_in,
            drive_net,
            dephasing_net,
            chain_population,
            bare_energy,
            drive_power.re,
            drive_power.im,
            dephasing_power.re,
            dephasing_power.im,
            positivity.trace_error,
            positivity.hermiticity_error,
        ];
        Ok(Row {
            time,
            envelope: drive_envelope(time, config(spec.gamma_site()).tau),
            energy: work.energy,
            work: work.ergotropy,
            diagonal_work,
            coherence_work: work.ergotropy - diagonal_work,
            coherence_l1,
            usable: ratio(work.ergotropy, work.energy),
            drive_in,
            drive_net,
            dephasing_net,
            w_over_ein: ratio(work.ergotropy, drive_in),
            sites,
            chain_population,
            load_populations,
            bare_energy,
            drive_power: drive_power.re,
            drive_power_imag: drive_power.im,
            dephasing_power: dephasing_power.re,
            dephasing_power_imag: dephasing_power.im,
            trace_error: positivity.trace_error,
            herm_error: positivity.hermiticity_error,
            min_eigenvalue: positivity.selected_minimum,
            ledger: bare_energy - bare0 - drive_net - dephasing_net,
            finite: values.iter().all(|x| x.is_finite()) && all_elements_finite(rho),
            reduced_trace_error: (load.trace() - C64::new(1.0, 0.0)).norm(),
        })
    }

    fn run_validation(
    ) -> Result<(ValidationRun, Vec<(String, String, String, bool)>), Box<dyn std::error::Error>>
    {
        let total_start = Instant::now();
        let construction_start = Instant::now();
        let spec = N7_SPEC;
        let params = ModelParams::default();
        let ops = build_operators_for_chain(&params, spec.n)?;
        let gammas = spec.gammas();
        let kernel_start = Instant::now();
        let kernel = DiagonalDephasingKernel::new(spec.n, params.load_dim, &gammas)?;
        let kernel_seconds = kernel_start.elapsed().as_secs_f64();
        let mut rho = ComplexMatrix::zeros(spec.dim, spec.dim);
        rho[(0, 0)] = C64::new(1.0, 0.0);
        let construction = construction_checks(spec, &gammas, &ops, &rho, &kernel)?;
        if let Some(failed) = construction.iter().find(|x| !x.3) {
            return Err(format!("construction failed: {} observed {}", failed.0, failed.1).into());
        }
        let construction_seconds = construction_start.elapsed().as_secs_f64();
        let (working_set_before, _) = process_memory();
        let bare0 = expectation(&rho, &ops.h_total).re;
        let mut drive_in = 0.0;
        let mut drive_net = 0.0;
        let mut dephasing_net = 0.0;
        let (power0, dephasing0) =
            instantaneous_powers(&rho, 0.0, &ops, &kernel, spec.gamma_site())?;
        let mut previous_power = power0.re;
        let mut previous_dephasing = dephasing0.re;
        let mut rows = Vec::with_capacity(1001);
        let mut eigen = Vec::with_capacity(1001);
        let diag_start = Instant::now();
        let positivity0 = evaluate_positivity(&rho, 0.0);
        rows.push(diagnose_validation(
            spec,
            &rho,
            0.0,
            &ops,
            &params,
            drive_in,
            drive_net,
            dephasing_net,
            bare0,
            power0,
            dephasing0,
            &positivity0,
        )?);
        eigen.push(positivity0);
        let mut diagnostics_seconds = diag_start.elapsed().as_secs_f64();
        let mut propagation_seconds = 0.0;
        for step in 0..4000 {
            let time = step as f64 * DT;
            let start = Instant::now();
            rho = rk4_step(&rho, time, &ops, &kernel, spec.gamma_site())?;
            propagation_seconds += start.elapsed().as_secs_f64();
            if (step + 1) % SAVE_STEPS == 0 {
                let now = (step + 1) as f64 * DT;
                let start = Instant::now();
                let (power, dephasing) =
                    instantaneous_powers(&rho, now, &ops, &kernel, spec.gamma_site())?;
                drive_net += 0.5 * SAVE_INTERVAL * (previous_power + power.re);
                drive_in += 0.5 * SAVE_INTERVAL * (previous_power.max(0.0) + power.re.max(0.0));
                dephasing_net += 0.5 * SAVE_INTERVAL * (previous_dephasing + dephasing.re);
                previous_power = power.re;
                previous_dephasing = dephasing.re;
                let positivity = evaluate_positivity(&rho, now);
                if !positivity.input_finite
                    || positivity.trace_error > TRACE_TOL
                    || positivity.hermiticity_error > HERM_TOL
                    || positivity.solver_failure
                    || !positivity.positivity_pass
                {
                    return Err(format!(
                        "quality stop at t={now}: solver={} min={}",
                        positivity.selected_solver,
                        vfmt(positivity.selected_minimum)
                    )
                    .into());
                }
                rows.push(diagnose_validation(
                    spec,
                    &rho,
                    now,
                    &ops,
                    &params,
                    drive_in,
                    drive_net,
                    dephasing_net,
                    bare0,
                    power,
                    dephasing,
                    &positivity,
                )?);
                eigen.push(positivity);
                diagnostics_seconds += start.elapsed().as_secs_f64();
                if rows.len() % 100 == 1 {
                    println!("N7_fixed_total_validation progress t={now:.2} saved={} propagation={:.1}s diagnostics={:.1}s", rows.len(), propagation_seconds, diagnostics_seconds);
                }
            }
        }
        let total_seconds = total_start.elapsed().as_secs_f64();
        let (working_set_after, peak_working_set) = process_memory();
        Ok((
            ValidationRun {
                rows,
                eigen,
                construction_seconds,
                kernel_construction_seconds: kernel_seconds,
                propagation_seconds,
                diagnostics_seconds,
                total_seconds,
                working_set_before,
                working_set_after,
                peak_working_set,
            },
            construction,
        ))
    }

    fn write_validation_timeseries(rows: &[Row]) -> Result<(), Box<dyn std::error::Error>> {
        let mut w = BufWriter::new(File::create(V_OUTPUTS[1])?);
        writeln!(w, "condition,chain_length,noise_normalization,noisy_site_count,gamma_phi_per_site,gamma_phi_total,hilbert_dimension,time,Omega,drive_envelope,load_energy,load_ergotropy,load_diagonal_ergotropy,load_coherence_ergotropy,load_coherence_l1,usable_fraction,drive_energy_in,drive_energy_net,dephasing_energy_net,W_over_Ein,total_chain_population,load_population_0,load_population_1,load_population_2,load_top_level_population,bare_network_energy,drive_power,dephasing_power,trace_error,hermiticity_error,min_eigenvalue,energy_ledger_residual")?;
        for r in rows {
            writeln!(w, "N7_fixed_total_validation,7,fixed_total_gamma_1p5,7,{},{},384,{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}", vfmt(GAMMA_SITE), vfmt(TOTAL_GAMMA), vfmt(r.time), vfmt(OMEGA), vfmt(r.envelope), vfmt(r.energy), vfmt(r.work), vfmt(r.diagonal_work), vfmt(r.coherence_work), vfmt(r.coherence_l1), vfmt(r.usable), vfmt(r.drive_in), vfmt(r.drive_net), vfmt(r.dephasing_net), vfmt(r.w_over_ein), vfmt(r.chain_population), vfmt(r.load_populations[0]), vfmt(r.load_populations[1]), vfmt(r.load_populations[2]), vfmt(r.load_populations[2]), vfmt(r.bare_energy), vfmt(r.drive_power), vfmt(r.dephasing_power), vfmt(r.trace_error), vfmt(r.herm_error), vfmt(r.min_eigenvalue), vfmt(r.ledger))?;
        }
        Ok(())
    }

    fn write_eigen_diagnostics(
        rows: &[PositivityDiagnostic],
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut w = BufWriter::new(File::create(V_OUTPUTS[2])?);
        writeln!(w, "condition,time,input_finite,trace_error,hermiticity_error,hermitianization_correction_norm,primary_solver,primary_all_eigenvalues_finite,primary_minimum_eigenvalue,primary_eigenvalue_sum_trace_difference,fallback_attempted,fallback_solver,fallback_all_eigenvalues_finite,fallback_minimum_eigenvalue,fallback_max_eigenvalue_imaginary_part,fallback_eigenvalue_sum_trace_difference,selected_solver,selected_minimum_eigenvalue,fallback_used,positivity_pass,solver_failure")?;
        for x in rows {
            writeln!(w, "N7_fixed_total_validation,{},{},{},{},{},symmetric_eigen,{},{},{},{},complex_schur,{},{},{},{},{},{},{},{},{}", vfmt(x.time), x.input_finite, vfmt(x.trace_error), vfmt(x.hermiticity_error), vfmt(x.correction_norm), x.primary.all_finite, vfmt(x.primary.minimum), vfmt(x.primary.sum_trace_difference), x.fallback.attempted, x.fallback.all_finite, vfmt(x.fallback.minimum), vfmt(x.fallback.max_imag), vfmt(x.fallback.sum_trace_difference), x.selected_solver, vfmt(x.selected_minimum), x.fallback_used, x.positivity_pass, x.solver_failure)?;
        }
        Ok(())
    }

    fn parse_csv(
        path: &str,
    ) -> Result<(Vec<String>, Vec<Vec<String>>), Box<dyn std::error::Error>> {
        let mut lines = BufReader::new(File::open(path)?).lines();
        let header = lines
            .next()
            .ok_or("empty CSV")??
            .split(',')
            .map(str::to_string)
            .collect();
        let rows = lines
            .map(|line| Ok(line?.split(',').map(str::to_string).collect()))
            .collect::<Result<Vec<_>, Box<dyn std::error::Error>>>()?;
        Ok((header, rows))
    }

    fn trajectory_comparison(
        rows: &[Row],
    ) -> Result<Vec<ComparisonRow>, Box<dyn std::error::Error>> {
        let (header, all) = parse_csv("fixed_total_noise_timeseries.csv")?;
        let condition = header
            .iter()
            .position(|x| x == "condition")
            .ok_or("missing condition")?;
        let old: Vec<&Vec<String>> = all
            .iter()
            .filter(|x| x[condition] == "N7_fixed_total_noise")
            .collect();
        if old.len() != rows.len() {
            return Err(format!("old N7 rows={} new={}", old.len(), rows.len()).into());
        }
        let specs: Vec<(&str, Box<dyn Fn(&Row) -> f64>)> = vec![
            ("load_energy", Box::new(|r| r.energy)),
            ("load_ergotropy", Box::new(|r| r.work)),
            ("load_coherence_l1", Box::new(|r| r.coherence_l1)),
            ("usable_fraction", Box::new(|r| r.usable)),
            ("drive_energy_in", Box::new(|r| r.drive_in)),
            ("drive_energy_net", Box::new(|r| r.drive_net)),
            ("dephasing_energy_net", Box::new(|r| r.dephasing_net)),
            ("total_chain_population", Box::new(|r| r.chain_population)),
            ("load_population_0", Box::new(|r| r.load_populations[0])),
            ("load_population_1", Box::new(|r| r.load_populations[1])),
            ("load_population_2", Box::new(|r| r.load_populations[2])),
            ("bare_network_energy", Box::new(|r| r.bare_energy)),
            ("drive_power", Box::new(|r| r.drive_power)),
            ("dephasing_power", Box::new(|r| r.dephasing_power)),
            ("trace_error", Box::new(|r| r.trace_error)),
            ("hermiticity_error", Box::new(|r| r.herm_error)),
            ("energy_ledger_residual", Box::new(|r| r.ledger)),
        ];
        let mut output = Vec::new();
        for (metric, value) in specs {
            let column = header
                .iter()
                .position(|x| x == metric)
                .ok_or_else(|| format!("missing {metric}"))?;
            let mut max_difference = 0.0_f64;
            let mut differing = 0;
            for (new, old) in rows.iter().zip(&old) {
                let old_value: f64 = old[column].parse()?;
                let new_value = value(new);
                let difference = if old_value.is_nan() && new_value.is_nan() {
                    0.0
                } else {
                    (old_value - new_value).abs()
                };
                max_difference = max_difference.max(difference);
                if !difference.is_finite() || difference > TRAJECTORY_TOL {
                    differing += 1;
                }
            }
            output.push(ComparisonRow {
                metric: metric.into(),
                compared_rows: rows.len(),
                max_difference,
                differing_rows: differing,
                pass: differing == 0,
            });
        }
        Ok(output)
    }

    fn write_trajectory_comparison(
        rows: &[ComparisonRow],
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut w = BufWriter::new(File::create(V_OUTPUTS[4])?);
        writeln!(
            w,
            "metric,compared_rows,max_absolute_difference,tolerance,differing_rows,pass"
        )?;
        for x in rows {
            writeln!(
                w,
                "{},{},{},{},{},{}",
                x.metric,
                x.compared_rows,
                vfmt(x.max_difference),
                vfmt(TRAJECTORY_TOL),
                x.differing_rows,
                x.pass
            )?;
        }
        Ok(())
    }

    fn summary_counts(eigen: &[PositivityDiagnostic]) -> (usize, usize, usize, usize, usize) {
        let primary_success = eigen.iter().filter(|x| x.primary.pass).count();
        let primary_failure = eigen.len() - primary_success;
        let fallback_attempt = eigen.iter().filter(|x| x.fallback.attempted).count();
        let fallback_success = eigen
            .iter()
            .filter(|x| x.fallback.attempted && x.fallback.pass)
            .count();
        let solver_failure = eigen.iter().filter(|x| x.solver_failure).count();
        (
            primary_success,
            primary_failure,
            fallback_attempt,
            fallback_success,
            solver_failure,
        )
    }

    fn write_validation_summary(
        run: &ValidationRun,
        summary: &Summary,
        final_status: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let (ps, pf, fa, fs, sf) = summary_counts(&run.eigen);
        let worst = run
            .eigen
            .iter()
            .map(|x| x.selected_minimum)
            .fold(f64::INFINITY, f64::min);
        let max_imag = run
            .eigen
            .iter()
            .filter(|x| x.fallback.attempted)
            .map(|x| x.fallback.max_imag)
            .fold(0.0, f64::max);
        let max_diff = run
            .eigen
            .iter()
            .map(|x| x.selected_sum_trace_difference)
            .fold(0.0, f64::max);
        let mut w = BufWriter::new(File::create(V_OUTPUTS[3])?);
        writeln!(w,"condition,chain_length,gamma_phi_per_site,gamma_phi_total,saved_time_points,primary_success_count,primary_failure_count,fallback_attempt_count,fallback_success_count,solver_failure_count,worst_selected_minimum_eigenvalue,max_fallback_eigenvalue_imaginary_part,max_selected_sum_trace_difference,E_at_t10,W_at_t10,usable_fraction_at_t10,W_over_Ein_at_t10,W_max,t_at_W_max,energy_arrival_time,ergotropy_arrival_time,E_time_area_0_to_t10,W_time_area_0_to_t10,endpoint_peak_classification,final_status")?;
        writeln!(w,"N7_fixed_total_validation,7,{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",vfmt(GAMMA_SITE),vfmt(TOTAL_GAMMA),run.rows.len(),ps,pf,fa,fs,sf,vfmt(worst),vfmt(max_imag),vfmt(max_diff),vfmt(summary.endpoint.energy),vfmt(summary.endpoint.work),vfmt(summary.endpoint.usable),vfmt(summary.endpoint.w_over_ein),vfmt(summary.w_max.work),vfmt(summary.w_max.time),vfmt(summary.arrivals[0].time),vfmt(summary.arrivals[1].time),vfmt(summary.e_area),vfmt(summary.w_area),summary.peak_class,final_status)?;
        Ok(())
    }

    fn write_validation_performance(run: &ValidationRun) -> Result<(), Box<dyn std::error::Error>> {
        let mut w = BufWriter::new(File::create(V_OUTPUTS[6])?);
        writeln!(w,"condition,chain_length,hilbert_dimension,rk4_steps,saved_time_points,construction_seconds,kernel_construction_seconds,propagation_seconds,diagnostics_seconds,total_seconds,working_set_before_bytes,working_set_after_bytes,peak_working_set_bytes")?;
        writeln!(
            w,
            "N7_fixed_total_validation,7,384,4000,1001,{},{},{},{},{},{},{},{}",
            vfmt(run.construction_seconds),
            vfmt(run.kernel_construction_seconds),
            vfmt(run.propagation_seconds),
            vfmt(run.diagnostics_seconds),
            vfmt(run.total_seconds),
            run.working_set_before,
            run.working_set_after,
            run.peak_working_set
        )?;
        Ok(())
    }

    fn read_metric_value(
        path: &str,
        metric: &str,
        value_column: &str,
        condition_value: Option<&str>,
    ) -> Result<f64, Box<dyn std::error::Error>> {
        let (header, rows) = parse_csv(path)?;
        let metric_i = header.iter().position(|x| x == "metric");
        let condition_i = header.iter().position(|x| x == "condition");
        let value_i = header
            .iter()
            .position(|x| x == value_column)
            .ok_or("missing value column")?;
        for row in rows {
            let metric_match = metric_i.is_none_or(|i| row[i] == metric);
            let condition_match =
                condition_value.is_none_or(|wanted| condition_i.is_some_and(|i| row[i] == wanted));
            if metric_match && condition_match {
                return Ok(row[value_i].parse()?);
            }
        }
        Err(format!("value not found {path} {metric} {value_column}").into())
    }

    fn write_final_comparison(summary: &Summary) -> Result<(f64, f64), Box<dyn std::error::Error>> {
        let n3 = read_metric_value(
            "fixed_total_noise_length_comparison.csv",
            "W_max",
            "N3_value",
            None,
        )?;
        let n5 = read_metric_value(
            "fixed_total_noise_summary.csv",
            "",
            "W_max",
            Some("N5_fixed_total_noise"),
        )?;
        let n7 = summary.w_max.work;
        let ratio_75 = n7 / n5;
        let ratio_73 = n7 / n3;
        let class = if ratio_75 < 0.9 {
            "N7 substantially lower"
        } else if ratio_75 <= 1.1 {
            "similar within 10 percent descriptive band"
        } else {
            "N7 substantially higher"
        };
        let mut w = BufWriter::new(File::create(V_OUTPUTS[7])?);
        writeln!(w,"condition,chain_length,W_max,source,N7_over_N5,N7_over_N3,N7_minus_N5,N7_minus_N3,N7_vs_N5_classification")?;
        writeln!(
            w,
            "N3_fixed_total_noise,3,{},existing_9c_reference,NaN,NaN,NaN,NaN,not_applicable",
            vfmt(n3)
        )?;
        writeln!(
            w,
            "N5_fixed_total_noise,5,{},existing_9c_reference,NaN,NaN,NaN,NaN,not_applicable",
            vfmt(n5)
        )?;
        writeln!(
            w,
            "N7_fixed_total_validation,7,{},validation_rerun,{},{},{},{},{}",
            vfmt(n7),
            vfmt(ratio_75),
            vfmt(ratio_73),
            vfmt(n7 - n5),
            vfmt(n7 - n3),
            class
        )?;
        Ok((ratio_75, ratio_73))
    }

    fn make_unit_checks() -> Result<Vec<UnitCheck>, Box<dyn std::error::Error>> {
        let matrix = |diag: &[f64]| {
            ComplexMatrix::from_diagonal(&nalgebra::DVector::from_iterator(
                diag.len(),
                diag.iter().map(|x| C64::new(*x, 0.0)),
            ))
        };
        let normal = evaluate_positivity(&matrix(&[0.25, 0.75]), 0.0);
        let pure = evaluate_positivity(&matrix(&[1.0, 0.0]), 0.0);
        let tiny = evaluate_positivity(&matrix(&[1.0 + 1e-12, -1e-12]), 0.0);
        let negative = evaluate_positivity(&matrix(&[1.0 + 1e-6, -1e-6]), 0.0);
        let failed = SolverAttempt {
            attempted: true,
            all_finite: false,
            minimum: f64::NAN,
            max_imag: 0.0,
            sum_trace_difference: f64::NAN,
            pass: false,
        };
        let fallback = SolverAttempt {
            attempted: true,
            all_finite: true,
            minimum: -1e-12,
            max_imag: 0.0,
            sum_trace_difference: 1e-15,
            pass: true,
        };
        let selection = select_attempts(failed.clone(), fallback);
        let both = select_attempts(failed.clone(), failed);
        let finite_matrix = matrix(&[1.0, 0.0]);
        let independent = all_elements_finite(&finite_matrix);
        let params = ModelParams::default();
        let ops = build_operators_for_chain(&params, N)?;
        let kernel = DiagonalDephasingKernel::new(N, LOAD_DIM, &vec![GAMMA_SITE; N])?;
        let mut rho = ComplexMatrix::zeros(DIM, DIM);
        rho[(0, 0)] = C64::new(1.0, 0.0);
        for step in 0..8 {
            rho = rk4_step(&rho, step as f64 * DT, &ops, &kernel, GAMMA_SITE)?;
        }
        let regression = evaluate_positivity(&rho, 0.02);
        Ok(vec![
            UnitCheck {
                name: "normal_2x2_hermitian",
                observed: format!(
                    "solver={} min={}",
                    normal.selected_solver,
                    vfmt(normal.selected_minimum)
                ),
                expected: "primary success minimum 0.25",
                pass: normal.selected_solver == "symmetric_eigen"
                    && (normal.selected_minimum - 0.25).abs() < 1e-15,
            },
            UnitCheck {
                name: "pure_state",
                observed: vfmt(pure.selected_minimum),
                expected: "minimum approximately 0 and positivity PASS",
                pass: pure.positivity_pass && pure.selected_minimum.abs() < 1e-15,
            },
            UnitCheck {
                name: "tiny_negative",
                observed: vfmt(tiny.selected_minimum),
                expected: "-1e-12 accepted",
                pass: tiny.positivity_pass,
            },
            UnitCheck {
                name: "out_of_tolerance_negative",
                observed: vfmt(negative.selected_minimum),
                expected: "-1e-6 rejected",
                pass: !negative.positivity_pass && !negative.solver_failure,
            },
            UnitCheck {
                name: "primary_failure_fallback_success",
                observed: selection.0.into(),
                expected: "complex_schur_fallback selected",
                pass: selection.0 == "complex_schur_fallback" && selection.3 && !selection.4,
            },
            UnitCheck {
                name: "both_solvers_fail",
                observed: both.4.to_string(),
                expected: "solver_failure true",
                pass: both.4,
            },
            UnitCheck {
                name: "nan_inf_formatter",
                observed: format!(
                    "{}|{}|{}",
                    vfmt(f64::NAN),
                    vfmt(f64::INFINITY),
                    vfmt(f64::NEG_INFINITY)
                ),
                expected: "NaN|+Inf|-Inf",
                pass: vfmt(f64::NAN) == "NaN"
                    && vfmt(f64::INFINITY) == "+Inf"
                    && vfmt(f64::NEG_INFINITY) == "-Inf",
            },
            UnitCheck {
                name: "eigenvalue_sum_trace_predicate",
                observed: format!(
                    "{}|{}|{}",
                    eigenvalue_sum_trace_pass(9.992e-16, 1e-10),
                    eigenvalue_sum_trace_pass(1e-8, 1e-10),
                    eigenvalue_sum_trace_pass(f64::NAN, 1e-10)
                ),
                expected: "true|false|false",
                pass: eigenvalue_sum_trace_pass(9.992e-16, 1e-10)
                    && !eigenvalue_sum_trace_pass(1e-8, 1e-10)
                    && !eigenvalue_sum_trace_pass(f64::NAN, 1e-10),
            },
            UnitCheck {
                name: "state_finite_independent_of_solver",
                observed: independent.to_string(),
                expected: "state finite remains true",
                pass: independent,
            },
            UnitCheck {
                name: "t002_regression",
                observed: format!(
                    "primary={} selected={} min={}",
                    regression.primary.pass,
                    regression.selected_solver,
                    vfmt(regression.selected_minimum)
                ),
                expected: "primary failure fallback success minimum near -3.237e-24",
                pass: !regression.primary.pass
                    && regression.selected_solver == "complex_schur_fallback"
                    && regression.positivity_pass
                    && (regression.selected_minimum - (-3.2372313550345495e-24)).abs() < 1e-20,
            },
        ])
    }

    fn write_unit_checks(checks: &[UnitCheck]) -> Result<(), Box<dyn std::error::Error>> {
        let mut w = BufWriter::new(File::create(V_OUTPUTS[0])?);
        writeln!(w, "check,observed,expected,status")?;
        for x in checks {
            writeln!(
                w,
                "{},{},{},{}",
                x.name,
                x.observed,
                x.expected,
                if x.pass { "PASS" } else { "FAIL" }
            )?;
        }
        Ok(())
    }

    fn validation_checks(
        run: &ValidationRun,
        summary: &Summary,
        construction: &[(String, String, String, bool)],
        comparison: &[ComparisonRow],
        units: &[UnitCheck],
    ) -> Vec<(String, String, String, String, bool)> {
        let mut out: Vec<_> = construction
            .iter()
            .map(|x| {
                (
                    "construction".into(),
                    x.0.clone(),
                    x.1.clone(),
                    x.2.clone(),
                    x.3,
                )
            })
            .collect();
        let (ps, pf, fa, fs, sf) = summary_counts(&run.eigen);
        let state_finite = run.rows.iter().all(|r| r.finite);
        let trace = run.rows.iter().map(|r| r.trace_error).fold(0.0, f64::max);
        let herm = run.rows.iter().map(|r| r.herm_error).fold(0.0, f64::max);
        let minimum = run
            .eigen
            .iter()
            .map(|x| x.selected_minimum)
            .fold(f64::INFINITY, f64::min);
        let max_imag = run
            .eigen
            .iter()
            .filter(|x| x.fallback.attempted)
            .map(|x| x.fallback.max_imag)
            .fold(0.0, f64::max);
        let max_sum = run
            .eigen
            .iter()
            .map(|x| x.selected_sum_trace_difference)
            .fold(0.0, f64::max);
        let ledger = run.rows.iter().map(|r| r.ledger.abs()).fold(0.0, f64::max);
        let reduced = run
            .rows
            .iter()
            .map(|r| r.reduced_trace_error)
            .fold(0.0, f64::max);
        let w_bound = run
            .rows
            .iter()
            .all(|r| r.work >= -1e-10 && r.work <= r.energy + 1e-10);
        let pop = run.rows.iter().all(|r| {
            r.sites
                .iter()
                .chain(r.load_populations.iter())
                .all(|x| *x >= -1e-10 && *x <= 1.0 + 1e-10)
        });
        let trajectory = comparison.iter().all(|x| x.pass);
        let unit_pass = units.iter().all(|x| x.pass);
        let rows = vec![
            (
                "unit",
                "required_unit_tests",
                unit_pass.to_string(),
                "true",
                unit_pass,
            ),
            (
                "execution",
                "saved_time_points",
                run.rows.len().to_string(),
                "1001",
                run.rows.len() == 1001,
            ),
            (
                "state",
                "all_physical_values_finite",
                state_finite.to_string(),
                "true",
                state_finite,
            ),
            ("state", "trace", vfmt(trace), "<=1e-8", trace <= TRACE_TOL),
            (
                "state",
                "hermiticity",
                vfmt(herm),
                "<=1e-8",
                herm <= HERM_TOL,
            ),
            (
                "solver",
                "primary_accounting",
                format!("success={ps} failure={pf}"),
                "sum=1001",
                ps + pf == 1001,
            ),
            (
                "solver",
                "fallback_success",
                format!("attempt={fa} success={fs}"),
                "all attempts successful",
                fa == fs,
            ),
            (
                "solver",
                "solver_failure_count",
                sf.to_string(),
                "0",
                sf == 0,
            ),
            (
                "solver",
                "selected_positivity",
                vfmt(minimum),
                ">=-1e-8",
                minimum >= -POS_TOL && run.eigen.iter().all(|x| x.positivity_pass),
            ),
            (
                "solver",
                "fallback_max_imag",
                vfmt(max_imag),
                "<=1e-10",
                max_imag <= EIG_IMAG_TOL,
            ),
            (
                "solver",
                "selected_sum_trace_difference",
                vfmt(max_sum),
                "finite and <=1e-10",
                eigenvalue_sum_trace_pass(max_sum, SUM_TRACE_TOL),
            ),
            (
                "state",
                "energy_ledger",
                vfmt(ledger),
                "<=5e-5",
                ledger <= LEDGER_TOL,
            ),
            (
                "state",
                "load_reduced_trace",
                vfmt(reduced),
                "<=1e-8",
                reduced <= TRACE_TOL,
            ),
            ("state", "W_le_E", w_bound.to_string(), "true", w_bound),
            ("state", "population_bounds", pop.to_string(), "true", pop),
            (
                "comparison",
                "saved_9c_trajectory",
                trajectory.to_string(),
                "all metrics <=1e-12",
                trajectory,
            ),
            (
                "summary",
                "endpoint_peak_classification",
                summary.peak_class.clone(),
                "finite classification",
                !summary.peak_class.is_empty(),
            ),
        ];
        out.extend(
            rows.into_iter()
                .map(|(a, b, c, d, e)| (a.into(), b.into(), c, d.into(), e)),
        );
        out
    }

    fn write_validation_checks_csv(
        rows: &[(String, String, String, String, bool)],
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut w = BufWriter::new(File::create(V_OUTPUTS[5])?);
        writeln!(w, "stage,check,observed,expected,status")?;
        for x in rows {
            writeln!(
                w,
                "{},{},{},{},{}",
                x.0,
                x.1,
                x.2,
                x.3,
                if x.4 { "PASS" } else { "FAIL" }
            )?;
        }
        Ok(())
    }

    fn at_time<'a>(rows: &'a [PositivityDiagnostic], time: f64) -> &'a PositivityDiagnostic {
        rows.iter().find(|x| (x.time - time).abs() < 1e-12).unwrap()
    }

    fn write_validation_report(
        run: &ValidationRun,
        summary: &Summary,
        ratio75: f64,
        ratio73: f64,
        final_status: &str,
        all_checks: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let (ps, pf, fa, fs, sf) = summary_counts(&run.eigen);
        let t1 = at_time(&run.eigen, 0.01);
        let t2 = at_time(&run.eigen, 0.02);
        let t3 = at_time(&run.eigen, 0.03);
        let worst = run
            .eigen
            .iter()
            .map(|x| x.selected_minimum)
            .fold(f64::INFINITY, f64::min);
        let max_imag = run
            .eigen
            .iter()
            .filter(|x| x.fallback.attempted)
            .map(|x| x.fallback.max_imag)
            .fold(0.0, f64::max);
        let max_sum = run
            .eigen
            .iter()
            .map(|x| x.selected_sum_trace_difference)
            .fold(0.0, f64::max);
        let sections=vec![
            ("1. 目的","N=7 fixed-total-noiseをrobust positivity診断でt=10まで再検証した。".into()),
            ("2. 元9cの停止理由","SymmetricEigenの非有限出力をstate finite失敗と混同しnumerical_issue_stopとなった。".into()),
            ("3. diagnosticで判明したこと","rhoは有限、traceとHermiticityは正常、Schurはt=0.02でminimum -3.237e-24を返した。".into()),
            ("4. 判定規則の問題","state finitenessとsolver finitenessを独立に判定する必要があった。".into()),
            ("5. robust solver policy","Hermitian化SymmetricEigenをprimaryとし、不合格時だけHermitian化Complex Schurへfallbackした。".into()),
            ("6. state finiteとsolver finiteの分離","primary失敗だけでは物理状態をFAILにせず、両solver失敗時だけsolver_failureとした。".into()),
            ("7. unit tests","必須10検査をunit testとruntime CSVで確認した。".into()),
            ("8. N=7再検証条件","N=7、total gamma=1.5、gamma_site=1.5/7、dt=0.0025、4000 RK4 steps、1001保存点。物理模型は不変。".into()),
            ("9. 実行時間",format!("total {:.3}s、propagation {:.3}s、diagnostics {:.3}s。",run.total_seconds,run.propagation_seconds,run.diagnostics_seconds)),
            ("10. primary solver結果",format!("成功{ps}時刻、失敗{pf}時刻。")),
            ("11. fallback結果",format!("attempt {fa}、success {fs}、両solver失敗 {sf}。max imag={:.3e}。",max_imag)),
            ("12. t=0.01診断",format!("selected={} minimum={} fallback={}。",t1.selected_solver,vfmt(t1.selected_minimum),t1.fallback_used)),
            ("13. t=0.02診断",format!("selected={} minimum={} fallback={}。",t2.selected_solver,vfmt(t2.selected_minimum),t2.fallback_used)),
            ("14. t=0.03診断",format!("selected={} minimum={} fallback={}。",t3.selected_solver,vfmt(t3.selected_minimum),t3.fallback_used)),
            ("15. 全1001時刻のpositivity",format!("worst selected minimum={}、solver_failure={}。",vfmt(worst),sf)),
            ("16. traceとHermiticity",format!("max trace={:.3e}、max Hermiticity={:.3e}。",summary.max_trace,summary.max_herm)),
            ("17. energy ledger",format!("max abs ledger={:.3e}。",summary.max_ledger)),
            ("18. 既存9c trajectoryとの一致",format!("全指定物理量を1001時刻、許容値1e-12で比較。all checks={}。",all_checks)),
            ("19. t=10結果",format!("E10={} W10={} usable={} W/Ein={}。",vfmt(summary.endpoint.energy),vfmt(summary.endpoint.work),vfmt(summary.endpoint.usable),vfmt(summary.endpoint.w_over_ein))),
            ("20. W最大値",format!("Wmax={} at t={}、E at Wmax={}、usable={}。",vfmt(summary.w_max.work),vfmt(summary.w_max.time),vfmt(summary.w_max.energy),vfmt(summary.w_max.usable))),
            ("21. N=3/N=5/N=7比較",format!("N7/N5={}、N7/N3={}。N7/N5は{}。",vfmt(ratio75),vfmt(ratio73),if ratio75<0.9{"N7 substantially lower"}else if ratio75<=1.1{"similar within 10 percent descriptive band"}else{"N7 substantially higher"})),
            ("22. fixed-per-site比較との違い","今回はtotal gammaを1.5へ固定した比較であり、siteごと0.5固定とは総雑音が異なる。".into()),
            ("23. 中心結果",format!("判定 **{final_status}**。fallbackは診断層だけで、物理時間発展を変更していない。")),
            ("24. 直接確認できたこと","N=7 t=10 trajectory、全保存点のselected positivity、fallback成功、既存9c物理量再現を確認した。".into()),
            ("25. 確認できていないこと","別dt、別gamma、t>10、N>7、別外部solver crateは未確認。".into()),
            ("26. 主張してはいけないこと","10%帯は統計的有意差でなく、距離だけの因果証明でもない。".into()),
            ("27. Milestone 9cの正式判定",format!("**{final_status}**。物理比較結果を正式採用{}。",if final_status=="completed_comparison_with_fallback_diagnostic"{"してよい"}else{"しない"})),
            ("28. 生成ファイル一覧",V_OUTPUTS.iter().map(|x|format!("- `{x}`")).collect::<Vec<_>>().join("\n")),
        ];
        let mut w = BufWriter::new(File::create(V_OUTPUTS[8])?);
        writeln!(
            w,
            "# Milestone 9c validation: robust positivity diagnostic\n"
        )?;
        for (s, b) in sections {
            writeln!(w, "## {s}\n\n{b}\n")?;
        }
        writeln!(w, "selected eigenvalue sum-trace maximum={:.3e}。", max_sum)?;
        Ok(())
    }

    pub fn run() -> Result<(), Box<dyn std::error::Error>> {
        for path in V_OUTPUTS {
            if std::path::Path::new(path).exists() {
                return Err(format!("refusing to overwrite {path}").into());
            }
        }
        for path in [
            "fixed_total_noise_timeseries.csv",
            "fixed_total_noise_summary.csv",
            "fixed_total_noise_length_comparison.csv",
        ] {
            if !std::path::Path::new(path).exists() {
                return Err(format!("missing reference {path}").into());
            }
        }
        let units = make_unit_checks()?;
        if let Some(x) = units.iter().find(|x| !x.pass) {
            return Err(format!("unit precheck failed {}", x.name).into());
        }
        write_unit_checks(&units)?;
        let (run, construction) = run_validation()?;
        let summary = summarize(&run.rows);
        let comparison = trajectory_comparison(&run.rows)?;
        write_trajectory_comparison(&comparison)?;
        let preliminary = validation_checks(&run, &summary, &construction, &comparison, &units);
        let all_checks = preliminary.iter().all(|x| x.4);
        let (_, _, _, _, solver_failures) = summary_counts(&run.eigen);
        let trajectory_ok = comparison.iter().all(|x| x.pass);
        let state_ok = run.rows.iter().all(|r| {
            r.finite
                && r.trace_error <= TRACE_TOL
                && r.herm_error <= HERM_TOL
                && r.ledger.abs() <= LEDGER_TOL
        });
        let final_status = if solver_failures > 0 {
            "solver_failure_stop"
        } else if run.eigen.iter().any(|x| x.selected_minimum < -POS_TOL) {
            "positivity_violation_stop"
        } else if !trajectory_ok {
            "trajectory_mismatch_stop"
        } else if !state_ok {
            "state_quality_stop"
        } else if all_checks {
            "completed_comparison_with_fallback_diagnostic"
        } else {
            "state_quality_stop"
        };
        write_validation_timeseries(&run.rows)?;
        write_eigen_diagnostics(&run.eigen)?;
        write_validation_summary(&run, &summary, final_status)?;
        write_validation_checks_csv(&preliminary)?;
        write_validation_performance(&run)?;
        let (r75, r73) = write_final_comparison(&summary)?;
        write_validation_report(&run, &summary, r75, r73, final_status, all_checks)?;
        println!(
            "completed status={final_status} checks={all_checks} E10={} W10={} Wmax={} fallback={}",
            vfmt(summary.endpoint.energy),
            vfmt(summary.endpoint.work),
            vfmt(summary.w_max.work),
            summary_counts(&run.eigen).2
        );
        if final_status != "completed_comparison_with_fallback_diagnostic" {
            return Err(format!("validation stop: {final_status}").into());
        }
        Ok(())
    }

    #[cfg(test)]
    mod robust_tests {
        use super::*;
        fn diag(values: &[f64]) -> PositivityDiagnostic {
            let m = ComplexMatrix::from_diagonal(&nalgebra::DVector::from_iterator(
                values.len(),
                values.iter().map(|x| C64::new(*x, 0.0)),
            ));
            evaluate_positivity(&m, 0.0)
        }
        #[test]
        fn normal_2x2() {
            let x = diag(&[0.25, 0.75]);
            assert_eq!(x.selected_solver, "symmetric_eigen");
            assert!((x.selected_minimum - 0.25).abs() < 1e-15);
        }
        #[test]
        fn pure_state() {
            let x = diag(&[1., 0.]);
            assert!(x.positivity_pass);
            assert!(x.selected_minimum.abs() < 1e-15);
        }
        #[test]
        fn tiny_negative_passes() {
            assert!(diag(&[1.0 + 1e-12, -1e-12]).positivity_pass);
        }
        #[test]
        fn large_negative_fails() {
            let x = diag(&[1.0 + 1e-6, -1e-6]);
            assert!(!x.positivity_pass);
            assert!(!x.solver_failure);
        }
        #[test]
        fn fallback_selection() {
            let p = SolverAttempt {
                attempted: true,
                all_finite: false,
                minimum: f64::NAN,
                max_imag: 0.,
                sum_trace_difference: f64::NAN,
                pass: false,
            };
            let f = SolverAttempt {
                attempted: true,
                all_finite: true,
                minimum: -1e-12,
                max_imag: 0.,
                sum_trace_difference: 1e-15,
                pass: true,
            };
            let x = select_attempts(p, f);
            assert_eq!(x.0, "complex_schur_fallback");
            assert!(x.3 && !x.4);
        }
        #[test]
        fn both_fail() {
            let p = SolverAttempt {
                attempted: true,
                all_finite: false,
                minimum: f64::NAN,
                max_imag: 0.,
                sum_trace_difference: f64::NAN,
                pass: false,
            };
            assert!(select_attempts(p.clone(), p).4);
        }
        #[test]
        fn formatter_distinguishes() {
            assert_eq!(vfmt(f64::NAN), "NaN");
            assert_eq!(vfmt(f64::INFINITY), "+Inf");
            assert_eq!(vfmt(f64::NEG_INFINITY), "-Inf");
        }
        #[test]
        fn sum_trace_predicate() {
            assert!(eigenvalue_sum_trace_pass(9.992e-16, 1e-10));
            assert!(!eigenvalue_sum_trace_pass(1e-8, 1e-10));
            assert!(!eigenvalue_sum_trace_pass(f64::NAN, 1e-10));
        }
        #[test]
        fn state_finite_independent() {
            let m = ComplexMatrix::identity(2, 2) * C64::new(0.5, 0.);
            assert!(all_elements_finite(&m));
            let failed = SolverAttempt {
                attempted: true,
                all_finite: false,
                minimum: f64::NAN,
                max_imag: 0.,
                sum_trace_difference: f64::NAN,
                pass: false,
            };
            assert!(!failed.pass);
        }
        #[test]
        fn t002_regression() {
            let p = ModelParams::default();
            let o = build_operators_for_chain(&p, N).unwrap();
            let k = DiagonalDephasingKernel::new(N, LOAD_DIM, &vec![GAMMA_SITE; N]).unwrap();
            let mut r = ComplexMatrix::zeros(DIM, DIM);
            r[(0, 0)] = C64::new(1., 0.);
            for s in 0..8 {
                r = rk4_step(&r, s as f64 * DT, &o, &k, GAMMA_SITE).unwrap();
            }
            let x = evaluate_positivity(&r, 0.02);
            assert!(!x.primary.pass);
            assert_eq!(x.selected_solver, "complex_schur_fallback");
            assert!(x.positivity_pass);
            assert!((x.selected_minimum + 3.2372313550345495e-24).abs() < 1e-20);
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    validation::run()
}
