use crate::config::Config;
use crate::rng::Rng;
use crate::xsec::{CollisionKind, CrossSection, CrossSections};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

const QE: f64 = 1.602_176_634e-19;
const EPS0: f64 = 8.854_187_8128e-12;
const KB: f64 = 1.380_649e-23;
const ME: f64 = 9.109_383_7015e-31;
const AMU: f64 = 1.660_539_066_60e-27;
const M_HE_PLUS: f64 = 4.002_602 * AMU;

#[derive(Default, Clone, Debug)]
struct IntervalStats {
    electron_elastic: u64,
    electron_excitation: u64,
    electron_ionization: u64,
    ion_backscatter: u64,
    ion_isotropic: u64,
    electron_absorbed_left: u64,
    electron_absorbed_right: u64,
    electron_reflected_left: u64,
    electron_reflected_right: u64,
    ion_absorbed_left: u64,
    ion_absorbed_right: u64,
    secondary_emitted: u64,
    fn_emitted_left: u64,
    fn_emitted_right: u64,
    emission_clipped: u64,
}

#[derive(Clone, Debug)]
struct Particles {
    x: Vec<f64>,
    vx: Vec<f64>,
    vy: Vec<f64>,
    vz: Vec<f64>,
}

impl Particles {
    fn with_capacity(capacity: usize) -> Self {
        Self {
            x: Vec::with_capacity(capacity),
            vx: Vec::with_capacity(capacity),
            vy: Vec::with_capacity(capacity),
            vz: Vec::with_capacity(capacity),
        }
    }

    fn len(&self) -> usize {
        self.x.len()
    }

    fn push(&mut self, x: f64, vx: f64, vy: f64, vz: f64) {
        self.x.push(x);
        self.vx.push(vx);
        self.vy.push(vy);
        self.vz.push(vz);
    }

    fn swap_remove(&mut self, index: usize) {
        self.x.swap_remove(index);
        self.vx.swap_remove(index);
        self.vy.swap_remove(index);
        self.vz.swap_remove(index);
    }

    fn speed2(&self, index: usize) -> f64 {
        self.vx[index] * self.vx[index]
            + self.vy[index] * self.vy[index]
            + self.vz[index] * self.vz[index]
    }

    fn energy_ev(&self, index: usize, mass_kg: f64) -> f64 {
        0.5 * mass_kg * self.speed2(index) / QE
    }

    fn set_velocity(&mut self, index: usize, vx: f64, vy: f64, vz: f64) {
        self.vx[index] = vx;
        self.vy[index] = vy;
        self.vz[index] = vz;
    }

    fn scale_velocity(&mut self, index: usize, scale: f64) {
        self.vx[index] *= scale;
        self.vy[index] *= scale;
        self.vz[index] *= scale;
    }

    fn mean_energy_ev(&self, mass_kg: f64) -> f64 {
        if self.len() == 0 {
            return 0.0;
        }
        let sum: f64 = (0..self.len()).map(|i| self.energy_ev(i, mass_kg)).sum();
        sum / self.len() as f64
    }
}

pub struct Simulation {
    cfg: Config,
    xs: CrossSections,
    rng: Rng,
    electrons: Particles,
    ions: Particles,
    rho: Vec<f64>,
    phi: Vec<f64>,
    efield: Vec<f64>,
    poisson_c: Vec<f64>,
    poisson_d: Vec<f64>,
    dx: f64,
    gas_density_m3: f64,
    neutral_temperature_ev: f64,
    macro_weight: f64,
    time_s: f64,
    stats: IntervalStats,
    diagnostics: BufWriter<File>,
    fields: BufWriter<File>,
}

impl Simulation {
    pub fn new(cfg: Config, xs: CrossSections) -> Result<Self, String> {
        fs::create_dir_all(&cfg.output_dir)
            .map_err(|e| format!("failed to create output dir {}: {e}", cfg.output_dir))?;
        let diagnostic_path = PathBuf::from(&cfg.output_dir).join("diagnostics.csv");
        let field_path = PathBuf::from(&cfg.output_dir).join("fields.csv");
        let diagnostics = BufWriter::new(
            File::create(&diagnostic_path)
                .map_err(|e| format!("failed to create {}: {e}", diagnostic_path.display()))?,
        );
        let fields = BufWriter::new(
            File::create(&field_path)
                .map_err(|e| format!("failed to create {}: {e}", field_path.display()))?,
        );

        let nodes = cfg.cells + 1;
        let initial_capacity = cfg.initial_particles_per_species.min(cfg.max_particles);
        let dx = cfg.length_m / cfg.cells as f64;
        let gas_density_m3 = cfg.pressure_pa / (KB * cfg.gas_temperature_k);
        let neutral_temperature_ev = KB * cfg.gas_temperature_k / QE;
        let macro_weight = cfg.macro_weight();

        let mut sim = Self {
            cfg: cfg.clone(),
            xs,
            rng: Rng::new(cfg.rng_seed),
            electrons: Particles::with_capacity(initial_capacity),
            ions: Particles::with_capacity(initial_capacity),
            rho: vec![0.0; nodes],
            phi: vec![0.0; nodes],
            efield: vec![0.0; nodes],
            poisson_c: vec![0.0; cfg.cells.saturating_sub(1)],
            poisson_d: vec![0.0; cfg.cells.saturating_sub(1)],
            dx,
            gas_density_m3,
            neutral_temperature_ev,
            macro_weight,
            time_s: 0.0,
            stats: IntervalStats::default(),
            diagnostics,
            fields,
        };
        sim.initialize_particles();
        Ok(sim)
    }

    pub fn run(&mut self) -> Result<(), String> {
        self.write_headers()?;
        self.solve_fields();
        self.write_fields(0)?;
        self.write_diagnostics(0)?;

        for step in 1..=self.cfg.steps {
            self.time_s = step as f64 * self.cfg.dt_s;
            self.solve_fields();
            if step % self.cfg.field_interval == 0 {
                self.write_fields(step)?;
            }

            self.push_electrons();
            self.push_ions();
            self.handle_electron_walls();
            self.handle_ion_walls();
            self.collide_electrons();
            self.collide_ions();
            self.emit_fowler_nordheim();

            if step % self.cfg.diagnostic_interval == 0 {
                self.write_diagnostics(step)?;
            }
        }

        self.diagnostics
            .flush()
            .map_err(|e| format!("failed to flush diagnostics: {e}"))?;
        self.fields
            .flush()
            .map_err(|e| format!("failed to flush fields: {e}"))?;
        Ok(())
    }

    pub fn summary(&self) -> String {
        format!(
            "cells={}, dx={:.4e} m, dt={:.4e} s, gas_density={:.4e} m^-3, macro_weight={:.4e}, electron_xsec={}, ion_xsec={}",
            self.cfg.cells,
            self.dx,
            self.cfg.dt_s,
            self.gas_density_m3,
            self.macro_weight,
            self.xs.electron.len(),
            self.xs.ion.len()
        )
    }

    fn initialize_particles(&mut self) {
        let n = self.cfg.initial_particles_per_species;
        for i in 0..n {
            let x = self.cfg.length_m * (i as f64 + self.rng.uniform()) / n as f64;
            let (vx, vy, vz) =
                maxwellian_velocity(self.cfg.initial_electron_temperature_ev, ME, &mut self.rng);
            self.electrons.push(x, vx, vy, vz);
            let (vx, vy, vz) = maxwellian_velocity(
                self.cfg.initial_ion_temperature_ev,
                M_HE_PLUS,
                &mut self.rng,
            );
            self.ions.push(x, vx, vy, vz);
        }
    }

    fn write_headers(&mut self) -> Result<(), String> {
        writeln!(
            self.diagnostics,
            "step,time_s,electrons_macro,ions_macro,ne_m3,ni_m3,mean_e_energy_ev,mean_i_energy_ev,left_voltage_v,right_voltage_v,electron_elastic,electron_excitation,electron_ionization,ion_backscatter,ion_isotropic,e_abs_left,e_abs_right,e_ref_left,e_ref_right,ion_abs_left,ion_abs_right,secondary_emitted,fn_left,fn_right,emission_clipped"
        )
        .map_err(|e| format!("failed to write diagnostics header: {e}"))?;
        writeln!(
            self.fields,
            "step,time_s,node,x_m,phi_v,efield_v_m,rho_c_m3"
        )
        .map_err(|e| format!("failed to write field header: {e}"))?;
        Ok(())
    }

    fn write_diagnostics(&mut self, step: usize) -> Result<(), String> {
        let volume = self.cfg.length_m * self.cfg.area_m2;
        let ne = self.electrons.len() as f64 * self.macro_weight / volume;
        let ni = self.ions.len() as f64 * self.macro_weight / volume;
        let s = self.stats.clone();
        writeln!(
            self.diagnostics,
            "{},{:.9e},{},{},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            step,
            self.time_s,
            self.electrons.len(),
            self.ions.len(),
            ne,
            ni,
            self.electrons.mean_energy_ev(ME),
            self.ions.mean_energy_ev(M_HE_PLUS),
            self.cfg.left_voltage(self.time_s),
            self.cfg.right_voltage(self.time_s),
            s.electron_elastic,
            s.electron_excitation,
            s.electron_ionization,
            s.ion_backscatter,
            s.ion_isotropic,
            s.electron_absorbed_left,
            s.electron_absorbed_right,
            s.electron_reflected_left,
            s.electron_reflected_right,
            s.ion_absorbed_left,
            s.ion_absorbed_right,
            s.secondary_emitted,
            s.fn_emitted_left,
            s.fn_emitted_right,
            s.emission_clipped
        )
        .map_err(|e| format!("failed to write diagnostics: {e}"))?;
        self.stats = IntervalStats::default();
        Ok(())
    }

    fn write_fields(&mut self, step: usize) -> Result<(), String> {
        for i in 0..=self.cfg.cells {
            writeln!(
                self.fields,
                "{},{:.9e},{},{:.9e},{:.9e},{:.9e},{:.9e}",
                step,
                self.time_s,
                i,
                i as f64 * self.dx,
                self.phi[i],
                self.efield[i],
                self.rho[i]
            )
            .map_err(|e| format!("failed to write fields: {e}"))?;
        }
        Ok(())
    }

    fn solve_fields(&mut self) {
        self.deposit_charge();
        self.solve_poisson();
        self.compute_electric_field();
    }

    fn deposit_charge(&mut self) {
        self.rho.fill(0.0);
        let electron_charge = -QE * self.macro_weight;
        let ion_charge = QE * self.macro_weight;
        deposit_species(
            &mut self.rho,
            &self.electrons,
            electron_charge,
            self.dx,
            self.cfg.cells,
            self.cfg.area_m2,
        );
        deposit_species(
            &mut self.rho,
            &self.ions,
            ion_charge,
            self.dx,
            self.cfg.cells,
            self.cfg.area_m2,
        );
    }

    fn solve_poisson(&mut self) {
        let n = self.cfg.cells;
        self.phi[0] = self.cfg.left_voltage(self.time_s);
        self.phi[n] = self.cfg.right_voltage(self.time_s);
        if n <= 1 {
            return;
        }

        let m = n - 1;
        let dx2_over_eps0 = self.dx * self.dx / EPS0;
        for row in 0..m {
            let i = row + 1;
            let a = if row == 0 { 0.0 } else { -1.0 };
            let c = if row + 1 == m { 0.0 } else { -1.0 };
            let mut d = self.rho[i] * dx2_over_eps0;
            if row == 0 {
                d += self.phi[0];
            }
            if row + 1 == m {
                d += self.phi[n];
            }
            let denom = if row == 0 {
                2.0
            } else {
                2.0 - a * self.poisson_c[row - 1]
            };
            self.poisson_c[row] = c / denom;
            self.poisson_d[row] = if row == 0 {
                d / denom
            } else {
                (d - a * self.poisson_d[row - 1]) / denom
            };
        }
        self.phi[n - 1] = self.poisson_d[m - 1];
        for row in (0..m - 1).rev() {
            let i = row + 1;
            self.phi[i] = self.poisson_d[row] - self.poisson_c[row] * self.phi[i + 1];
        }
    }

    fn compute_electric_field(&mut self) {
        let n = self.cfg.cells;
        self.efield[0] = -(self.phi[1] - self.phi[0]) / self.dx;
        for i in 1..n {
            self.efield[i] = -(self.phi[i + 1] - self.phi[i - 1]) / (2.0 * self.dx);
        }
        self.efield[n] = -(self.phi[n] - self.phi[n - 1]) / self.dx;
    }

    fn push_electrons(&mut self) {
        let qm_dt = -QE / ME * self.cfg.dt_s;
        for i in 0..self.electrons.len() {
            let e = interpolate_nodes(
                self.electrons.x[i],
                self.dx,
                self.cfg.cells,
                self.cfg.length_m,
                &self.efield,
            );
            self.electrons.vx[i] += qm_dt * e;
            self.electrons.x[i] += self.electrons.vx[i] * self.cfg.dt_s;
        }
    }

    fn push_ions(&mut self) {
        let qm_dt = QE / M_HE_PLUS * self.cfg.dt_s;
        for i in 0..self.ions.len() {
            let e = interpolate_nodes(
                self.ions.x[i],
                self.dx,
                self.cfg.cells,
                self.cfg.length_m,
                &self.efield,
            );
            self.ions.vx[i] += qm_dt * e;
            self.ions.x[i] += self.ions.vx[i] * self.cfg.dt_s;
        }
    }

    fn handle_electron_walls(&mut self) {
        let mut i = 0;
        let scale = self.cfg.electron_reflection_energy_factor.sqrt();
        while i < self.electrons.len() {
            let hit_left = self.electrons.x[i] < 0.0;
            let hit_right = self.electrons.x[i] > self.cfg.length_m;
            if hit_left || hit_right {
                if self.rng.uniform() < self.cfg.electron_reflection_probability {
                    if hit_left {
                        self.electrons.x[i] = wall_epsilon(self.dx, self.cfg.length_m);
                        self.electrons.vx[i] = self.electrons.vx[i].abs();
                        self.stats.electron_reflected_left += 1;
                    } else {
                        self.electrons.x[i] =
                            self.cfg.length_m - wall_epsilon(self.dx, self.cfg.length_m);
                        self.electrons.vx[i] = -self.electrons.vx[i].abs();
                        self.stats.electron_reflected_right += 1;
                    }
                    self.electrons.scale_velocity(i, scale);
                    i += 1;
                } else {
                    if hit_left {
                        self.stats.electron_absorbed_left += 1;
                    } else {
                        self.stats.electron_absorbed_right += 1;
                    }
                    self.electrons.swap_remove(i);
                }
            } else {
                i += 1;
            }
        }
    }

    fn handle_ion_walls(&mut self) {
        let mut i = 0;
        while i < self.ions.len() {
            let hit_left = self.ions.x[i] < 0.0;
            let hit_right = self.ions.x[i] > self.cfg.length_m;
            if hit_left || hit_right {
                if hit_left {
                    self.stats.ion_absorbed_left += 1;
                } else {
                    self.stats.ion_absorbed_right += 1;
                }
                self.ions.swap_remove(i);
                let emitted = self.rng.sample_count(self.cfg.ion_secondary_yield);
                for _ in 0..emitted {
                    if self.add_wall_electron(hit_left, self.cfg.secondary_electron_energy_ev) {
                        self.stats.secondary_emitted += 1;
                    }
                }
            } else {
                i += 1;
            }
        }
    }

    fn collide_electrons(&mut self) {
        let initial_len = self.electrons.len();
        for i in 0..initial_len {
            if i >= self.electrons.len() {
                break;
            }
            let energy_ev = self.electrons.energy_ev(i, ME);
            let speed = self.electrons.speed2(i).sqrt();
            if speed <= 0.0 {
                continue;
            }
            let Some((kind, threshold_ev, total_sigma)) =
                choose_collision(&self.xs.electron, energy_ev, &mut self.rng)
            else {
                continue;
            };
            let probability =
                1.0 - (-self.gas_density_m3 * speed * total_sigma * self.cfg.dt_s).exp();
            if self.rng.uniform() >= probability {
                continue;
            }

            match kind {
                CollisionKind::ElectronElastic => {
                    let (vx, vy, vz) = isotropic_velocity(energy_ev, ME, &mut self.rng);
                    self.electrons.set_velocity(i, vx, vy, vz);
                    self.stats.electron_elastic += 1;
                }
                CollisionKind::ElectronExcitation => {
                    let new_energy = (energy_ev - threshold_ev).max(0.0);
                    let (vx, vy, vz) = isotropic_velocity(new_energy, ME, &mut self.rng);
                    self.electrons.set_velocity(i, vx, vy, vz);
                    self.stats.electron_excitation += 1;
                }
                CollisionKind::ElectronIonization => {
                    let remaining = (energy_ev - threshold_ev).max(0.0);
                    let split = self.rng.uniform();
                    let primary_energy = remaining * split;
                    let secondary_energy = remaining - primary_energy;
                    let x = self.electrons.x[i];
                    let (vx, vy, vz) = isotropic_velocity(primary_energy, ME, &mut self.rng);
                    self.electrons.set_velocity(i, vx, vy, vz);
                    let (svx, svy, svz) = isotropic_velocity(secondary_energy, ME, &mut self.rng);
                    self.add_electron(x, svx, svy, svz);
                    let (ivx, ivy, ivz) =
                        maxwellian_velocity(self.neutral_temperature_ev, M_HE_PLUS, &mut self.rng);
                    self.add_ion(x, ivx, ivy, ivz);
                    self.stats.electron_ionization += 1;
                }
                CollisionKind::IonBackscatter | CollisionKind::IonIsotropic => {}
            }
        }
    }

    fn collide_ions(&mut self) {
        if self.xs.ion.is_empty() {
            return;
        }
        for i in 0..self.ions.len() {
            let energy_ev = self.ions.energy_ev(i, M_HE_PLUS);
            let speed = self.ions.speed2(i).sqrt();
            if speed <= 0.0 {
                continue;
            }
            let Some((kind, _threshold_ev, total_sigma)) =
                choose_collision(&self.xs.ion, energy_ev, &mut self.rng)
            else {
                continue;
            };
            let probability =
                1.0 - (-self.gas_density_m3 * speed * total_sigma * self.cfg.dt_s).exp();
            if self.rng.uniform() >= probability {
                continue;
            }

            match kind {
                CollisionKind::IonBackscatter => {
                    let (vx, vy, vz) =
                        maxwellian_velocity(self.neutral_temperature_ev, M_HE_PLUS, &mut self.rng);
                    self.ions.set_velocity(i, vx, vy, vz);
                    self.stats.ion_backscatter += 1;
                }
                CollisionKind::IonIsotropic => {
                    let (vx, vy, vz) = isotropic_velocity(energy_ev, M_HE_PLUS, &mut self.rng);
                    self.ions.set_velocity(i, vx, vy, vz);
                    self.stats.ion_isotropic += 1;
                }
                CollisionKind::ElectronElastic
                | CollisionKind::ElectronExcitation
                | CollisionKind::ElectronIonization => {}
            }
        }
    }

    fn emit_fowler_nordheim(&mut self) {
        if !self.cfg.enable_fn_emission {
            return;
        }
        if self.efield[0] < 0.0 {
            let count = self.fowler_nordheim_count(self.efield[0].abs());
            for _ in 0..count {
                if self.add_wall_electron(true, self.cfg.fn_electron_energy_ev) {
                    self.stats.fn_emitted_left += 1;
                }
            }
        }
        let n = self.cfg.cells;
        if self.efield[n] > 0.0 {
            let count = self.fowler_nordheim_count(self.efield[n].abs());
            for _ in 0..count {
                if self.add_wall_electron(false, self.cfg.fn_electron_energy_ev) {
                    self.stats.fn_emitted_right += 1;
                }
            }
        }
    }

    fn fowler_nordheim_count(&mut self, field_v_m: f64) -> usize {
        let current_density = fowler_nordheim_current_density(
            field_v_m,
            self.cfg.fn_work_function_ev,
            self.cfg.fn_field_enhancement,
        );
        let expected = current_density * self.cfg.fn_emission_area_m2 * self.cfg.dt_s
            / (QE * self.macro_weight);
        let count = self.rng.sample_count(expected);
        if count > self.cfg.max_emitted_per_step {
            self.stats.emission_clipped += (count - self.cfg.max_emitted_per_step) as u64;
            self.cfg.max_emitted_per_step
        } else {
            count
        }
    }

    fn add_wall_electron(&mut self, left_wall: bool, energy_ev: f64) -> bool {
        let x = if left_wall {
            wall_epsilon(self.dx, self.cfg.length_m)
        } else {
            self.cfg.length_m - wall_epsilon(self.dx, self.cfg.length_m)
        };
        let sign_x = if left_wall { 1.0 } else { -1.0 };
        let (vx, vy, vz) = hemisphere_velocity(energy_ev, ME, sign_x, &mut self.rng);
        self.add_electron(x, vx, vy, vz)
    }

    fn add_electron(&mut self, x: f64, vx: f64, vy: f64, vz: f64) -> bool {
        if self.electrons.len() < self.cfg.max_particles {
            self.electrons.push(x, vx, vy, vz);
            true
        } else {
            self.stats.emission_clipped += 1;
            false
        }
    }

    fn add_ion(&mut self, x: f64, vx: f64, vy: f64, vz: f64) {
        if self.ions.len() < self.cfg.max_particles {
            self.ions.push(x, vx, vy, vz);
        }
    }
}

fn choose_collision(
    processes: &[CrossSection],
    energy_ev: f64,
    rng: &mut Rng,
) -> Option<(CollisionKind, f64, f64)> {
    let total_sigma: f64 = processes.iter().map(|p| p.sigma(energy_ev)).sum();
    if total_sigma <= 0.0 || !total_sigma.is_finite() {
        return None;
    }
    let target = rng.uniform() * total_sigma;
    let mut acc = 0.0;
    for process in processes {
        acc += process.sigma(energy_ev);
        if target <= acc {
            return Some((process.kind, process.threshold_ev, total_sigma));
        }
    }
    processes
        .last()
        .map(|p| (p.kind, p.threshold_ev, total_sigma))
}

fn deposit_species(
    rho: &mut [f64],
    particles: &Particles,
    charge_c: f64,
    dx: f64,
    cells: usize,
    area_m2: f64,
) {
    let cell_volume = area_m2 * dx;
    for i in 0..particles.len() {
        let x = particles.x[i].clamp(0.0, dx * cells as f64);
        let s = (x / dx).min(cells as f64 - f64::EPSILON);
        let cell = s.floor() as usize;
        let f = s - cell as f64;
        rho[cell] += charge_c * (1.0 - f) / cell_volume;
        rho[cell + 1] += charge_c * f / cell_volume;
    }
}

fn interpolate_nodes(x: f64, dx: f64, cells: usize, length_m: f64, values: &[f64]) -> f64 {
    if x <= 0.0 {
        return values[0];
    }
    if x >= length_m {
        return values[cells];
    }
    let s = x / dx;
    let cell = s.floor() as usize;
    let f = s - cell as f64;
    values[cell] * (1.0 - f) + values[cell + 1] * f
}

fn maxwellian_velocity(temp_ev: f64, mass_kg: f64, rng: &mut Rng) -> (f64, f64, f64) {
    if temp_ev <= 0.0 {
        return (0.0, 0.0, 0.0);
    }
    let sigma = (QE * temp_ev / mass_kg).sqrt();
    (
        sigma * rng.normal(),
        sigma * rng.normal(),
        sigma * rng.normal(),
    )
}

fn isotropic_velocity(energy_ev: f64, mass_kg: f64, rng: &mut Rng) -> (f64, f64, f64) {
    let speed = (2.0 * QE * energy_ev.max(0.0) / mass_kg).sqrt();
    let mu = 2.0 * rng.uniform() - 1.0;
    let azimuth = std::f64::consts::TAU * rng.uniform();
    let radial = (1.0 - mu * mu).max(0.0).sqrt();
    (
        speed * mu,
        speed * radial * azimuth.cos(),
        speed * radial * azimuth.sin(),
    )
}

fn hemisphere_velocity(
    energy_ev: f64,
    mass_kg: f64,
    sign_x: f64,
    rng: &mut Rng,
) -> (f64, f64, f64) {
    let speed = (2.0 * QE * energy_ev.max(0.0) / mass_kg).sqrt();
    let mu = rng.uniform().sqrt();
    let azimuth = std::f64::consts::TAU * rng.uniform();
    let radial = (1.0 - mu * mu).max(0.0).sqrt();
    (
        sign_x * speed * mu,
        speed * radial * azimuth.cos(),
        speed * radial * azimuth.sin(),
    )
}

fn fowler_nordheim_current_density(field_v_m: f64, work_function_ev: f64, beta: f64) -> f64 {
    if field_v_m <= 0.0 || beta <= 0.0 || work_function_ev <= 0.0 {
        return 0.0;
    }
    let local_field = beta * field_v_m;
    let exponent = -6.830_890e9 * work_function_ev.powf(1.5) / local_field;
    if exponent < -745.0 {
        return 0.0;
    }
    1.541_434e-6 * local_field * local_field / work_function_ev * exponent.exp()
}

fn wall_epsilon(dx: f64, length_m: f64) -> f64 {
    (1.0e-6 * dx).max(1.0e-12 * length_m)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fn_current_increases_with_field() {
        let j1 = fowler_nordheim_current_density(1.0e8, 4.5, 1.0);
        let j2 = fowler_nordheim_current_density(2.0e8, 4.5, 1.0);
        assert!(j2 > j1);
    }

    #[test]
    fn node_interpolation_mid_cell() {
        let values = [0.0, 10.0, 20.0];
        let x = interpolate_nodes(0.25, 0.5, 2, 1.0, &values);
        assert!((x - 5.0).abs() < 1.0e-12);
    }
}
