use std::fs;
use std::path::Path;

#[derive(Clone, Debug)]
pub struct Config {
    pub cross_section_file: String,
    pub output_dir: String,
    pub length_m: f64,
    pub cells: usize,
    pub area_m2: f64,
    pub dt_s: f64,
    pub steps: usize,
    pub diagnostic_interval: usize,
    pub field_interval: usize,
    pub pressure_pa: f64,
    pub gas_temperature_k: f64,
    pub initial_density_m3: f64,
    pub initial_particles_per_species: usize,
    pub particle_weight: f64,
    pub initial_electron_temperature_ev: f64,
    pub initial_ion_temperature_ev: f64,
    pub left_voltage_dc_v: f64,
    pub right_voltage_dc_v: f64,
    pub rf_voltage_amplitude_v: f64,
    pub rf_frequency_hz: f64,
    pub ion_secondary_yield: f64,
    pub electron_reflection_probability: f64,
    pub electron_reflection_energy_factor: f64,
    pub secondary_electron_energy_ev: f64,
    pub enable_fn_emission: bool,
    pub fn_work_function_ev: f64,
    pub fn_field_enhancement: f64,
    pub fn_emission_area_m2: f64,
    pub fn_electron_energy_ev: f64,
    pub max_particles: usize,
    pub max_emitted_per_step: usize,
    pub rng_seed: u64,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            cross_section_file: "xsec/Cross section.txt".to_string(),
            output_dir: "output".to_string(),
            length_m: 0.03,
            cells: 128,
            area_m2: 1.0,
            dt_s: 2.0e-12,
            steps: 20_000,
            diagnostic_interval: 100,
            field_interval: 1000,
            pressure_pa: 20.0,
            gas_temperature_k: 300.0,
            initial_density_m3: 1.0e15,
            initial_particles_per_species: 20_000,
            particle_weight: 0.0,
            initial_electron_temperature_ev: 3.0,
            initial_ion_temperature_ev: 0.03,
            left_voltage_dc_v: -200.0,
            right_voltage_dc_v: 0.0,
            rf_voltage_amplitude_v: 0.0,
            rf_frequency_hz: 13.56e6,
            ion_secondary_yield: 0.10,
            electron_reflection_probability: 0.10,
            electron_reflection_energy_factor: 1.0,
            secondary_electron_energy_ev: 2.0,
            enable_fn_emission: true,
            fn_work_function_ev: 4.5,
            fn_field_enhancement: 1.0,
            fn_emission_area_m2: 1.0,
            fn_electron_energy_ev: 2.0,
            max_particles: 2_000_000,
            max_emitted_per_step: 100_000,
            rng_seed: 1,
        }
    }
}

impl Config {
    pub fn from_file(path: &Path) -> Result<Self, String> {
        let text = fs::read_to_string(path)
            .map_err(|e| format!("failed to read config {}: {e}", path.display()))?;
        let mut cfg = Self::default();
        cfg.apply_text(&text)?;
        cfg.validate()?;
        Ok(cfg)
    }

    pub fn write_default(path: &Path) -> Result<(), String> {
        fs::write(path, default_config_text())
            .map_err(|e| format!("failed to write config {}: {e}", path.display()))
    }

    pub fn macro_weight(&self) -> f64 {
        if self.particle_weight > 0.0 {
            self.particle_weight
        } else {
            self.initial_density_m3 * self.length_m * self.area_m2
                / self.initial_particles_per_species as f64
        }
    }

    pub fn left_voltage(&self, time_s: f64) -> f64 {
        self.left_voltage_dc_v
            + self.rf_voltage_amplitude_v
                * (std::f64::consts::TAU * self.rf_frequency_hz * time_s).sin()
    }

    pub fn right_voltage(&self, _time_s: f64) -> f64 {
        self.right_voltage_dc_v
    }

    fn apply_text(&mut self, text: &str) -> Result<(), String> {
        for (idx, raw) in text.lines().enumerate() {
            let mut line = raw;
            if let Some(pos) = line.find('#') {
                line = &line[..pos];
            }
            if let Some(pos) = line.find(';') {
                line = &line[..pos];
            }
            let line = line.trim();
            if line.is_empty() || line.starts_with('[') {
                continue;
            }
            let Some((key, value)) = line.split_once('=') else {
                return Err(format!("config line {} is not key=value: {raw}", idx + 1));
            };
            self.apply_pair(key.trim(), value.trim(), idx + 1)?;
        }
        Ok(())
    }

    fn apply_pair(&mut self, key: &str, value: &str, line: usize) -> Result<(), String> {
        let value = unquote(value);
        match key {
            "cross_section_file" => self.cross_section_file = value.to_string(),
            "output_dir" => self.output_dir = value.to_string(),
            "length_m" => self.length_m = parse_f64(key, value, line)?,
            "cells" => self.cells = parse_usize(key, value, line)?,
            "area_m2" => self.area_m2 = parse_f64(key, value, line)?,
            "dt_s" => self.dt_s = parse_f64(key, value, line)?,
            "steps" => self.steps = parse_usize(key, value, line)?,
            "diagnostic_interval" => self.diagnostic_interval = parse_usize(key, value, line)?,
            "field_interval" => self.field_interval = parse_usize(key, value, line)?,
            "pressure_pa" => self.pressure_pa = parse_f64(key, value, line)?,
            "gas_temperature_k" => self.gas_temperature_k = parse_f64(key, value, line)?,
            "initial_density_m3" => self.initial_density_m3 = parse_f64(key, value, line)?,
            "initial_particles_per_species" => {
                self.initial_particles_per_species = parse_usize(key, value, line)?
            }
            "particle_weight" => self.particle_weight = parse_f64(key, value, line)?,
            "initial_electron_temperature_ev" => {
                self.initial_electron_temperature_ev = parse_f64(key, value, line)?
            }
            "initial_ion_temperature_ev" => {
                self.initial_ion_temperature_ev = parse_f64(key, value, line)?
            }
            "left_voltage_dc_v" => self.left_voltage_dc_v = parse_f64(key, value, line)?,
            "right_voltage_dc_v" => self.right_voltage_dc_v = parse_f64(key, value, line)?,
            "rf_voltage_amplitude_v" => self.rf_voltage_amplitude_v = parse_f64(key, value, line)?,
            "rf_frequency_hz" => self.rf_frequency_hz = parse_f64(key, value, line)?,
            "ion_secondary_yield" => self.ion_secondary_yield = parse_f64(key, value, line)?,
            "electron_reflection_probability" => {
                self.electron_reflection_probability = parse_f64(key, value, line)?
            }
            "electron_reflection_energy_factor" => {
                self.electron_reflection_energy_factor = parse_f64(key, value, line)?
            }
            "secondary_electron_energy_ev" => {
                self.secondary_electron_energy_ev = parse_f64(key, value, line)?
            }
            "enable_fn_emission" => self.enable_fn_emission = parse_bool(key, value, line)?,
            "fn_work_function_ev" => self.fn_work_function_ev = parse_f64(key, value, line)?,
            "fn_field_enhancement" => self.fn_field_enhancement = parse_f64(key, value, line)?,
            "fn_emission_area_m2" => self.fn_emission_area_m2 = parse_f64(key, value, line)?,
            "fn_electron_energy_ev" => self.fn_electron_energy_ev = parse_f64(key, value, line)?,
            "max_particles" => self.max_particles = parse_usize(key, value, line)?,
            "max_emitted_per_step" => self.max_emitted_per_step = parse_usize(key, value, line)?,
            "rng_seed" => self.rng_seed = parse_u64(key, value, line)?,
            _ => return Err(format!("unknown config key on line {line}: {key}")),
        }
        Ok(())
    }

    pub fn validate(&self) -> Result<(), String> {
        require(self.length_m > 0.0, "length_m must be > 0")?;
        require(self.cells >= 2, "cells must be >= 2")?;
        require(self.area_m2 > 0.0, "area_m2 must be > 0")?;
        require(self.dt_s > 0.0, "dt_s must be > 0")?;
        require(
            self.diagnostic_interval > 0,
            "diagnostic_interval must be > 0",
        )?;
        require(self.field_interval > 0, "field_interval must be > 0")?;
        require(self.pressure_pa >= 0.0, "pressure_pa must be >= 0")?;
        require(
            self.gas_temperature_k > 0.0,
            "gas_temperature_k must be > 0",
        )?;
        require(
            self.initial_particles_per_species > 0,
            "initial_particles_per_species must be > 0",
        )?;
        require(
            self.initial_density_m3 >= 0.0,
            "initial_density_m3 must be >= 0",
        )?;
        require(self.particle_weight >= 0.0, "particle_weight must be >= 0")?;
        require(
            self.initial_electron_temperature_ev >= 0.0,
            "initial_electron_temperature_ev must be >= 0",
        )?;
        require(
            self.initial_ion_temperature_ev >= 0.0,
            "initial_ion_temperature_ev must be >= 0",
        )?;
        require(
            self.ion_secondary_yield >= 0.0,
            "ion_secondary_yield must be >= 0",
        )?;
        require(
            (0.0..=1.0).contains(&self.electron_reflection_probability),
            "electron_reflection_probability must be in [0, 1]",
        )?;
        require(
            (0.0..=1.0).contains(&self.electron_reflection_energy_factor),
            "electron_reflection_energy_factor must be in [0, 1]",
        )?;
        require(
            self.secondary_electron_energy_ev >= 0.0,
            "secondary_electron_energy_ev must be >= 0",
        )?;
        require(
            self.fn_work_function_ev > 0.0,
            "fn_work_function_ev must be > 0",
        )?;
        require(
            self.fn_field_enhancement > 0.0,
            "fn_field_enhancement must be > 0",
        )?;
        require(
            self.fn_emission_area_m2 >= 0.0,
            "fn_emission_area_m2 must be >= 0",
        )?;
        require(
            self.fn_electron_energy_ev >= 0.0,
            "fn_electron_energy_ev must be >= 0",
        )?;
        require(self.max_particles > 0, "max_particles must be > 0")?;
        require(
            self.initial_particles_per_species <= self.max_particles,
            "initial_particles_per_species must be <= max_particles",
        )?;
        Ok(())
    }
}

fn require(ok: bool, message: &str) -> Result<(), String> {
    if ok { Ok(()) } else { Err(message.to_string()) }
}

fn unquote(value: &str) -> &str {
    let value = value.trim();
    if value.len() >= 2
        && ((value.starts_with('"') && value.ends_with('"'))
            || (value.starts_with('\'') && value.ends_with('\'')))
    {
        &value[1..value.len() - 1]
    } else {
        value
    }
}

fn parse_f64(key: &str, value: &str, line: usize) -> Result<f64, String> {
    value
        .parse::<f64>()
        .map_err(|e| format!("invalid f64 for {key} on line {line}: {value} ({e})"))
}

fn parse_usize(key: &str, value: &str, line: usize) -> Result<usize, String> {
    value
        .parse::<usize>()
        .map_err(|e| format!("invalid usize for {key} on line {line}: {value} ({e})"))
}

fn parse_u64(key: &str, value: &str, line: usize) -> Result<u64, String> {
    value
        .parse::<u64>()
        .map_err(|e| format!("invalid u64 for {key} on line {line}: {value} ({e})"))
}

fn parse_bool(key: &str, value: &str, line: usize) -> Result<bool, String> {
    match value.to_ascii_lowercase().as_str() {
        "true" | "1" | "yes" | "on" => Ok(true),
        "false" | "0" | "no" | "off" => Ok(false),
        _ => Err(format!("invalid bool for {key} on line {line}: {value}")),
    }
}

pub fn default_config_text() -> &'static str {
    r#"# He PIC/MCC sample configuration.
# Lines are key=value. Units are SI unless the key name says _ev.

cross_section_file = "xsec/Cross section.txt"
output_dir = "output"

length_m = 0.03
cells = 128
area_m2 = 1.0
dt_s = 2.0e-12
steps = 20000
diagnostic_interval = 100
field_interval = 1000

pressure_pa = 20.0
gas_temperature_k = 300.0

initial_density_m3 = 1.0e15
initial_particles_per_species = 20000
# particle_weight <= 0 means auto: n0 * length * area / particles.
particle_weight = 0.0
initial_electron_temperature_ev = 3.0
initial_ion_temperature_ev = 0.03

left_voltage_dc_v = -200.0
right_voltage_dc_v = 0.0
rf_voltage_amplitude_v = 0.0
rf_frequency_hz = 13.56e6

ion_secondary_yield = 0.10
electron_reflection_probability = 0.10
electron_reflection_energy_factor = 1.0
secondary_electron_energy_ev = 2.0

enable_fn_emission = true
fn_work_function_ev = 4.5
fn_field_enhancement = 1.0
fn_emission_area_m2 = 1.0
fn_electron_energy_ev = 2.0

max_particles = 2000000
max_emitted_per_step = 100000
rng_seed = 1
"#
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_default_config() {
        let mut cfg = Config::default();
        cfg.apply_text(default_config_text()).unwrap();
        cfg.validate().unwrap();
        assert_eq!(cfg.cells, 128);
        assert!(cfg.macro_weight() > 0.0);
    }
}
