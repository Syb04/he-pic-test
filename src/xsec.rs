use std::fs;
use std::path::Path;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CollisionKind {
    ElectronElastic,
    ElectronExcitation,
    ElectronIonization,
    IonBackscatter,
    IonIsotropic,
}

#[derive(Clone, Debug)]
pub struct CrossSection {
    pub kind: CollisionKind,
    pub name: String,
    pub threshold_ev: f64,
    energy_ev: Vec<f64>,
    sigma_m2: Vec<f64>,
}

#[derive(Clone, Debug)]
pub struct CrossSections {
    pub electron: Vec<CrossSection>,
    pub ion: Vec<CrossSection>,
}

impl CrossSection {
    pub fn sigma(&self, energy_ev: f64) -> f64 {
        if self.energy_ev.is_empty() || !energy_ev.is_finite() {
            return 0.0;
        }
        if energy_ev <= self.energy_ev[0] {
            return if self.energy_ev[0] == 0.0 {
                self.sigma_m2[0]
            } else {
                0.0
            };
        }
        let last = self.energy_ev.len() - 1;
        if energy_ev >= self.energy_ev[last] {
            return self.sigma_m2[last];
        }
        match self
            .energy_ev
            .binary_search_by(|probe| probe.total_cmp(&energy_ev))
        {
            Ok(i) => self.sigma_m2[i],
            Err(i) => {
                let lo = i - 1;
                let hi = i;
                let x0 = self.energy_ev[lo];
                let x1 = self.energy_ev[hi];
                let y0 = self.sigma_m2[lo];
                let y1 = self.sigma_m2[hi];
                if x1 == x0 {
                    y0
                } else {
                    y0 + (y1 - y0) * (energy_ev - x0) / (x1 - x0)
                }
            }
        }
    }
}

impl CrossSections {
    pub fn load_lxcat_txt(path: &Path) -> Result<Self, String> {
        let text = fs::read_to_string(path)
            .map_err(|e| format!("failed to read cross sections {}: {e}", path.display()))?;
        let mut electron = Vec::new();
        let mut ion = Vec::new();
        let mut header: Vec<String> = Vec::new();
        let mut data: Vec<(f64, f64)> = Vec::new();
        let mut in_table = false;

        for line in text.lines() {
            if is_dash_line(line) {
                if in_table {
                    if let Some(process) = build_process(&header, &data) {
                        match process.kind {
                            CollisionKind::ElectronElastic
                            | CollisionKind::ElectronExcitation
                            | CollisionKind::ElectronIonization => electron.push(process),
                            CollisionKind::IonBackscatter | CollisionKind::IonIsotropic => {
                                ion.push(process)
                            }
                        }
                    }
                    data.clear();
                    header.clear();
                    in_table = false;
                } else {
                    in_table = true;
                    data.clear();
                }
                continue;
            }

            if in_table {
                if let Some(pair) = parse_energy_sigma(line) {
                    data.push(pair);
                }
            } else {
                header.push(line.to_string());
            }
        }

        if electron.is_empty() {
            return Err("no electron cross sections found in LXCat txt file".to_string());
        }
        Ok(Self { electron, ion })
    }
}

fn build_process(header: &[String], data: &[(f64, f64)]) -> Option<CrossSection> {
    if data.is_empty() {
        return None;
    }
    let process_line = header
        .iter()
        .rev()
        .find(|line| line.trim_start().starts_with("PROCESS:"))?;
    let species_line = header
        .iter()
        .rev()
        .find(|line| line.trim_start().starts_with("SPECIES:"))
        .map(String::as_str)
        .unwrap_or("");
    let process = process_line.trim();
    let lower = process.to_ascii_lowercase();
    let species = species_line.to_ascii_lowercase();

    let kind = if species.contains("e / he") || process.contains("E + He") {
        if lower.contains("ionization") {
            CollisionKind::ElectronIonization
        } else if lower.contains("excitation") {
            CollisionKind::ElectronExcitation
        } else if lower.contains("elastic") {
            CollisionKind::ElectronElastic
        } else {
            return None;
        }
    } else if species.contains("he^+ / he") || process.contains("He+ + He") {
        if lower.contains("backscat") {
            CollisionKind::IonBackscatter
        } else if lower.contains("isotropic") {
            CollisionKind::IonIsotropic
        } else {
            return None;
        }
    } else {
        return None;
    };

    let threshold_ev = match kind {
        CollisionKind::ElectronExcitation | CollisionKind::ElectronIonization => {
            parse_threshold_ev(header).unwrap_or(data[0].0)
        }
        _ => 0.0,
    };
    let (energy_ev, sigma_m2): (Vec<_>, Vec<_>) = data.iter().copied().unzip();
    let name = process
        .strip_prefix("PROCESS:")
        .unwrap_or(process)
        .trim()
        .to_string();

    Some(CrossSection {
        kind,
        name,
        threshold_ev,
        energy_ev,
        sigma_m2,
    })
}

fn parse_threshold_ev(header: &[String]) -> Option<f64> {
    for line in header.iter().rev() {
        let trimmed = line.trim_start();
        if trimmed.starts_with("PARAM.:") && trimmed.contains("E =") {
            return extract_numbers(trimmed).into_iter().next();
        }
    }
    None
}

fn is_dash_line(line: &str) -> bool {
    let trimmed = line.trim();
    trimmed.len() >= 5 && trimmed.chars().all(|c| c == '-')
}

fn parse_energy_sigma(line: &str) -> Option<(f64, f64)> {
    let nums = extract_numbers(line);
    if nums.len() >= 2 {
        Some((nums[0], nums[1]))
    } else {
        None
    }
}

fn extract_numbers(text: &str) -> Vec<f64> {
    let mut out = Vec::new();
    let mut token = String::new();
    for ch in text.chars().chain(std::iter::once(' ')) {
        if ch.is_ascii_digit() || matches!(ch, '.' | '+' | '-' | 'e' | 'E') {
            token.push(ch);
        } else {
            if token.chars().any(|c| c.is_ascii_digit()) {
                if let Ok(value) = token.parse::<f64>() {
                    out.push(value);
                }
            }
            token.clear();
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_provided_lxcat_file() {
        let xs = CrossSections::load_lxcat_txt(Path::new("xsec/Cross section.txt")).unwrap();
        assert_eq!(xs.electron.len(), 4);
        assert_eq!(xs.ion.len(), 2);
        assert!(
            xs.electron
                .iter()
                .any(|p| p.kind == CollisionKind::ElectronIonization)
        );
        assert!(
            xs.ion
                .iter()
                .any(|p| p.kind == CollisionKind::IonBackscatter)
        );
    }

    #[test]
    fn interpolation_is_finite() {
        let xs = CrossSections::load_lxcat_txt(Path::new("xsec/Cross section.txt")).unwrap();
        let sigma = xs.electron[0].sigma(10.0);
        assert!(sigma.is_finite());
        assert!(sigma > 0.0);
    }
}
