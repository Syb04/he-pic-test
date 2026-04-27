mod config;
mod rng;
mod sim;
mod xsec;

use crate::config::Config;
use crate::sim::Simulation;
use crate::xsec::CrossSections;
use std::env;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("error: {err}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<(), String> {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.iter().any(|a| a == "-h" || a == "--help") {
        print_help();
        return Ok(());
    }
    if args.first().map(String::as_str) == Some("--write-config") {
        let path = args
            .get(1)
            .map(PathBuf::from)
            .unwrap_or_else(|| PathBuf::from("config.ini"));
        Config::write_default(&path)?;
        println!("wrote {}", path.display());
        return Ok(());
    }

    let config_path = pick_config_path(args.first().map(String::as_str));
    let cfg = Config::from_file(&config_path)?;
    let xsec_path = PathBuf::from(&cfg.cross_section_file);
    let xs = CrossSections::load_lxcat_txt(&xsec_path)?;
    println!("config: {}", config_path.display());
    println!("cross sections: {}", xsec_path.display());
    for process in xs.electron.iter().chain(xs.ion.iter()) {
        println!(
            "  {:?}: {} (threshold {:.4} eV)",
            process.kind, process.name, process.threshold_ev
        );
    }

    let mut sim = Simulation::new(cfg, xs)?;
    println!("{}", sim.summary());
    let start = Instant::now();
    sim.run()?;
    println!("done in {:.3} s", start.elapsed().as_secs_f64());
    Ok(())
}

fn pick_config_path(arg: Option<&str>) -> PathBuf {
    if let Some(path) = arg {
        return PathBuf::from(path);
    }
    if Path::new("config.ini").exists() {
        PathBuf::from("config.ini")
    } else {
        PathBuf::from("config.example.ini")
    }
}

fn print_help() {
    println!(
        "He 1D PIC/MCC simulator\n\nUSAGE:\n  he_pic_mcc [config.ini]\n  he_pic_mcc --write-config [config.ini]\n\nIf no config is supplied, config.ini is used when present, otherwise config.example.ini."
    );
}
