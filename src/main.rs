mod combiner;

use std::{fs::read_to_string, io, path::PathBuf};

use chemrust_misctools::{copy_potentials, to_xsd_scripts, write_lsf_script};
use chemrust_parser::CellParser;
use combiner::{combined_cell_name, generate_combined_cell, generate_seed_file};
use glob::glob;

fn main() -> Result<(), io::Error> {
    let hcp_export_loc = "hcp_SiNx_growth";
    let potential_loc = "./Potentials";
    let si_lattice = CellParser::new(&read_to_string("./SiNx_am/SiNx_am.cell").unwrap())
        .to_lattice_cart()
        .to_positions()
        .build_lattice();
    glob("hcp_growth/**/*.cell")
        .expect("Glob pattern not found")
        .into_iter()
        .try_for_each(|cell_path| {
            let cell_path = cell_path.unwrap();
            if cell_path
                .file_stem()
                .unwrap()
                .to_str()
                .unwrap()
                .contains("DOS")
            {
                return Ok(());
            }
            let combined_name = combined_cell_name(&cell_path);
            let combined_cell = generate_combined_cell(&si_lattice, &cell_path);
            generate_seed_file(combined_cell, &combined_name, hcp_export_loc, potential_loc)
        })?;
    let fcc_export_loc = "fcc_SiNx_growth";
    glob("fcc_growth/**/*.cell")
        .expect("Glob pattern not found")
        .into_iter()
        .try_for_each(|cell_path| {
            let cell_path = cell_path.unwrap();
            if cell_path
                .file_stem()
                .unwrap()
                .to_str()
                .unwrap()
                .contains("DOS")
            {
                return Ok(());
            }
            let combined_name = combined_cell_name(&cell_path);
            let combined_cell = generate_combined_cell(&si_lattice, &cell_path);
            generate_seed_file(combined_cell, &combined_name, fcc_export_loc, potential_loc)
        })?;
    glob("*_SiNx_growth/**/*.cell")
        .expect("Glob pattern not found")
        .into_iter()
        .try_for_each(|cell_path| {
            let cell_path = cell_path.unwrap();
            if cell_path
                .file_stem()
                .unwrap()
                .to_str()
                .unwrap()
                .contains("DOS")
            {
                return Ok(());
            }
            write_lsf_script(&cell_path, 18)?;
            copy_potentials(&cell_path, &PathBuf::new().join(potential_loc))
        })?;
    to_xsd_scripts("*_SiNx_growth").unwrap();
    Ok(())
}
