use std::{
    fs::{self, read_to_string},
    io,
    path::Path,
};

use chemrust_core::data::{Atom, LatticeModel};
use chemrust_formats::{
    castep_param::{BandStructureParam, GeomOptParam},
    seed_writer::{MetalMethodsControl, SeedWriter},
    Cell, StructureFile,
};
use chemrust_parser::CellParser;
use nalgebra::Translation3;

fn translation_matrix() -> Translation3<f64> {
    Translation3::new(2.67872 - 0.707254, 9.97325 - 8.57500, 5.931993 - 1.01750)
}

fn combine_ni_sinx(ni_lattice: &LatticeModel, si_lattice: &LatticeModel) -> LatticeModel {
    let copy_of_ni_lattice = ni_lattice.clone();
    let mut ni_atoms: Vec<Atom> = copy_of_ni_lattice.atoms().to_vec();
    ni_atoms.iter_mut().for_each(|atom| {
        atom.set_cartesian_coord(translation_matrix().transform_point(&atom.cartesian_coord()))
    });
    let mut new_lattice = si_lattice.clone();
    new_lattice.append_atom(&mut ni_atoms);
    new_lattice
}

pub fn generate_combined_cell<P: AsRef<Path>>(
    si_lattice: &LatticeModel,
    ni_cell_path: &P,
) -> StructureFile<Cell> {
    let cell_text = read_to_string(ni_cell_path).unwrap();
    let ni_lattice = CellParser::new(&cell_text)
        .to_lattice_cart()
        .to_positions()
        .build_lattice();
    let new_lattice = combine_ni_sinx(&ni_lattice, si_lattice);
    StructureFile::new(new_lattice)
}

pub fn combined_cell_name<P: AsRef<Path>>(ni_cell_path: &P) -> String {
    let cell_name = ni_cell_path.as_ref().file_stem().unwrap().to_str().unwrap();
    format!("{cell_name}_SiNx_am")
}

pub fn generate_seed_file(
    cell_file: StructureFile<Cell>,
    cell_name: &str,
    export_loc: &str,
    potential_loc: &str,
) -> Result<(), io::Error> {
    let geom_seed_writer = SeedWriter::<GeomOptParam>::build(cell_file)
        .with_seed_name(cell_name)
        .with_potential_loc(potential_loc)
        .with_export_loc(export_loc)
        .build_edft(false);
    geom_seed_writer.write_seed_files()?;
    copy_smcastep_extension(&geom_seed_writer)?;
    let bs_writer: SeedWriter<BandStructureParam> = geom_seed_writer.into();
    bs_writer.write_seed_files()?;
    Ok(())
}

fn copy_smcastep_extension(writer: &SeedWriter<GeomOptParam>) -> Result<(), io::Error> {
    let dest_dir = writer.create_export_dir()?;
    let with_seed_name = format!("SMCastep_Extension_{}.xms", writer.seed_name());
    let dest_path = dest_dir.join(&with_seed_name);
    let cwd = env!("CARGO_MANIFEST_DIR");
    if !dest_path.exists() {
        fs::copy(&format!("{}/SMCastep_Extension.xms", cwd), dest_path)?;
    }
    Ok(())
}
