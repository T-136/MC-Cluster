use super::read_and_write;
use chemfiles::{Atom, Frame, Trajectory, UnitCell};
use std::collections::HashMap;

pub struct GridStructure {
    pub nn: HashMap<u32, [u32; super::CN], fnv::FnvBuildHasher>,
    pub nn_pair_no_intersec:
        HashMap<u64, [[u32; super::NN_PAIR_NO_INTERSEC_NUMBER]; 2], fnv::FnvBuildHasher>,
    pub xsites_positions: Vec<[f64; 3]>,
    pub unit_cell: [f64; 3],
}

impl GridStructure {
    pub fn new(
        (
            pairlist_file,
            // n_pairlist_file: String,
            nn_pair_no_int_file,
            atom_sites,
            bulk_file_name,
        ): (String, String, String, String),
    ) -> GridStructure {
        let nn = read_and_write::read_nn(&pairlist_file);
        let nn_pair_no_intersec = read_and_write::read_nn_pair_no_intersec(&nn_pair_no_int_file);

        // let bulk = Poscar::from_path(bulk_file_name).unwrap_or_else(|err| {
        //     panic!(
        //         "Could not parse '{:?}': {}",
        //         stringify!(bulk_file_name),
        //         err
        //     )
        // });

        let mut trajectory = Trajectory::open(bulk_file_name, 'r').unwrap();
        let mut frame = Frame::new();
        trajectory.read(&mut frame).unwrap();
        // frame.cell();
        let unit_cell = frame.cell().lengths();
        println!("unitcell {:?}", unit_cell);
        // let unit_cell = [
        //     unit_cell_size[0][0] * super::GRID_SIZE[0] as f64,
        //     unit_cell_size[1][1] * super::GRID_SIZE[1] as f64,
        //     unit_cell_size[2][2] * super::GRID_SIZE[2] as f64,
        // ];

        let xsites_positions = read_and_write::read_atom_sites(&atom_sites, nn.len() as u32);

        GridStructure {
            nn,
            nn_pair_no_intersec,
            xsites_positions,
            unit_cell,
        }
    }
}
