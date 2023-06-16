use anyhow;
use chemfiles::{Atom, Frame, Trajectory, UnitCell};
use rand;
use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use rand_chacha::ChaCha20Rng;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::{cmp, eprint, fs, println};
use vasp_poscar::Poscar;

mod listdict;
mod read_files;
mod setup;
mod sim;
mod sim_struct;

pub use sim::Results;

pub fn find_simulation_with_lowest_energy(folder: String) -> anyhow::Result<()> {
    // let mut folder_with_lowest_e: PathBuf = PathBuf::new();
    let mut lowest_e: f64 = f64::INFINITY;

    for _ in 0..2 {
        let paths = fs::read_dir(&folder).unwrap();
        for path in paths {
            let ok_path = match path {
                Ok(ok_path) => ok_path,
                Err(e) => {
                    eprintln!("{:?}", e);
                    continue;
                }
            };

            if !ok_path.path().is_dir() {
                continue;
            }
            let folder = fs::read_dir(ok_path.path().as_path())?;
            for folder_entry in folder {
                let file = match folder_entry {
                    Ok(path2) => {
                        if !path2.path().is_dir() {
                            path2
                        } else {
                            println!("unexpected folder");
                            continue;
                        }
                    }
                    Err(e) => {
                        eprintln!("{:?}", e);
                        continue;
                    }
                };
                if !file.path().to_str().unwrap().ends_with(".json") {
                    continue;
                }

                let file = fs::File::open(file.path()).unwrap();
                let reader = BufReader::new(file);
                let res: Result<Results, serde_json::Error> = serde_json::from_reader(reader);
                match res {
                    Ok(res) => {
                        if res.lowest_energy_struct.energy <= lowest_e {
                            // folder_with_lowest_e = path2.path();
                            lowest_e = res.lowest_energy_struct.energy;
                        } else {
                            println!("{:?}", res.lowest_energy_struct.energy);
                            println!("{:?}", ok_path.path())
                            // fs::remove_dir_all(ok_path.path())?;
                        }
                    }
                    Err(e) => {
                        eprintln!("{:?}", format!("{:?} in folder {:?}", e, ok_path.path()));
                        // fs::remove_dir_all(ok_path.path())?;
                    }
                }
            }
        }
    }

    Ok(())
}

// #[derive(Clone)]
// pub struct Simulation {
//     niter: u64,
//     number_all_atoms: u32,
//     occ: Vec<u8>,
//     onlyocc: HashSet<u32>,
//     cn: Vec<usize>,
//     former_energy_dict: HashMap<u32, i64>,
//     possible_moves: listdict::ListDict,
//     total_energy_1000: i64,
//     nn: HashMap<u32, [u32; 12], fnv::FnvBuildHasher>,
//     nn_pair: HashMap<u64, [u32; 20], fnv::FnvBuildHasher>,
//     // nnn_pair: HashMap<u64, [u32; 74], fnv::FnvBuildHasher>,
//     // mut xyz: Vec<[f64; 3]>,
//     xsites_positions: Vec<[f64; 3]>,
//     unit_cell: UnitCell,
//     cn_dict: [u32; 13],
//     save_folder: String,
//     trajectory_frequency: Option<u64>,
//     last_frames_trajectory: Option<u64>,
//     start_temperature: Option<f64>,
//     temperature: f64,
//     cn_dict_sections: Vec<HashMap<u8, f64>>,
//     energy_sections_list: Vec<f64>,
//     optimization_cut_off_perc: f64,
//     unique_levels: HashMap<BTreeMap<u8, u32>, (i64, u64)>,
// }

impl sim_struct::Simulation {
    pub fn new(
        niter: u64,
        nsites: u32,
        input_file: Option<String>,
        atoms_input: Option<u32>,
        temperature: f64,
        start_temperature: Option<f64>,
        save_folder_name: String,
        pairlist_file: String,
        nn_pairlist_file: String,
        // nnn_pairlist_file: String,
        atom_sites: String,
        trajectory_frequency: Option<u64>,
        last_frames_trajectory: Option<u64>,
        bulk_file_name: String,
        repetition: usize,
        optimization_cut_off_perc: f64,
    ) -> sim_struct::Simulation {
        let nn = read_files::read_nn(&pairlist_file);
        let nn_pair = read_files::read_nn_pairlists(&nn_pairlist_file);
        // let nnn_pair = read_files::read_nnn_pairlists(&nnn_pairlist_file);

        let bulk = Poscar::from_path(bulk_file_name).unwrap_or_else(|err| {
            panic!(
                "Could not parse '{:?}': {}",
                stringify!(bulk_file_name),
                err
            )
        });
        let unit_cell_size = bulk.unscaled_lattice_vectors();
        let unit_cell = UnitCell::new([
            unit_cell_size[0][0] * 15.,
            unit_cell_size[1][1] * 15.,
            unit_cell_size[2][2] * 15.,
        ]);
        let mut cn_dict: [u32; 13] = [0; 13];

        let xsites_positions = read_files::read_atom_sites(&atom_sites, nsites);
        let (occ, onlyocc, number_all_atoms) = if input_file.is_some() {
            let xyz = read_files::read_sample(&input_file.unwrap());
            let (occ, onlyocc) = setup::occ_onlyocc_from_xyz(&xyz, nsites, &xsites_positions);
            let number_of_atoms: u32 = onlyocc.len() as u32;
            (occ, onlyocc, number_of_atoms)
        } else if atoms_input.is_some() {
            let number_of_atom = atoms_input.unwrap();
            let (occ, onlyocc) =
                setup::create_input_cluster(&atoms_input.unwrap(), &xsites_positions, &nn, nsites);
            (occ, onlyocc, number_of_atom)
        } else {
            panic!("gib input atoms or input file");
        };
        let mut cn: Vec<usize> = Vec::with_capacity(nsites as usize);
        for o in 0..nsites {
            let mut neighbors: u8 = 0;
            for o1 in nn[&o].iter() {
                if occ[*o1 as usize] == 1 {
                    // cn.entry(o).and_modify(|x| *x += 1).or_insert(1);
                    neighbors += 1;
                }
            }
            cn.push(neighbors as usize);
            if occ[o as usize] == 1 {
                cn_dict[cn[o as usize] as usize] += 1;
            };
        }
        let mut former_energy_dict: HashMap<u32, i64> = HashMap::with_capacity(nsites as usize);
        let mut total_energy_1000: i64 = 0;
        let mut possible_moves: listdict::ListDict = listdict::ListDict::new();
        for o in onlyocc.iter() {
            let energy_1000: i64 = sim::energy_calculation(o, &cn);
            total_energy_1000 += energy_1000;
            former_energy_dict.insert(o.clone(), energy_1000);

            for u in &nn[o] {
                if occ[*u as usize] == 0 {
                    // if possible_moves.contains(o.clone(), u.clone()) {
                    //     continue;
                    // }
                    if cn[*o as usize] < 12 && cn[*u as usize] > 1 {
                        possible_moves.add_item(o.clone(), u.clone())
                    }
                }
            }
        }

        let simulation_folder_name = if let Some(start_temp) = start_temperature {
            std::format!(
                "{}-{}K_{}I_{}A",
                start_temp,
                temperature,
                niter,
                onlyocc.len()
            )
        } else {
            std::format!("{}K_{}I_{}A", temperature, niter, onlyocc.len())
        };

        let mut sub_folder = save_folder_name + &simulation_folder_name;

        sub_folder = std::format!("{}_{}", sub_folder, repetition);

        fs::create_dir(&sub_folder).unwrap_or_else(|error| {
            println!(
                "could not create folder: {} with error: {}",
                sub_folder, error
            )
        });

        let cn_dict_sections = Vec::new();
        let energy_sections_list = Vec::new();
        let unique_levels = HashMap::new();

        sim_struct::Simulation {
            niter,
            number_all_atoms,
            occ,
            onlyocc,
            cn,
            former_energy_dict,
            possible_moves,
            total_energy_1000,
            nn,
            nn_pair,
            // nnn_pair,
            xsites_positions,
            unit_cell,
            cn_dict,
            save_folder: sub_folder,
            trajectory_frequency,
            last_frames_trajectory,
            start_temperature,
            temperature,
            cn_dict_sections,
            energy_sections_list,
            optimization_cut_off_perc,
            unique_levels,
        }
    }

    pub fn run(&mut self, mut amount_unique_levels: i32) -> Results {
        let mut rng_choose = ChaCha20Rng::from_entropy();
        let choose_seed: [u8; 32] = rng_choose.get_seed();

        let mut rng_e_number = ChaCha20Rng::from_entropy();
        let e_number_seed: [u8; 32] = rng_e_number.get_seed();

        let seed = sim::Seed {
            rust: "used rust".to_string(),
            choose_seed,
            e_number_seed,
        };

        let mut trajectory: Option<Trajectory>;

        if let Some(_) = self.trajectory_frequency {
            trajectory = Some(
                Trajectory::open(self.save_folder.clone() + "/total_time_traj.xyz", 'w').unwrap(),
            );
        } else {
            trajectory = None;
        }

        let mut trajectory_last_frames: Option<Trajectory>;
        if let Some(i) = self.last_frames_trajectory {
            trajectory_last_frames = Some(
                Trajectory::open(
                    self.save_folder.clone() + &format!("/last_{i}_frames_trajectory.xyz"),
                    'w',
                )
                .unwrap(),
            );
        } else {
            trajectory_last_frames = None;
        }

        let mut lowest_energy_struct: sim::LowestEnergy = sim::LowestEnergy {
            energy: f64::INFINITY,
            cn_total: HashMap::new(),
            empty_cn: HashMap::new(),
            iiter: 0,
        };
        let mut temp_energy_section: i64 = 0;
        let mut temp_cn_dict_section: [u32; 13] = [0; 13];

        // for k in 1..13 {
        //     temp_cn_dict_section.insert(k, 0);
        // }

        let start_struct = self.create_start_struct();

        if self.niter == 0 {
            self.save_lowest_energy(
                &0,
                &mut lowest_energy_struct,
                // &mut trajectory_lowest_energy,
            );
        }

        for iiter in 0..self.niter {
            if iiter % 100000 == 0 {
                println!(
                    "iteration {}; {}%",
                    iiter,
                    (iiter as f64 / self.niter as f64 * 100.)
                );
            }
            let (move_from, move_to) = self.possible_moves.choose_random_item(&mut rng_choose);

            self.perform_move(move_from, move_to);

            let mut total_temp_energy: i64 = self.total_energy_1000.clone();

            self.temp_energy_calculation(move_from, move_to, &mut total_temp_energy);

            if self.is_acceptance_criteria_fulfilled(total_temp_energy, &mut rng_e_number, iiter) {
                self.accept_move(total_temp_energy, move_from, move_to);
            } else {
                self.perform_move(move_to, move_from);
            }

            self.save_lowest_energy(
                &iiter,
                &mut lowest_energy_struct,
                // &mut trajectory_lowest_energy,
            );

            self.write_trajectorys(
                &iiter,
                // &mut xyz,
                &mut trajectory,
                &mut trajectory_last_frames,
            );
            temp_energy_section = self.save_sections(
                &iiter,
                temp_energy_section,
                &mut temp_cn_dict_section,
                &mut amount_unique_levels,
            );
        }
        self.return_results(start_struct, lowest_energy_struct, seed)
        // Results {
        //     start,
        //     lowest_energy_struct,
        //     number_all_atoms: self.number_all_atoms,
        //     energy_section_list: self.energy_sections_list.clone(),
        //     cn_dict_sections: self.cn_dict_sections.clone(),
        //     seed,
        //     unique_levels: self.unique_levels.clone(),
        // }
    }

    pub fn write_exp_file(&self, exp: &Results) {
        let mut file = File::create(self.save_folder.clone() + "/exp_file.json").unwrap();
        file.write_all(serde_json::to_string_pretty(exp).unwrap().as_bytes())
            .unwrap();
    }
}
