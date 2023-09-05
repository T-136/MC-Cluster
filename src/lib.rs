use anyhow;
use chemfiles::{Atom, Frame, Trajectory, UnitCell};
use rand;
use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use rand::rngs::SmallRng;
use std::collections::hash_map::Entry;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::{cmp, eprint, fs, println};
use vasp_poscar::Poscar;

// mod build_pairs;
mod listdict;
mod read_files;
mod setup;
mod sim;

pub use sim::Results;

const CN: usize = 12;
const NN_PAIR_NUMBER: usize = 20;
const AMOUNT_SECTIONS: usize = 10000;
const SAVE_TH: u64 = 1000;

const GRID_SIZE: [u32; 3] = [15, 15, 15];

const VARIANT_A: bool = true;

#[derive(Clone)]
pub struct Simulation {
    niter: u64,
    number_all_atoms: u32,
    occ: Vec<u8>,
    onlyocc: HashSet<u32, fnv::FnvBuildHasher>,
    cn_metal: Vec<usize>,
    possible_moves: listdict::ListDict,
    total_energy_1000: i64,
    nn: HashMap<u32, [u32; CN], fnv::FnvBuildHasher>,
    // nn_pair: HashMap<u64, [u32; NN_PAIR_NUMBER], fnv::FnvBuildHasher>,
    // nnn_pair: HashMap<u64, [u32; 74], fnv::FnvBuildHasher>,
    // mut xyz: Vec<[f64; 3]>,
    xsites_positions: Vec<[f64; 3]>,
    unit_cell: UnitCell,
    cn_dict: [u32; CN + 1],
    save_folder: String,
    last_traj_frequency: u64,
    last_frames_trajectory_amount: Option<u64>,
    start_temperature: Option<f64>,
    temperature: f64,
    cn_dict_sections: Vec<HashMap<u8, f64>>,
    energy_sections_list: Vec<f64>,
    optimization_cut_off_fraction: Vec<u64>,
    unique_levels: HashMap<BTreeMap<u8, u32>, (i64, u64)>,
}

impl Simulation {
    pub fn new(
        niter: u64,
        input_file: Option<String>,
        atoms_input: Option<u32>,
        temperature: f64,
        start_temperature: Option<f64>,
        save_folder_name: String,
        pairlist_file: String,
        nn_pairlist_file: String,
        // nnn_pairlist_file: String,
        atom_sites: String,
        last_traj_frequency: u64,
        last_frames_trajectory_amount: Option<u64>,
        bulk_file_name: String,
        repetition: usize,
        optimization_cut_off_fraction: Vec<u64>,
    ) -> Simulation {
        let nsites: u32 = GRID_SIZE[0] * GRID_SIZE[1] * GRID_SIZE[2] * 4;
        let nn = read_files::read_nn(&pairlist_file);
        // let nn_pair = read_files::read_nn_pairlists(&nn_pairlist_file);
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
            unit_cell_size[0][0] * GRID_SIZE[0] as f64,
            unit_cell_size[1][1] * GRID_SIZE[1] as f64,
            unit_cell_size[2][2] * GRID_SIZE[2] as f64,
        ]);
        let mut cn_dict: [u32; CN + 1] = [0; CN + 1];

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
        let mut cn_metal: Vec<usize> = Vec::with_capacity(nsites as usize);
        for o in 0..nsites {
            let mut neighbors: u8 = 0;
            for o1 in nn[&o].iter() {
                if occ[*o1 as usize] == 1 {
                    // cn.entry(o).and_modify(|x| *x += 1).or_insert(1);
                    neighbors += 1;
                }
            }
            cn_metal.push(neighbors as usize);
            if occ[o as usize] == 1 {
                cn_dict[cn_metal[o as usize] as usize] += 1;
            };
        }
        // let mut former_energy_dict: fnv::FnvHashMap<u32, i64> =
        //     fnv::FnvHashMap::with_capacity_and_hasher(nsites as usize, Default::default());
        let mut total_energy_1000: i64 = 0;
        let mut possible_moves: listdict::ListDict = listdict::ListDict::new();
        for o in onlyocc.iter() {
            let energy_1000: i64 = sim::energy_calculation(o, &cn_metal);
            total_energy_1000 += energy_1000;
            // former_energy_dict.insert(o.clone(), energy_1000);

            for u in &nn[o] {
                if occ[*u as usize] == 0 {
                    // >1 so that atoms cant leave the cluster
                    // <x cant move if all neighbors are occupied
                    if cn_metal[*o as usize] < CN && cn_metal[*u as usize] > 1 {
                        possible_moves.add_item(o.clone(), u.clone())
                    }
                }
            }
        }

        let simulation_folder_name = match start_temperature {
            Some(start_temp) => {
                std::format!(
                    "{}-{}K_{}I_{}A",
                    start_temp,
                    temperature,
                    niter,
                    onlyocc.len()
                )
            }
            None => std::format!("{}K_{}I_{}A", temperature, niter, onlyocc.len()),
        };

        // let mut sub_folder = start_temperature
        //     .map_or(
        //         std::format!("{}K_{}I_{}A", temperature, niter, onlyocc.len()),
        //         |start_temp| {
        //             std::format!(
        //                 "{}-{}K_{}I_{}A",
        //                 start_temp,
        //                 temperature,
        //                 niter,
        //                 onlyocc.len()
        //             )
        //         },
        //     )
        //     .map(&save_folder_name);

        let mut sub_folder = save_folder_name + &simulation_folder_name;

        sub_folder = std::format!("{}_{}", sub_folder, repetition);

        fs::create_dir(&sub_folder).unwrap_or_else(|error| {
            println!(
                "could not create folder: {} with error: {}",
                sub_folder, error
            )
        });

        let cn_dict_sections = Vec::with_capacity(AMOUNT_SECTIONS);
        let energy_sections_list = Vec::with_capacity(AMOUNT_SECTIONS);
        let unique_levels = HashMap::new();

        Simulation {
            niter,
            number_all_atoms,
            occ,
            onlyocc,
            cn_metal,
            possible_moves,
            total_energy_1000,
            nn,
            // nn_pair,
            // nnn_pair,
            xsites_positions,
            unit_cell,
            cn_dict,
            save_folder: sub_folder,
            last_traj_frequency,
            last_frames_trajectory_amount,
            start_temperature,
            temperature,
            cn_dict_sections,
            energy_sections_list,
            optimization_cut_off_fraction,
            unique_levels,
        }
    }

    pub fn run(&mut self, mut amount_unique_levels: i32) -> Results {
        let mut rng_choose = SmallRng::from_entropy();
        // let choose_seed: [u8; 32] = rng_choose.get_seed();
        //

        // let mut rng_e_number = SmallRng::from_entropy();
        // let e_number_seed: [u8; 32] = rng_e_number.get_seed();
        let cut_off_perc = self.optimization_cut_off_fraction[0] as f64
            / self.optimization_cut_off_fraction[1] as f64;

        let seed = sim::Seed {
            rust: "used rust".to_string(),
            choose_seed: [0; 32],
            e_number_seed: [0; 32],
        };

        let mut trajectory_last_frames: Option<Trajectory> =
            self.last_frames_trajectory_amount.clone().map(|i| {
                let range = i * self.last_traj_frequency;
                Trajectory::open(
                    self.save_folder.clone() + &format!("/last_{range}_frames.xyz"),
                    'w',
                )
                .unwrap()
            });

        // if let Some(i) = self.last_frames_trajectory_amount {
        //     let range = i * self.last_traj_frequency;
        //     Some(
        //         Trajectory::open(
        //             self.save_folder.clone() + &format!("/last_{range}_frames.xyz"),
        //             'w',
        //         )
        //         .unwrap(),
        //     )
        // } else {
        //     None
        // };

        let mut lowest_energy_struct: sim::LowestEnergy = sim::LowestEnergy {
            energy: f64::INFINITY,
            cn_total: HashMap::new(),
            empty_cn: HashMap::new(),
            iiter: 0,
        };
        let mut temp_energy_section: i64 = 0;
        let mut temp_cn_dict_section: [u64; CN + 1] = [0; CN + 1];

        let start_energy = self.total_energy_1000 as f64 / 1000.;

        let mut start_cn_dict = HashMap::new();

        for (k, v) in self.cn_dict.iter().enumerate() {
            start_cn_dict.insert(k as u8, *v);
        }

        let start: sim::Start = sim::Start {
            start_energy: start_energy,
            start_cn: start_cn_dict,
        };

        if self.niter == 0 {
            self.save_lowest_energy(
                &0,
                &mut lowest_energy_struct,
                // &mut trajectory_lowest_energy,
            );
        }
        let section_size: u64 = self.niter / AMOUNT_SECTIONS as u64;
        println!("section_size: {}", section_size);
        println!("SAVE_TH: {}", SAVE_TH);
        println!("niter: {}", self.niter);

        for iiter in 0..self.niter {
            if iiter % section_size == 0 {
                println!(
                    "iteration {}; {}%",
                    iiter,
                    (iiter as f64 / self.niter as f64 * 100.)
                );
                // println!("{:?}", self.possible_moves.len());
            }
            let (move_from, move_to) = self.possible_moves.choose_random_item(&mut rng_choose);

            let energy1000_diff = self.energy_diff(move_from, move_to);

            if VARIANT_A
                && iiter * self.optimization_cut_off_fraction[1]
                    == self.niter * self.optimization_cut_off_fraction[0]
            {
                for o in 0..self.cn_metal.len() {
                    let mut neighbors: u8 = 0;
                    for o1 in self.nn[&(o as u32)].iter() {
                        if self.occ[*o1 as usize] == 1 {
                            // cn.entry(o).and_modify(|x| *x += 1).or_insert(1);
                            neighbors += 1;
                        }
                    }
                    self.cn_metal.push(neighbors as usize);
                    if self.occ[o as usize] == 1 {
                        self.cn_dict[self.cn_metal[o as usize] as usize] += 1;
                    };
                }
            }

            if self.is_acceptance_criteria_fulfilled(
                energy1000_diff,
                &mut rng_choose,
                iiter,
                cut_off_perc,
            ) {
                self.perform_move(move_from, move_to, energy1000_diff, iiter);
                self.update_possible_moves(move_from, move_to)
            }

            if iiter * self.optimization_cut_off_fraction[1]
                >= self.niter * self.optimization_cut_off_fraction[0]
            {
                self.save_lowest_energy(
                    &iiter,
                    &mut lowest_energy_struct,
                    // &mut trajectory_lowest_energy,
                )
            }

            self.cond_last_frams_traj_write(
                &iiter,
                // &mut xyz,
                // &mut trajectory,
                &mut trajectory_last_frames,
            );
            temp_energy_section = self.save_sections(
                &iiter,
                temp_energy_section,
                &mut temp_cn_dict_section,
                &mut amount_unique_levels,
                section_size,
            );
        }

        Results {
            start,
            lowest_energy_struct,
            number_all_atoms: self.number_all_atoms,
            energy_section_list: self.energy_sections_list.clone(),
            cn_dict_sections: self.cn_dict_sections.clone(),
            seed,
            unique_levels: self.unique_levels.clone(),
        }
    }

    pub fn write_exp_file(&self, exp: &Results) {
        let mut file = File::create(self.save_folder.clone() + "/exp_file.json").unwrap();
        file.write_all(serde_json::to_string_pretty(exp).unwrap().as_bytes())
            .unwrap();
    }

    fn save_sections(
        &mut self,
        iiter: &u64,
        mut temp_energy_section_1000: i64,
        temp_cn_dict_section: &mut [u64; CN + 1],
        amount_unique_levels: &mut i32,
        section_size: u64,
    ) -> i64 {
        if (iiter + 1) % SAVE_TH == 0 {
            temp_energy_section_1000 += self.total_energy_1000;

            temp_cn_dict_section
                .iter_mut()
                .enumerate()
                .for_each(|(i, v)| *v += self.cn_dict[i] as u64);
        }

        if *amount_unique_levels != 0 {
            let mut cn_hash_map = HashMap::new();
            for (i, v) in self.cn_dict.into_iter().enumerate() {
                cn_hash_map.insert(i as u8, v);
            }

            if *iiter * self.optimization_cut_off_fraction[0]
                >= self.niter * self.optimization_cut_off_fraction[0]
            {
                let cn_btree: BTreeMap<_, _> = cn_hash_map.into_iter().collect();
                match self.unique_levels.entry(cn_btree) {
                    Entry::Occupied(mut entry) => {
                        let (_, x) = entry.get_mut();
                        *x += 1;
                    }
                    Entry::Vacant(entry) => {
                        if *amount_unique_levels == 1 {
                            eprint!("amount_unique_levels reached");
                        }
                        *amount_unique_levels -= 1;
                        entry.insert((self.total_energy_1000 / 1000, 1));
                    }
                }
                // if *amount_unique_levels == 1 {
                //     eprint!("amount_unique_levels reached");
                // }
                // self.unique_levels
                //     .entry(cn_btree)
                //     .and_modify(|(_, x)| *x += 1)
                //     .or_insert_with(|| {
                //         *amount_unique_levels -= 1;
                //         (self.total_energy_1000 / 1000, 1)
                //     });
            }
        }
        if (iiter + 1) % section_size == 0 {
            self.energy_sections_list
                .push(temp_energy_section_1000 as f64 / (section_size / SAVE_TH) as f64 / 1000.);
            temp_energy_section_1000 = 0;
            // temp_energy_section.clear();

            let mut section: HashMap<u8, f64> = HashMap::new();
            for (k, list) in temp_cn_dict_section.iter_mut().enumerate() {
                section.insert(k as u8, *list as f64 / section_size as f64);
                // list.clear();
                *list = 0;
            }
            assert_eq!(temp_cn_dict_section, &mut [0_u64; CN + 1]);
            self.cn_dict_sections.push(section.clone())
        }
        temp_energy_section_1000
    }

    fn save_lowest_energy(&mut self, iiter: &u64, lowest_energy_struct: &mut sim::LowestEnergy) {
        if &lowest_energy_struct.energy > &(self.total_energy_1000 as f64 / 1000.) {
            let mut empty_neighbor_cn: HashMap<u8, u32> = HashMap::new();
            let empty_set: HashSet<&u32> =
                HashSet::from_iter(self.possible_moves.iter().map(|(_, empty)| empty));
            for empty in empty_set {
                if self.cn_metal[*empty as usize] > 3 {
                    empty_neighbor_cn
                        .entry(self.cn_metal[*empty as usize] as u8)
                        .and_modify(|x| *x += 1)
                        .or_insert(1);
                }
            }
            lowest_energy_struct.empty_cn = empty_neighbor_cn;
            lowest_energy_struct.energy = self.total_energy_1000.clone() as f64 / 1000.;
            lowest_energy_struct.iiter = iiter.clone();

            let mut cn_hash_map: HashMap<u8, u32> = HashMap::new();
            for (i, v) in self.cn_dict.into_iter().enumerate() {
                cn_hash_map.insert(i as u8, v);
            }
            lowest_energy_struct.cn_total = cn_hash_map;

            let mut trajectory_lowest_energy =
                Trajectory::open(self.save_folder.clone() + "/lowest_energy.xyz", 'w').unwrap();

            self.write_traj(&mut trajectory_lowest_energy);
        };
    }

    fn cond_last_frams_traj_write(
        &self,
        iiter: &u64,
        trajectory_last_frames_option: &mut Option<Trajectory>,
    ) {
        if let Some(traj_last_frames) = trajectory_last_frames_option {
            if (self.niter - iiter
                <= self.last_frames_trajectory_amount.unwrap() * self.last_traj_frequency)
                && ((self.niter - iiter) % self.last_traj_frequency == 0)
            {
                self.write_traj(traj_last_frames);
            }
        }
    }

    fn write_traj(&self, trajectory: &mut Trajectory) {
        let mut xyz: Vec<[f64; 3]> = Vec::new();
        for (j, ii) in self.onlyocc.iter().enumerate() {
            xyz.insert(j, self.xsites_positions[ii.clone() as usize]);
        }
        let mut frame = Frame::new();
        frame.set_cell(&self.unit_cell);

        for atom in xyz.into_iter() {
            frame.add_atom(&Atom::new("Pt"), [atom[0], atom[1], atom[2]], None);
        }

        trajectory
            .write(&mut frame)
            .unwrap_or_else(|x| eprintln!("{}", x));
    }

    fn calculate_current_temp(&self, iiter: u64, cut_off_perc: f64) -> f64 {
        let heating_temp = 5500.;
        if self.start_temperature.is_some() {
            if (iiter + 1) as f64 <= self.niter as f64 * cut_off_perc {
                heating_temp
                    - ((iiter + 1) as f64 / (self.niter as f64 * cut_off_perc))
                        * (heating_temp - self.start_temperature.unwrap())
            } else {
                self.start_temperature.unwrap()
                    - ((iiter + 1) as f64 - self.niter as f64 * cut_off_perc)
                        / (self.niter as f64 * (1. - cut_off_perc))
                        * (self.start_temperature.unwrap() - self.temperature)
            }
        } else {
            if (iiter + 1) as f64 <= self.niter as f64 * cut_off_perc {
                heating_temp
                    - (iiter as f64 / (self.niter as f64 * cut_off_perc))
                        * (heating_temp - self.temperature)
            } else {
                self.temperature
            }
        }
    }

    fn is_acceptance_criteria_fulfilled(
        &mut self,
        energy1000_diff: i64,
        rng_e_number: &mut SmallRng,
        iiter: u64,
        cut_off_perc: f64,
    ) -> bool {
        const KB: f64 = 8.6173324e-5;
        // if self.start_temperature.is_some() {
        if energy1000_diff < 0 {
            return true;
        }
        let acceptance_temp = self.calculate_current_temp(iiter, cut_off_perc);
        let between = Uniform::new_inclusive(0., 1.);
        // let delta_energy = proposed_energy - self.total_energy_1000
        let rand_value = between.sample(rng_e_number);
        (rand_value) < ((-energy1000_diff as f64 / 1000.) / (KB * acceptance_temp)).exp()
    }

    pub fn energy_diff(&self, move_from: u32, move_to: u32) -> i64 {
        const M_BETA: i64 = -0330;
        const M_ALPHA: i64 = 3960;
        (2 * ((self.cn_metal[move_to as usize] as i64 - 1) * M_BETA + M_ALPHA))
            - (2 * ((self.cn_metal[move_from as usize] as i64) * M_BETA + M_ALPHA))
    }

    fn perform_move(&mut self, move_from: u32, move_to: u32, energy1000_diff: i64, iiter: u64) {
        self.occ[move_to as usize] = self.occ[move_from as usize]; // covers different alloys also
        self.occ[move_from as usize] = 0;

        self.onlyocc.remove(&move_from);
        self.onlyocc.insert(move_to);

        // let nn_intersection: Vec<u32> = self.nn[&move_from]
        //     .into_iter()
        //     .filter(|x| self.nn[&move_to].contains(x))
        //     .collect();

        self.cn_dict[self.cn_metal[move_from as usize]] -= 1;
        for o in self.nn[&move_from] {
            // if !nn_intersection.contains(&o) {
            if VARIANT_A
                && iiter * self.optimization_cut_off_fraction[1]
                    >= self.niter * self.optimization_cut_off_fraction[0]
            {
                if self.occ[o as usize] == 1 && o != move_to {
                    self.cn_dict[self.cn_metal[o as usize]] -= 1;
                    self.cn_dict[self.cn_metal[o as usize] - 1] += 1;
                }
            }
            self.cn_metal[o as usize] -= 1;
            // }
        }
        for o in self.nn[&move_to] {
            // if !nn_intersection.contains(&o) {
            if VARIANT_A
                && iiter * self.optimization_cut_off_fraction[1]
                    >= self.niter * self.optimization_cut_off_fraction[0]
            {
                if self.occ[o as usize] == 1 && o != move_from {
                    self.cn_dict[self.cn_metal[o as usize]] -= 1;
                    self.cn_dict[self.cn_metal[o as usize] + 1] += 1;
                }
            }
            self.cn_metal[o as usize] += 1;
            // }
        }
        self.cn_dict[self.cn_metal[move_to as usize]] += 1;

        self.total_energy_1000 += energy1000_diff;
    }

    fn update_possible_moves(&mut self, move_from: u32, move_to: u32) {
        self.possible_moves.remove_item(move_from, move_to);
        for neighbor_atom in self.nn[&move_from] {
            if self.occ[neighbor_atom as usize] == 0 {
                self.possible_moves.remove_item(move_from, neighbor_atom);
            }
            if self.occ[neighbor_atom as usize] == 1 {
                // greater than one because of neighbor moving in this spot
                if self.cn_metal[move_from as usize] > 1 {
                    self.possible_moves.add_item(neighbor_atom, move_from);
                }
            }
        }

        for empty_neighbor in self.nn[&move_to] {
            if self.occ[empty_neighbor as usize] == 1 {
                self.possible_moves.remove_item(empty_neighbor, move_to);
            }
            if self.occ[empty_neighbor as usize] == 0 {
                // greater than one because of neighbor moving in this spot
                if self.cn_metal[empty_neighbor as usize] > 1 {
                    self.possible_moves.add_item(move_to, empty_neighbor);
                }
            }
        }
    }
}

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
