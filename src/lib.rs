use anyhow;
use chemfiles::{Atom, Frame, Trajectory, UnitCell};
use core::panic;
use csv::Writer;
use energy::EnergyInput;
use rand;
use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use rand::rngs::SmallRng;
use rayon::prelude::*;
use std::collections::hash_map::Entry;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::sync::Arc;
use std::{cmp, eprint, fs, println, usize};

pub mod energy;
mod grid_structure;
mod listdict;
mod read_and_write;
mod setup;
mod sim;

pub use grid_structure::GridStructure;
pub use sim::Results;

const CN: usize = 12;
// const GCN: usize = 54;
const GCN: usize = 145;
const NN_PAIR_NUMBER: usize = 20;
const NN_PAIR_NO_INTERSEC_NUMBER: usize = 7;
const NNN_PAIR_NO_INTERSEC_NUMBER: usize = 20;
const AMOUNT_SECTIONS: usize = 10000;
const SAVE_TH: u64 = 1000;

const GRID_SIZE: [u32; 3] = [9, 9, 9];

const SAVE_ENTIRE_SIM: bool = false;

#[derive(Clone)]
pub struct Simulation {
    niter: u64,
    number_all_atoms: u32,
    occ: Vec<u8>,
    onlyocc: HashSet<u32, fnv::FnvBuildHasher>,
    cn_metal: Vec<usize>,
    gcn_metal: Vec<usize>,
    nn_support: Option<Vec<u8>>,
    possible_moves: listdict::ListDict,
    total_energy_1000: i64,
    cn_dict: [u32; CN + 1],
    save_folder: String,
    start_temperature: Option<f64>,
    temperature: f64,
    cn_dict_sections: Vec<HashMap<u8, f64>>,
    energy_sections_list: Vec<f64>,
    optimization_cut_off_fraction: Vec<u64>,
    unique_levels: HashMap<BTreeMap<u8, u32>, (i64, u64)>,
    heat_map: Option<Vec<u64>>,
    snap_shot_sections: Option<Vec<Vec<u8>>>,
    heat_map_sections: Vec<Vec<u64>>,
    energy: EnergyInput,
    gridstructure: &'static GridStructure,
    support_e: i64,
}

impl Simulation {
    pub fn new(
        niter: u64,
        input_file: Option<String>,
        atoms_input: Option<u32>,
        temperature: f64,
        start_temperature: Option<f64>,
        save_folder_name: String,
        write_snap_shots: bool,
        is_heat_map: bool,
        repetition: usize,
        optimization_cut_off_fraction: Vec<u64>,
        energy: EnergyInput,
        support_indices: Option<Vec<u32>>,
        gridstructure: &'static GridStructure,
        support_e: i64,
    ) -> Simulation {
        // 12
        let nsites: u32 = GRID_SIZE[0] * GRID_SIZE[1] * GRID_SIZE[2] * 4;
        let mut cn_dict: [u32; CN + 1] = [0; CN + 1];
        let (occ, onlyocc, number_all_atoms, nn_support) = if input_file.is_some() {
            let xyz = read_and_write::read_sample(&input_file.unwrap());
            let (occ, onlyocc) =
                setup::occ_onlyocc_from_xyz(&xyz, nsites, &gridstructure.xsites_positions);
            let number_of_atoms: u32 = onlyocc.len() as u32;
            (occ, onlyocc, number_of_atoms, None)
        } else if atoms_input.is_some() {
            let number_of_atom = atoms_input.unwrap();
            let (occ, onlyocc, nn_support) = setup::create_input_cluster(
                &atoms_input.unwrap(),
                &gridstructure.xsites_positions,
                &gridstructure.nn,
                nsites,
                support_indices,
            );
            (occ, onlyocc, number_of_atom, nn_support)
        } else {
            panic!("gib input atoms or input file");
        };
        let mut cn_metal: Vec<usize> = Vec::with_capacity(nsites as usize);

        for o in 0..nsites {
            let mut neighbors: u8 = 0;
            for o1 in gridstructure.nn[&o].iter() {
                if occ[*o1 as usize] == 1 {
                    // cn.entry(o).and_modify(|x| *x += 1).or_insert(1);
                    neighbors += 1;
                }
            }
            cn_metal.push(neighbors as usize);
            if occ[o as usize] == 1 {
                cn_dict[cn_metal[o as usize]] += 1;
            };
        }
        let mut gcn_metal: Vec<usize> = Vec::with_capacity(nsites as usize);
        for o in 0..nsites {
            let mut gcn: usize = 0;
            for o1 in gridstructure.nn[&o].iter() {
                if occ[*o1 as usize] == 1 {
                    gcn += cn_metal[*o1 as usize];
                }
            }
            gcn_metal.push(gcn);
        }

        let mut total_energy_1000: i64 = 0;
        let mut possible_moves: listdict::ListDict = listdict::ListDict::new(GRID_SIZE);
        for o in onlyocc.iter() {
            let at_support = if let Some(nn_support) = &nn_support {
                nn_support[*o as usize]
            } else {
                0
            };
            match energy {
                EnergyInput::LinearCn(_) | EnergyInput::Cn(_) => {
                    total_energy_1000 += energy::energy_1000_calculation(
                        &energy,
                        cn_metal[*o as usize],
                        at_support,
                        support_e,
                    );
                }

                EnergyInput::LinearGcn(_) | EnergyInput::Gcn(_) => {
                    total_energy_1000 += energy::energy_1000_calculation(
                        &energy,
                        gcn_metal[*o as usize],
                        at_support,
                        support_e,
                    );
                }
            };
            // total_energy_1000 += energy::energy_1000_calculation(
            //     &energy,
            //     gcn_metal[*o as usize],
            //     at_support,
            //     support_e,
            // );

            for u in &gridstructure.nn[o] {
                if occ[*u as usize] == 0 {
                    // >1 so that atoms cant leave the cluster
                    // <x cant move if all neighbors are occupied
                    if cn_metal[*u as usize] > 1 {
                        possible_moves.add_item(*o, *u, None)
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

        let snap_shot_sections: Option<Vec<Vec<u8>>> = if write_snap_shots {
            Some(Vec::new())
        } else {
            None
        };

        let heat_map: Option<Vec<u64>> = if is_heat_map {
            Some(vec![0; nsites as usize])
        } else {
            None
        };

        let heat_map_sections: Vec<Vec<u64>> = Vec::new();

        Simulation {
            niter,
            number_all_atoms,
            occ,
            onlyocc,
            cn_metal,
            gcn_metal,
            nn_support,
            possible_moves,
            total_energy_1000,
            cn_dict,
            save_folder: sub_folder,
            start_temperature,
            temperature,
            cn_dict_sections,
            energy_sections_list,
            optimization_cut_off_fraction,
            unique_levels,
            snap_shot_sections,
            heat_map,
            heat_map_sections,
            energy,
            gridstructure,
            support_e,
        }
    }

    pub fn run(&mut self, mut amount_unique_levels: i32) -> Results {
        let mut rng_choose = SmallRng::from_entropy();

        let cut_off_perc = self.optimization_cut_off_fraction[0] as f64
            / self.optimization_cut_off_fraction[1] as f64;

        let mut lowest_energy_struct = sim::LowestEnergy::default();

        let mut temp_energy_section: i64 = 0;
        let mut temp_cn_dict_section: [u64; CN + 1] = [0; CN + 1];

        let start = sim::Start::new(self.total_energy_1000, &self.cn_dict);

        let mut lowest_e_onlyocc: HashSet<u32, fnv::FnvBuildHasher> =
            fnv::FnvHashSet::with_capacity_and_hasher(
                self.number_all_atoms as usize,
                Default::default(),
            );
        if self.niter == 0 {
            if let Some(x) = self.opt_save_lowest_energy(&0, &mut lowest_energy_struct) {
                lowest_e_onlyocc = x;
            };
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
                // println!("{:?}", self.cn_metal);
            }
            let is_recording_sections = iiter * self.optimization_cut_off_fraction[1]
                >= self.niter * self.optimization_cut_off_fraction[0];

            if !SAVE_ENTIRE_SIM
                && iiter * self.optimization_cut_off_fraction[1]
                    == self.niter * self.optimization_cut_off_fraction[0]
            {
                self.cn_dict.iter_mut().for_each(|x| {
                    *x = 0;
                });
                for o in 0..self.cn_metal.len() {
                    if self.occ[o] == 1 {
                        self.cn_dict[self.cn_metal[o]] += 1;
                    };
                }
            };

            let (move_from, move_to, _) =
                self.possible_moves.choose_random_item_mc(&mut rng_choose);

            let energy1000_diff = self.energy_change_by_move(move_from, move_to);

            if self.is_acceptance_criteria_fulfilled(
                energy1000_diff,
                &mut rng_choose,
                iiter,
                cut_off_perc,
            ) {
                self.perform_move(move_from, move_to, energy1000_diff, is_recording_sections);
                self.update_possible_moves(move_from, move_to);
                if let Some(map) = &mut self.heat_map {
                    map[move_to as usize] += 1;
                    map[move_from as usize] += 1;
                }
            }

            self.cond_snap_and_heat_map(&iiter);

            if iiter * self.optimization_cut_off_fraction[1]
                >= self.niter * self.optimization_cut_off_fraction[0]
            {
                if let Some(x) = self.opt_save_lowest_energy(&iiter, &mut lowest_energy_struct) {
                    lowest_e_onlyocc = x;
                };
            }

            if SAVE_ENTIRE_SIM || is_recording_sections {
                temp_energy_section = self.save_sections(
                    &iiter,
                    temp_energy_section,
                    &mut temp_cn_dict_section,
                    &mut amount_unique_levels,
                    section_size,
                );
            }
        }

        println!("heatmap section len: {:?}", self.heat_map_sections.len());

        read_and_write::write_occ_as_xyz(
            self.save_folder.clone(),
            lowest_e_onlyocc,
            &self.gridstructure.xsites_positions,
            &self.gridstructure.unit_cell,
            &self.occ,
        );

        if self.heat_map.is_some() {
            let mut wtr = Writer::from_path(self.save_folder.clone() + "/heat_map.csv").unwrap();
            for heat_section in &self.heat_map_sections {
                wtr.write_record(heat_section.iter().map(|x| x.to_string()))
                    .unwrap();
            }
            wtr.flush().unwrap();
        }

        if self.snap_shot_sections.is_some() {
            let mut wtr =
                Writer::from_path(self.save_folder.clone() + "/snap_shot_sections.csv").unwrap();
            if let Some(snap_shot_sections) = &self.snap_shot_sections {
                for heat_section in snap_shot_sections {
                    wtr.write_record(heat_section.iter().map(|x| x.to_string()))
                        .unwrap();
                }
            }
            wtr.flush().unwrap();
        }

        Results {
            start,
            lowest_energy_struct,
            number_all_atoms: self.number_all_atoms,
            energy_section_list: self.energy_sections_list.clone(),
            cn_dict_sections: self.cn_dict_sections.clone(),
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

            if *iiter * self.optimization_cut_off_fraction[1]
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
            }
        }
        if (iiter + 1) % section_size == 0 {
            self.energy_sections_list
                .push(temp_energy_section_1000 as f64 / (section_size / SAVE_TH) as f64 / 1000.);
            temp_energy_section_1000 = 0;

            let mut section: HashMap<u8, f64> = HashMap::new();
            for (k, list) in temp_cn_dict_section.iter_mut().enumerate() {
                section.insert(k as u8, *list as f64 / (section_size / SAVE_TH) as f64);
                *list = 0;
            }
            assert_eq!(temp_cn_dict_section, &mut [0_u64; CN + 1]);
            self.cn_dict_sections.push(section.clone())
        }
        temp_energy_section_1000
    }

    fn opt_save_lowest_energy(
        &mut self,
        iiter: &u64,
        lowest_energy_struct: &mut sim::LowestEnergy,
    ) -> Option<HashSet<u32, fnv::FnvBuildHasher>> {
        if lowest_energy_struct.energy > (self.total_energy_1000 as f64 / 1000.) {
            let mut empty_neighbor_cn: HashMap<u8, u32> = HashMap::new();
            let empty_set: HashSet<&u32> =
                HashSet::from_iter(self.possible_moves.iter().map(|(_, empty, _)| empty));
            for empty in empty_set {
                if self.cn_metal[*empty as usize] > 3 {
                    empty_neighbor_cn
                        .entry(self.cn_metal[*empty as usize] as u8)
                        .and_modify(|x| *x += 1)
                        .or_insert(1);
                }
            }
            lowest_energy_struct.empty_cn = empty_neighbor_cn;
            lowest_energy_struct.energy = self.total_energy_1000 as f64 / 1000.;
            lowest_energy_struct.iiter = *iiter;

            let mut cn_hash_map: HashMap<u8, u32> = HashMap::new();
            for (i, v) in self.cn_dict.into_iter().enumerate() {
                cn_hash_map.insert(i as u8, v);
            }
            lowest_energy_struct.cn_total = cn_hash_map;

            Some(self.onlyocc.clone())
        } else {
            None
        }
    }

    fn calculate_current_temp(&self, iiter: u64, cut_off_perc: f64) -> f64 {
        let heating_temp = 5000.;
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
        if energy1000_diff <= 0 {
            return true;
        }
        let acceptance_temp = self.calculate_current_temp(iiter, cut_off_perc);
        let between = Uniform::new_inclusive(0., 1.);
        let rand_value = between.sample(rng_e_number);
        (rand_value) < ((-energy1000_diff as f64 / 1000.) / (KB * acceptance_temp)).exp()
    }

    fn perform_move(
        &mut self,
        move_from: u32,
        move_to: u32,
        energy1000_diff: i64,
        is_recording_sections: bool,
    ) {
        self.occ[move_to as usize] = self.occ[move_from as usize]; // covers different alloys also
        self.occ[move_from as usize] = 0;

        self.onlyocc.remove(&move_from);
        self.onlyocc.insert(move_to);

        if SAVE_ENTIRE_SIM || is_recording_sections {
            self.cn_dict[self.cn_metal[move_from as usize]] -= 1;
        }
        // let (from_change, to_change) = self.no_int_from_move(move_from, move_to);
        for o in self.gridstructure.nn[&move_from] {
            if (SAVE_ENTIRE_SIM || is_recording_sections)
                && self.occ[o as usize] == 1
                && o != move_to
            {
                self.cn_dict[self.cn_metal[o as usize]] -= 1;
                self.cn_dict[self.cn_metal[o as usize] - 1] += 1;
            }
            self.cn_metal[o as usize] -= 1;
        }
        for o in self.gridstructure.nn[&move_to] {
            if (SAVE_ENTIRE_SIM || is_recording_sections)
                && self.occ[o as usize] == 1
                && o != move_from
            {
                self.cn_dict[self.cn_metal[o as usize]] -= 1;
                self.cn_dict[self.cn_metal[o as usize] + 1] += 1;
            }
            self.cn_metal[o as usize] += 1;
        }

        match self.energy {
            EnergyInput::LinearCn(_) | EnergyInput::Cn(_) => {}
            EnergyInput::LinearGcn(_) | EnergyInput::Gcn(_) => {
                let (from_change, to_change, intersect, is_reverse) = no_int_nnn_from_move(
                    move_from,
                    move_to,
                    &self.gridstructure.nnn_pair_no_intersec,
                );
                let (from_change_nn, to_change_nn) = no_int_nn_from_move(
                    move_from,
                    move_to,
                    &self.gridstructure.nn_pair_no_intersec,
                );
                for o in from_change_nn {
                    if o == move_to {
                        continue;
                    }
                    if self.occ[o as usize] == 1 {
                        self.gcn_metal[move_from as usize] -= 1;
                    }
                }
                self.gcn_metal[move_from as usize] += self.cn_metal[move_to as usize];
                for o in to_change_nn {
                    if o == move_from {
                        continue;
                    }
                    if self.occ[o as usize] == 1 {
                        self.gcn_metal[move_to as usize] += 1;
                    }
                }
                self.gcn_metal[move_to as usize] -= self.cn_metal[move_from as usize] - 1;
                for atom_and_neighbors in to_change {
                    for n in atom_and_neighbors.iter().skip(1) {
                        if n == &move_to {
                            self.gcn_metal[atom_and_neighbors[0] as usize] +=
                                self.cn_metal[*n as usize];
                            continue;
                        }
                        #[cfg(debug_assertions)]
                        if n == &move_from {
                            panic!(
                                "found move from {:?}, move_from {}",
                                atom_and_neighbors, move_from
                            );
                        }
                        if self.occ[*n as usize] == 1 {
                            self.gcn_metal[atom_and_neighbors[0] as usize] += 1;
                        }
                    }
                }
                for atom_and_neighbors in from_change {
                    // println!("gcn bef{:?}", self.gcn_metal[*atom as usize]);
                    for n in atom_and_neighbors.iter().skip(1) {
                        #[cfg(debug_assertions)]
                        if n == &move_to {
                            panic!("found move to");
                        }
                        if n == &move_from {
                            // println!("cn move from {:?}", self.cn_metal[*n as usize] - 1);
                            self.gcn_metal[atom_and_neighbors[0] as usize] -=
                                self.cn_metal[*n as usize] - 1;
                            continue;
                        }
                        if self.occ[*n as usize] == 1 {
                            // println!("-1",);
                            self.gcn_metal[atom_and_neighbors[0] as usize] -= 1;
                        }
                    }
                }
                for (atom, first_neighbors, second_neighbors, to_from_atoms) in intersect {
                    if !is_reverse {
                        for n in first_neighbors {
                            if self.occ[*n as usize] == 1 {
                                self.gcn_metal[*atom as usize] -= 1;
                            }
                        }
                        for n in second_neighbors {
                            if self.occ[*n as usize] == 1 {
                                self.gcn_metal[*atom as usize] += 1;
                            }
                        }
                    } else if is_reverse {
                        for n in second_neighbors {
                            if self.occ[*n as usize] == 1 {
                                self.gcn_metal[*atom as usize] -= 1;
                            }
                        }
                        for n in first_neighbors {
                            if self.occ[*n as usize] == 1 {
                                self.gcn_metal[*atom as usize] += 1;
                            }
                        }
                    }
                    for n in to_from_atoms {
                        if n == &move_to {
                            self.gcn_metal[*atom as usize] += self.cn_metal[*n as usize];
                            continue;
                        }
                        if n == &move_from {
                            self.gcn_metal[*atom as usize] -= self.cn_metal[*n as usize] - 1;
                            continue;
                        }
                        panic!("neither start nor end found");
                    }
                }
            }
        }
        if SAVE_ENTIRE_SIM || is_recording_sections {
            self.cn_dict[self.cn_metal[move_to as usize]] += 1;
        }

        self.total_energy_1000 += energy1000_diff;
    }

    fn energy_change_by_move(&self, move_from: u32, move_to: u32) -> i64 {
        let (from_at_support, to_at_support) = if let Some(nn_support) = &self.nn_support {
            let from_at_support = nn_support[move_from as usize];
            let to_at_support = nn_support[move_to as usize];
            (from_at_support, to_at_support)
        } else {
            (0, 0)
        };

        match self.energy {
            EnergyInput::LinearCn(energy_l_cn) => energy::energy_diff_l_cn(
                energy_l_cn,
                self.cn_metal[move_from as usize],
                self.cn_metal[move_to as usize] - 1,
                from_at_support,
                to_at_support,
                self.support_e,
            ),
            EnergyInput::Cn(energy_cn) => {
                let (from_change, to_change) = no_int_nn_from_move(
                    move_from,
                    move_to,
                    &self.gridstructure.nn_pair_no_intersec,
                );

                energy::energy_diff_cn(
                    energy_cn,
                    from_change
                        .iter()
                        .filter(|x| self.occ[**x as usize] == 1)
                        .map(|x| self.cn_metal[*x as usize]),
                    to_change
                        .iter()
                        .filter(|x| self.occ[**x as usize] == 1)
                        .map(|x| self.cn_metal[*x as usize]),
                    self.cn_metal[move_from as usize],
                    self.cn_metal[move_to as usize],
                    from_at_support,
                    to_at_support,
                    self.support_e,
                )
            }
            EnergyInput::LinearGcn(energy_l_gcn) => {
                let (from_change, to_change) = no_int_nn_from_move(
                    move_from,
                    move_to,
                    &self.gridstructure.nn_pair_no_intersec,
                );
                energy::energy_diff_l_gcn(
                    energy_l_gcn,
                    from_change
                        .iter()
                        .filter(|x| self.occ[**x as usize] == 1)
                        .map(|x| {
                            let mut gcn = 0;
                            for o in self.gridstructure.nn[x] {
                                if self.occ[o as usize] == 1 {
                                    gcn += self.cn_metal[o as usize]
                                }
                            }
                            gcn
                        }),
                    to_change
                        .iter()
                        .filter(|x| self.occ[**x as usize] == 1)
                        .map(|x| {
                            let mut gcn = 0;
                            for o in self.gridstructure.nn[x] {
                                if self.occ[o as usize] == 1 {
                                    gcn += self.cn_metal[o as usize]
                                }
                            }
                            gcn
                        }),
                    self.cn_metal[move_from as usize],
                    self.cn_metal[move_to as usize],
                )
            }
            EnergyInput::Gcn(energy_gcn) => {
                let (from_change, to_change, intersect, is_reverse) = no_int_nnn_from_move(
                    move_from,
                    move_to,
                    &self.gridstructure.nnn_pair_no_intersec,
                );

                let (from_change_nn, to_change_nn) = no_int_nn_from_move(
                    move_from,
                    move_to,
                    &self.gridstructure.nn_pair_no_intersec,
                );

                let mut cn_from = 0;
                to_change_nn
                    .iter()
                    .filter(|x| self.occ[**x as usize] == 1)
                    .for_each(|_| cn_from += 1);

                energy::energy_diff_gcn(
                    energy_gcn,
                    from_change
                        .iter()
                        .filter_map(|x| self.to_change_map(x, move_from, FromOrTo::From)),
                    to_change.iter().filter_map(|atom_and_neighbors| {
                        self.to_change_map(atom_and_neighbors, move_to, FromOrTo::To)
                    }),
                    intersect
                        .iter()
                        .filter(|atom_and_neighbors| self.occ[atom_and_neighbors.0 as usize] == 1)
                        .map(|atom_and_neighbors| {
                            self.map_intersec(atom_and_neighbors, move_to, move_from, is_reverse)
                        }),
                    self.gcn_metal[move_from as usize],
                    self.gcn_metal[move_to as usize] - self.cn_metal[move_from as usize] + cn_from,
                )
            }
        }
    }

    fn map_intersec(
        &self,
        atom_and_neighbors: &(u32, Vec<u32>, Vec<u32>, Vec<u32>),
        move_to: u32,
        move_from: u32,
        is_reverse: bool,
    ) -> (usize, usize) {
        let old_gcn = self.gcn_metal[atom_and_neighbors.0 as usize];
        let mut new_gcn = self.gcn_metal[atom_and_neighbors.0 as usize];
        let (atom, first_neighbors, second_neighbors, to_from_atoms) = atom_and_neighbors;
        first_neighbors
            .iter()
            .filter(|atom| {
                self.occ[**atom as usize] == 1 || **atom == move_to || **atom == move_from
            })
            .for_each(|_| {
                if !is_reverse {
                    new_gcn -= 1
                } else {
                    new_gcn += 1
                }
            });
        second_neighbors
            .iter()
            .filter(|atom| {
                self.occ[**atom as usize] == 1 || **atom == move_to || **atom == move_from
            })
            .for_each(|x| {
                if is_reverse {
                    new_gcn -= 1
                } else if !is_reverse {
                    new_gcn += 1
                }
            });
        to_from_atoms
            .iter()
            .filter(|atom| {
                self.occ[**atom as usize] == 1 || **atom == move_to || **atom == move_from
            })
            .for_each(|x| {
                if x == &move_to {
                    new_gcn += self.cn_metal[*x as usize] - 1
                } else if x == &move_from {
                    new_gcn -= self.cn_metal[*x as usize]
                }
            });
        (old_gcn, new_gcn)
    }

    fn to_change_map(
        &self,
        atom_and_neighbors: &Vec<u32>,
        move_from_or_to: u32,
        from_or_to: FromOrTo,
        // move_to: u32,
    ) -> Option<(usize, usize)> {
        if self.occ[*atom_and_neighbors.first().unwrap() as usize] == 1 {
            let old_gcn = self.gcn_metal[*atom_and_neighbors.first().unwrap() as usize];
            let mut new_gcn = self.gcn_metal[*atom_and_neighbors.first().unwrap() as usize];
            match from_or_to {
                FromOrTo::From => {
                    atom_and_neighbors
                        .iter()
                        .skip(1)
                        .filter(|x| self.occ[**x as usize] == 1 || **x == move_from_or_to)
                        .for_each(|x| {
                            if x == &move_from_or_to {
                                new_gcn -= self.cn_metal[*x as usize]
                            } else {
                                new_gcn -= 1
                            }
                        });
                }
                FromOrTo::To => {
                    atom_and_neighbors
                        .iter()
                        .skip(1)
                        .filter(|x| self.occ[**x as usize] == 1 || **x == move_from_or_to)
                        .for_each(|x| {
                            if x == &move_from_or_to {
                                new_gcn += self.cn_metal[*x as usize] - 1
                            } else {
                                new_gcn += 1
                            }
                        });
                }
            }
            Some((old_gcn, new_gcn))
        } else {
            None
        }
    }

    fn update_possible_moves(&mut self, move_from: u32, move_to: u32) {
        self.possible_moves.remove_item(move_from, move_to);
        for neighbor_atom in self.gridstructure.nn[&move_from] {
            if self.occ[neighbor_atom as usize] == 0 {
                self.possible_moves.remove_item(move_from, neighbor_atom);
            }
            if self.occ[neighbor_atom as usize] == 1 {
                // greater than one because of neighbor moving in this spot
                if self.cn_metal[move_from as usize] > 1 {
                    self.possible_moves.add_item(neighbor_atom, move_from, None);
                }
            }
        }

        for empty_neighbor in self.gridstructure.nn[&move_to] {
            if self.occ[empty_neighbor as usize] == 1 {
                self.possible_moves.remove_item(empty_neighbor, move_to);
            }
            if self.occ[empty_neighbor as usize] == 0 {
                // greater than one because of neighbor moving in this spot
                if self.cn_metal[empty_neighbor as usize] > 1 {
                    self.possible_moves.add_item(move_to, empty_neighbor, None);
                }
            }
        }
    }

    fn cond_snap_and_heat_map(&mut self, iiter: &u64) {
        const NUMBER_HEAT_MAP_SECTIONS: u64 = 200;

        if let Some(snap_shot_sections) = &mut self.snap_shot_sections {
            if (iiter + 1) % (self.niter / NUMBER_HEAT_MAP_SECTIONS) == 0 {
                if iiter == &0 {
                    return;
                }
                let mut t_vec = vec![0; self.occ.len()];
                t_vec
                    .iter_mut()
                    .enumerate()
                    .for_each(|(i, x)| *x = self.occ[i]);
                snap_shot_sections.push(t_vec);
            }
        }

        if let Some(heat_map) = &mut self.heat_map {
            if (iiter + 1) % (self.niter / NUMBER_HEAT_MAP_SECTIONS) == 0 {
                if iiter == &0 {
                    return;
                }
                let mut t_vec = vec![0; heat_map.len()];
                t_vec
                    .iter_mut()
                    .enumerate()
                    .for_each(|(i, x)| *x = heat_map[i]);
                self.heat_map_sections.push(t_vec);
            }
        }
    }
}
fn no_int_nn_from_move(
    move_from: u32,
    move_to: u32,
    nn_pair_no_intersec: &HashMap<u64, [[u32; NN_PAIR_NO_INTERSEC_NUMBER]; 2], fnv::FnvBuildHasher>,
) -> ([u32; 7], [u32; 7]) {
    let no_int = nn_pair_no_intersec[&(std::cmp::min(move_from, move_to) as u64
        + ((std::cmp::max(move_to, move_from) as u64) << 32))];
    if move_to > move_from {
        (no_int[0], no_int[1])
    } else {
        (no_int[1], no_int[0])
    }
}

fn no_int_nnn_from_move(
    move_from: u32,
    move_to: u32,
    nnn_pair_no_intersec: &HashMap<
        u64,
        (
            Vec<Vec<u32>>,
            Vec<Vec<u32>>,
            Vec<(u32, Vec<u32>, Vec<u32>, Vec<u32>)>,
        ),
        fnv::FnvBuildHasher,
    >,
) -> (
    &Vec<Vec<u32>>,
    &Vec<Vec<u32>>,
    &Vec<(u32, Vec<u32>, Vec<u32>, Vec<u32>)>,
    bool,
) {
    let (min, max, inter) = &nnn_pair_no_intersec[&(std::cmp::min(move_from, move_to) as u64
        + ((std::cmp::max(move_to, move_from) as u64) << 32))];
    if move_to > move_from {
        (&min, &max, &inter, false)
    } else {
        (&max, &min, &inter, true)
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

// #[cfg(test)]
// mod tests {
//     // Note this useful idiom: importing names from outer (for mod tests) scope.
//     use super::*;
//
//     #[test]
//     fn test_perform_move() {
//         fn file_paths(
//             grid_folder: String,
//         ) -> (String, String, String, String, String, String, String) {
//             (
//                 format!("{}/nearest_neighbor", grid_folder),
//                 format!("{}/next_nearest_neighbor", grid_folder),
//                 format!("{}/nn_pairlist", grid_folder),
//                 format!("{}/nnn_pairlist", grid_folder),
//                 format!("{}/atom_sites", grid_folder),
//                 format!("{}/nn_pair_no_intersec", grid_folder),
//                 format!("{}/nnn_gcn_no_intersec.json", grid_folder),
//             )
//         }
//         let (
//             pairlist_file,
//             n_pairlist_file,
//             nn_pairlist_file,
//             nnn_pairlist_file,
//             atom_sites,
//             nn_pair_no_int_file,
//             nnn_pair_no_int_file,
//         ) = file_paths("../666-pair".to_string());
//
//         let energy = EnergyInput::Gcn([
//             4752, 4719, 4686, 4653, 4620, 4587, 4554, 4521, 4488, 4455, 4422, 4389, 4356, 4323,
//             4290, 4257, 4224, 4191, 4158, 4125, 4092, 4059, 4026, 3993, 3960, 3927, 3894, 3861,
//             3828, 3795, 3762, 3729, 3696, 3663, 3630, 3597, 3564, 3531, 3498, 3465, 3432, 3399,
//             3366, 3333, 3300, 3267, 3234, 3201, 3168, 3135, 3102, 3069, 3036, 3003, 2970, 2937,
//             2904, 2871, 2838, 2805, 2772, 2739, 2706, 2673, 2640, 2607, 2574, 2541, 2508, 2475,
//             2442, 2409, 2376, 2343, 2310, 2277, 2244, 2211, 2178, 2145, 2112, 2079, 2046, 2013,
//             1980, 1947, 1914, 1881, 1848, 1815, 1782, 1749, 1716, 1683, 1650, 1617, 1584, 1551,
//             1518, 1485, 1452, 1419, 1386, 1353, 1320, 1287, 1254, 1221, 1188, 1155, 1122, 1089,
//             1056, 1023, 990, 957, 924, 891, 858, 825, 792, 759, 726, 693, 660, 627, 594, 561, 528,
//             495, 462, 429, 396, 363, 330, 297, 264, 231, 198, 165, 132, 99, 66, 33, 0,
//         ]);
//
//         let gridstructure = GridStructure::new(
//             pairlist_file,
//             n_pairlist_file,
//             nn_pair_no_int_file,
//             nnn_pair_no_int_file,
//             atom_sites,
//             String::from("../input_cluster/bulk.poscar"),
//         );
//
//         let mut sim = Simulation::new(
//             1000000000,
//             None,
//             Some(20),
//             300.,
//             None,
//             String::from("./sim/"),
//             false,
//             false,
//             0_usize,
//             vec![3, 4],
//             energy,
//             None,
//             &gridstructure,
//             0,
//         );
//         let (from, to, to2) = 'bar: {
//             for (from, to, _) in &sim.possible_moves.moves {
//                 for x in sim.gridstructure.nn[to] {
//                     if sim.gridstructure.nn[from].contains(&x) && sim.occ[x as usize] == 0 {
//                         break 'bar (*from, *to, x);
//                     }
//                 }
//             }
//             (0, 0, 0)
//         };
//         println!("{},{},{}", from, to, to2);
//         // let (from, to) = sim.possible_moves.moves[1];
//         let mut effected_gcn = Vec::new();
//         for o in sim.gridstructure.nn[&from] {
//             for x in sim.gridstructure.nn[&o] {
//                 if !effected_gcn.contains(&x) && sim.cn_metal[x as usize] != 0 {
//                     effected_gcn.push(x)
//                 }
//             }
//         }
//         for o in sim.gridstructure.nn[&to] {
//             for x in sim.gridstructure.nn[&o] {
//                 if !effected_gcn.contains(&x) && sim.cn_metal[x as usize] != 0 {
//                     effected_gcn.push(x)
//                 }
//             }
//         }
//         for o in sim.gridstructure.nn[&to2] {
//             for x in sim.gridstructure.nn[&o] {
//                 if !effected_gcn.contains(&x) && sim.cn_metal[x as usize] != 0 {
//                     effected_gcn.push(x)
//                 }
//             }
//         }
//         let old_gcn = sim
//             .gcn_metal
//             .clone()
//             .into_iter()
//             .enumerate()
//             .filter(|(i, _)| effected_gcn.contains(&(*i as u32)))
//             .map(|(i, x)| (i, x))
//             .collect::<Vec<(usize, usize)>>();
//
//         println!(
//             "nn_from: {:?} nn_to: {:?} nn_to2: {:?}",
//             sim.gridstructure.nn[&from], sim.gridstructure.nn[&to], sim.gridstructure.nn[&to2]
//         );
//         // println!("gcn before: {:?}", sim.gcn_metal);
//         assert!(sim.occ[to2 as usize] == 0, "{}", sim.occ[to2 as usize]);
//         assert!(sim.occ[to as usize] == 0, "{}", sim.occ[to as usize]);
//
//         println!(
//             "occ {:?}",
//             sim.occ
//                 .iter()
//                 .enumerate()
//                 .filter(|(_, x)| **x == 1)
//                 .map(|(i, _)| i)
//                 .collect::<Vec<usize>>()
//         );
//         sim.perform_move(from, to, 0, false);
//         // println!("gcn between: {:?}", sim.gcn_metal);
//         // println!("gcn between: ");
//         sim.perform_move(to, to2, 0, false);
//         // println!("gcn between2: ");
//         println!(
//             "{:?}",
//             sim.gcn_metal
//                 .clone()
//                 .into_iter()
//                 .enumerate()
//                 .filter(|(i, _)| effected_gcn.contains(&(*i as u32)))
//                 .map(|(i, x)| (i, x))
//                 .collect::<Vec<(usize, usize)>>()
//         );
//         sim.perform_move(to2, from, 0, false);
//         // println!("gcn after: {:?}", sim.gcn_metal);
//         let (from_change, to_change, intercet, is_reverse) =
//             no_int_nnn_from_move(from, to, &sim.gridstructure.nnn_pair_no_intersec);
//         let (from_change2, to_change2, intercet2, _) =
//             no_int_nnn_from_move(to, to2, &sim.gridstructure.nnn_pair_no_intersec);
//         let (from_change3, to_change3, intercet3, _) =
//             no_int_nnn_from_move(to2, from, &sim.gridstructure.nnn_pair_no_intersec);
//         println!("from: {}, to: {}, to2: {}\n", from, to, to2);
//         println!(
//             "from: {:?},\n to: {:?},\n inter: {:?} \n",
//             from_change, to_change, intercet
//         );
//         println!(
//             "from: {:?},\n to: {:?},\n inter: {:?} \n",
//             from_change2, to_change2, intercet2
//         );
//         println!(
//             "from: {:?},\n to: {:?},\n inter: {:?} \n",
//             from_change3, to_change3, intercet3
//         );
//         let new_gcn = sim
//             .gcn_metal
//             .clone()
//             .into_iter()
//             .enumerate()
//             .filter(|(i, _)| effected_gcn.contains(&(*i as u32)))
//             .map(|(i, x)| (i, x))
//             .collect::<Vec<(usize, usize)>>();
//         assert_eq!(old_gcn, new_gcn);
//     }
// }
enum FromOrTo {
    From,
    To,
}
