use anyhow;
use core::panic;
use csv::Writer;
use energy::EnergyInput;
use rand;
use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use rand::rngs::SmallRng;
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
const NN_PAIR_NUMBER: usize = 20;
const NN_PAIR_NO_INTERSEC_NUMBER: usize = 7;
const AMOUNT_SECTIONS: usize = 10000;
const SAVE_TH: u64 = 1000;

const GRID_SIZE: [u32; 3] = [20, 20, 20];

const SAVE_ENTIRE_SIM: bool = true;

#[derive(Clone, Default)]
pub struct AtomPosition {
    occ: u8,
    cn_metal: usize,
    nn_support: u8,
    nn: [u32; CN],
}

#[derive(Clone)]
pub struct Simulation {
    atom_pos: Vec<AtomPosition>,
    niter: u64,
    number_all_atoms: u32,
    onlyocc: HashSet<u32, fnv::FnvBuildHasher>,
    possible_moves: listdict::ListDict,
    total_energy_1000: i64,
    cn_dict: [u32; CN + 1],
    cn_dict_at_supp: [u32; CN + 1],
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
    is_supported: bool,
}

fn copy_nn_in_atoms_pos(
    atom_pos: &mut [AtomPosition],
    nn: &HashMap<u32, [u32; CN], fnv::FnvBuildHasher>,
) {
    for (i, atom) in atom_pos.iter_mut().enumerate() {
        atom.nn = nn.get(&(i as u32)).unwrap().clone()
    }
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
        //4
        //111: 12
        //hcp: 8
        let nsites: u32 = GRID_SIZE[0] * GRID_SIZE[1] * GRID_SIZE[2] * 12;
        let mut atom_pos: Vec<AtomPosition> = vec![AtomPosition::default(); nsites as usize];
        let mut cn_dict: [u32; CN + 1] = [0; CN + 1];
        let mut cn_dict_at_supp: [u32; CN + 1] = [0; CN + 1];
        let is_supported = if support_indices.is_some() {
            true
        } else {
            false
        };
        let (onlyocc, number_all_atoms) = if input_file.is_some() {
            let xyz = read_and_write::read_sample(&input_file.unwrap());
            let onlyocc = setup::occ_onlyocc_from_xyz(
                &mut atom_pos,
                &xyz,
                nsites,
                &gridstructure.xsites_positions,
            );
            let number_of_atoms: u32 = onlyocc.len() as u32;
            (onlyocc, number_of_atoms)
        } else if atoms_input.is_some() {
            let number_of_atom = atoms_input.unwrap();
            let onlyocc = setup::create_input_cluster(
                &mut atom_pos,
                &atoms_input.unwrap(),
                &gridstructure.xsites_positions,
                &gridstructure.nn,
                nsites,
                support_indices,
            );
            (onlyocc, number_of_atom)
        } else {
            panic!("gib input atoms or input file");
        };

        for o in 0..nsites {
            let mut neighbors: u8 = 0;
            for o1 in gridstructure.nn[&o].iter() {
                if atom_pos[*o1 as usize].occ == 1 {
                    // cn.entry(o).and_modify(|x| *x += 1).or_insert(1);
                    neighbors += 1;
                }
            }
            atom_pos[o as usize].cn_metal = neighbors as usize;
            if atom_pos[o as usize].occ == 1 {
                if is_supported {
                    if atom_pos[o as usize].nn_support == 1 {
                        cn_dict_at_supp[atom_pos[o as usize].cn_metal] += 1;
                    }
                }
                cn_dict[atom_pos[o as usize].cn_metal] += 1;
            };
        }

        let mut total_energy_1000: i64 = 0;
        let mut possible_moves: listdict::ListDict = listdict::ListDict::new(GRID_SIZE);
        for o in onlyocc.iter() {
            let temp_total_e = total_energy_1000;
            match energy {
                EnergyInput::LinearCn(_) | EnergyInput::Cn(_) => {
                    total_energy_1000 += energy::energy_1000_calculation(
                        &energy,
                        atom_pos[*o as usize].cn_metal,
                        atom_pos[*o as usize].nn_support,
                        support_e,
                    );
                }
            };
            println!(
                "at_supp: {} total_e: {} cn: {}",
                atom_pos[*o as usize].nn_support,
                total_energy_1000 - temp_total_e,
                atom_pos[*o as usize].cn_metal
            );

            for u in &gridstructure.nn[o] {
                if atom_pos[*u as usize].occ == 0 {
                    // >1 so that atoms cant leave the cluster
                    // <x cant move if all neighbors are occupied
                    if atom_pos[*u as usize].cn_metal > 1 {
                        possible_moves.add_item(*o, *u, None)
                    }
                }
            }
        }
        // panic!("nooo");

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

        copy_nn_in_atoms_pos(&mut atom_pos, &gridstructure.nn);
        Simulation {
            atom_pos,
            niter,
            number_all_atoms,
            onlyocc,
            possible_moves,
            total_energy_1000,
            cn_dict,
            cn_dict_at_supp,
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
            is_supported,
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
                    "total cn: {:?}",
                    self.atom_pos.iter().map(|x| x.cn_metal).sum::<usize>()
                );
                println!(
                    "iteration {}; {}%",
                    iiter,
                    (iiter as f64 / self.niter as f64 * 100.)
                );
                // println!("{:?}", self.atom_pos.cn_metal);
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
                self.cn_dict_at_supp.iter_mut().for_each(|x| {
                    *x = 0;
                });
                for o in 0..self.atom_pos.len() {
                    if self.atom_pos[o].occ == 1 {
                        self.update_cn_dict(o, self.atom_pos[o as usize].cn_metal, true);
                        // if let Some(nn_support) = self.atom_pos.nn_support {
                        //     if nn_support[o as usize] == 1 {
                        //         self.cn_dict_at_supp[self.atom_pos.cn_metal[o as usize]] += 1;
                        //     } else {
                        //         self.cn_dict[self.atom_pos.cn_metal[o as usize]] += 1;
                        //     }
                        // } else {
                        //     self.cn_dict[self.atom_pos.cn_metal[o as usize]] += 1;
                        // }
                        // self.cn_dict[self.atom_pos.cn_metal[o]] += 1;
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
            &self.atom_pos,
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
            let empty_neighbor_cn = self.count_empty_sites();
            lowest_energy_struct.empty_cn = empty_neighbor_cn;
            lowest_energy_struct.energy = self.total_energy_1000 as f64 / 1000.;
            lowest_energy_struct.iiter = *iiter;

            let mut cn_hash_map: HashMap<u8, u32> = HashMap::new();
            for (i, v) in self.cn_dict.into_iter().enumerate() {
                cn_hash_map.insert(i as u8, v);
            }
            lowest_energy_struct.cn_total = cn_hash_map;

            let mut cn_hash_map_at_supp: HashMap<u8, u32> = HashMap::new();
            for (i, v) in self.cn_dict_at_supp.into_iter().enumerate() {
                cn_hash_map_at_supp.insert(i as u8, v);
            }
            lowest_energy_struct.cn_dict_at_supp = cn_hash_map_at_supp;

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
        self.atom_pos[move_to as usize].occ = self.atom_pos[move_from as usize].occ; // covers different alloys also
        self.atom_pos[move_from as usize].occ = 0;

        self.onlyocc.remove(&move_from);
        self.onlyocc.insert(move_to);

        if SAVE_ENTIRE_SIM || is_recording_sections {
            // self.cn_dict[self.atom_pos.cn_metal[move_from as usize]] -= 1;
            self.update_cn_dict(
                move_from as usize,
                self.atom_pos[move_from as usize].cn_metal,
                false,
            );
        }
        // let (from_change, to_change) = self.no_int_from_move(move_from, move_to);
        for o in self.atom_pos[move_from as usize].nn {
            if (SAVE_ENTIRE_SIM || is_recording_sections)
                && self.atom_pos[o as usize].occ == 1
                && o != move_to
            {
                self.update_cn_dict(o as usize, self.atom_pos[o as usize].cn_metal, false);
                self.update_cn_dict(o as usize, self.atom_pos[o as usize].cn_metal - 1, true);
                // self.cn_dict[self.atom_pos.cn_metal[o as usize]] -= 1;
                // self.cn_dict[self.atom_pos.cn_metal[o as usize] - 1] += 1;
            }
            self.atom_pos[o as usize].cn_metal -= 1;
        }
        for o in self.atom_pos[move_to as usize].nn {
            if (SAVE_ENTIRE_SIM || is_recording_sections)
                && self.atom_pos[o as usize].occ == 1
                && o != move_from
            {
                self.update_cn_dict(o as usize, self.atom_pos[o as usize].cn_metal, false);
                self.update_cn_dict(o as usize, self.atom_pos[o as usize].cn_metal + 1, true);
                // self.cn_dict[self.atom_pos.cn_metal[o as usize]] -= 1;
                // self.cn_dict[self.atom_pos.cn_metal[o as usize] + 1] += 1;
            }
            self.atom_pos[o as usize].cn_metal += 1;
        }

        if SAVE_ENTIRE_SIM || is_recording_sections {
            self.update_cn_dict(
                move_to as usize,
                self.atom_pos[move_to as usize].cn_metal,
                true,
            );
            // self.cn_dict[self.atom_pos.cn_metal[move_to as usize]] += 1;
        }

        self.total_energy_1000 += energy1000_diff;
    }

    fn energy_change_by_move(&self, move_from: u32, move_to: u32) -> i64 {
        let (from_at_support, to_at_support) = if self.is_supported {
            let from_at_support = self.atom_pos[move_from as usize].nn_support;
            let to_at_support = self.atom_pos[move_to as usize].nn_support;
            (from_at_support, to_at_support)
        } else {
            (0, 0)
        };

        match &self.energy {
            EnergyInput::LinearCn(energy_l_cn) => energy::energy_diff_l_cn(
                energy_l_cn.complet_energy,
                self.atom_pos[move_from as usize].cn_metal,
                self.atom_pos[move_to as usize].cn_metal - 1,
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
                    from_change.iter().filter_map(|x| {
                        if self.atom_pos[*x as usize].occ == 1 {
                            Some((
                                self.atom_pos[*x as usize].cn_metal,
                                self.atom_pos[*x as usize].nn_support,
                            ))
                        } else {
                            None
                        }
                    }),
                    // .filter(|x| self.atom_pos[**x as usize].occ == 1)
                    // .map(|x| {
                    //     (
                    //         self.atom_pos[*x as usize].cn_metal,
                    //         self.atom_pos[*x as usize].nn_support,
                    //     )
                    // }),
                    to_change
                        .iter()
                        // .filter(|x| self.atom_pos[**x as usize].occ == 1)
                        .filter_map(|x| {
                            if self.atom_pos[*x as usize].occ == 1 {
                                Some((
                                    self.atom_pos[*x as usize].cn_metal,
                                    self.atom_pos[*x as usize].nn_support,
                                ))
                            } else {
                                None
                            }
                        }),
                    // .map(|x| (self.atom_pos.cn_metal[*x as usize], self.atom_pos.nn_support[*x as usize])),
                    self.atom_pos[move_from as usize].cn_metal,
                    self.atom_pos[move_to as usize].cn_metal,
                    from_at_support,
                    to_at_support,
                    self.support_e,
                )
            }
        }
    }

    fn update_possible_moves(&mut self, move_from: u32, move_to: u32) {
        self.possible_moves.remove_item(move_from, move_to);
        for neighbor_atom in self.atom_pos[move_from as usize].nn {
            if self.atom_pos[neighbor_atom as usize].occ == 0 {
                self.possible_moves.remove_item(move_from, neighbor_atom);
            }
            if self.atom_pos[neighbor_atom as usize].occ == 1 {
                // greater than one because of neighbor moving in this spot
                if self.atom_pos[move_from as usize].cn_metal > 1 {
                    self.possible_moves.add_item(neighbor_atom, move_from, None);
                }
            }
        }

        for empty_neighbor in self.atom_pos[move_to as usize].nn {
            if self.atom_pos[empty_neighbor as usize].occ == 1 {
                self.possible_moves.remove_item(empty_neighbor, move_to);
            }
            if self.atom_pos[empty_neighbor as usize].occ == 0 {
                // greater than one because of neighbor moving in this spot
                if self.atom_pos[empty_neighbor as usize].cn_metal > 1 {
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
                let mut t_vec = vec![0; self.atom_pos.len()];
                t_vec
                    .iter_mut()
                    .enumerate()
                    .for_each(|(i, x)| *x = self.atom_pos[i].occ);
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

    #[inline]
    fn update_cn_dict(&mut self, atom: usize, cn: usize, change_is_positiv: bool) {
        match change_is_positiv {
            true => {
                if self.is_supported {
                    if self.atom_pos[atom].nn_support == 1 {
                        self.cn_dict_at_supp[cn] += 1;
                    }
                }
                self.cn_dict[cn] += 1;
            }
            false => {
                if self.is_supported {
                    if self.atom_pos[atom].nn_support == 1 {
                        self.cn_dict_at_supp[cn] -= 1;
                    }
                }
                self.cn_dict[cn] -= 1;
            }
        }
    }

    pub fn count_empty_sites(&self) -> HashMap<String, u32> {
        // let empty_set: HashSet<&u32> =
        //     HashSet::from_iter(self.possible_moves.iter().map(|(_, empty, _)| empty));
        let mut empty_sites = HashSet::new();
        let mut empty_sites_distribution: HashMap<String, u32> = HashMap::new();
        for atom in self.onlyocc.iter() {
            // if occupied != &1 {
            //     continue;
            // }
            for neigbor in self.atom_pos[*atom as usize].nn.iter() {
                if self.atom_pos[*neigbor as usize].occ == 0 {
                    empty_sites.insert(neigbor);
                }
            }
        }
        for site in empty_sites.into_iter() {
            let mut neigbors_count = 0_u32;
            for neighbor in self.atom_pos[*site as usize].nn {
                if self.atom_pos[neighbor as usize].occ == 1 {
                    neigbors_count += 1;
                }
            }
            if neigbors_count == 5 {
                let mut cn_ten_count = 0;
                let mut cn_seven_count = 0;
                for outer in self.atom_pos[*site as usize].nn {
                    if self.atom_pos[outer as usize].cn_metal == 10 {
                        cn_ten_count += 1;
                    }
                    if self.atom_pos[outer as usize].cn_metal == 7 {
                        cn_seven_count += 1;
                    }
                }
                if cn_ten_count >= 1 {
                    empty_sites_distribution
                        .entry("b5a".to_string())
                        .and_modify(|counter| *counter += 1)
                        .or_insert(1);
                }
            }
            empty_sites_distribution
                .entry(neigbors_count.to_string())
                .and_modify(|counter| *counter += 1)
                .or_insert(1);
        }
        empty_sites_distribution
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

enum FromOrTo {
    From,
    To,
}
