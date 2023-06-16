use chemfiles::{Atom, Frame, Trajectory, UnitCell};
use rand;
use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use rand_chacha::ChaCha20Rng;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::{cmp, eprint, fs, println};

pub use super::sim::{Results, Seed};
use super::{listdict, sim};

#[derive(Clone)]
pub struct Simulation {
    pub niter: u64,
    number_all_atoms: u32,
    occ: Vec<u8>,
    onlyocc: HashSet<u32>,
    cn: Vec<usize>,
    former_energy_dict: HashMap<u32, i64>,
    pub possible_moves: listdict::ListDict,
    pub total_energy_1000: i64,
    nn: HashMap<u32, [u32; 12], fnv::FnvBuildHasher>,
    nn_pair: HashMap<u64, [u32; 20], fnv::FnvBuildHasher>,
    // nnn_pair: HashMap<u64, [u32; 74], fnv::FnvBuildHasher>,
    // mut xyz: Vec<[f64; 3]>,
    xsites_positions: Vec<[f64; 3]>,
    unit_cell: UnitCell,
    cn_dict: [u32; 13],
    pub save_folder: String,
    trajectory_frequency: Option<u64>,
    last_frames_trajectory: Option<u64>,
    start_temperature: Option<f64>,
    temperature: f64,
    cn_dict_sections: Vec<HashMap<u8, f64>>,
    energy_sections_list: Vec<f64>,
    optimization_cut_off_perc: f64,
    unique_levels: HashMap<BTreeMap<u8, u32>, (i64, u64)>,
}

impl Simulation {
    pub fn create_start_struct(&self) -> sim::Start {
        let start_energy = self.total_energy_1000 as f64 / 1000.;

        let mut start_cn_dict = HashMap::new();

        for (k, v) in self.cn_dict.iter().enumerate() {
            start_cn_dict.insert(k as u8, *v);
        }

        sim::Start {
            start_energy: start_energy,
            start_cn: start_cn_dict,
        }
    }
    pub fn is_acceptance_criteria_fulfilled(
        &mut self,
        proposed_energy: i64,
        rng_e_number: &mut ChaCha20Rng,
        iiter: u64,
    ) -> bool {
        const KB: f64 = 8.6173324e-5;
        // if self.start_temperature.is_some() {
        if proposed_energy < self.total_energy_1000 {
            return true;
        }
        let acceptance_temp = self.calculate_current_temp(iiter);
        let between = Uniform::new_inclusive(0., 1.);
        let delta_energy = proposed_energy - self.total_energy_1000;
        let rand_value = between.sample(rng_e_number);
        (rand_value) < ((-delta_energy as f64 / 1000.) / (KB * acceptance_temp)).exp()
        // } else {
        //     if &proposed_energy < &self.total_energy {
        //         return true;
        //     }
        //     let acceptance_temp = 5000.
        //         - (iiter as f64 / self.niter as f64)
        //             * (self.start_temperature.unwrap() - self.temperature);
        //     let between = Uniform::new_inclusive(0., 1.);
        //     let delta_energy = proposed_energy - self.total_energy;
        //     let rand_value = between.sample(rng_e_number);
        //     if (rand_value) < (-delta_energy / (KB * acceptance_temp)).exp() {
        //         return true;
        //     } else {
        //         return false;
        //     }
        // }
    }
    pub fn return_results(
        &self,
        start: sim::Start,
        lowest_energy_struct: sim::LowestEnergy,
        seed: Seed,
    ) -> Results {
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
    pub fn save_sections(
        &mut self,
        iiter: &u64,
        mut temp_energy_section_1000: i64,
        temp_cn_dict_section: &mut [u32; 13],
        amount_unique_levels: &mut i32,
    ) -> i64 {
        const SECTION_SIZE: u64 = 100000;
        temp_energy_section_1000 += self.total_energy_1000;

        temp_cn_dict_section
            .iter_mut()
            .enumerate()
            .for_each(|(i, v)| *v += self.cn_dict[i]);

        let mut cn_hash_map = HashMap::new();
        for (i, v) in self.cn_dict.into_iter().enumerate() {
            cn_hash_map.insert(i as u8, v);
        }

        if *amount_unique_levels != 0 {
            if *iiter as f64 >= self.niter as f64 * self.optimization_cut_off_perc {
                let cn_btree: BTreeMap<_, _> = cn_hash_map.into_iter().collect();
                match self.unique_levels.entry(cn_btree) {
                    std::collections::hash_map::Entry::Occupied(mut entry) => {
                        let (_, x) = entry.get_mut();
                        *x += 1;
                    }
                    std::collections::hash_map::Entry::Vacant(entry) => {
                        if *amount_unique_levels == 1 {
                            eprint!("amount_unique_levels reached");
                        }
                        *amount_unique_levels -= 1;
                        entry.insert((self.total_energy_1000 / 1000, 1));
                    }
                }
                self.unique_levels
                    .entry(cn_btree)
                    .and_modify(|(_, x)| *x += 1)
                    .or_insert_with(|| {
                        if *amount_unique_levels == 1 {
                            eprint!("amount_unique_levels reached");
                        }
                        *amount_unique_levels -= 1;
                        (self.total_energy_1000 / 1000, 1)
                    });
            }
        }
        if (iiter + 1) % SECTION_SIZE == 0 {
            self.energy_sections_list
                .push(temp_energy_section_1000 as f64 / SECTION_SIZE as f64 / 1000.);
            temp_energy_section_1000 = 0;
            // temp_energy_section.clear();

            let mut section: HashMap<u8, f64> = HashMap::new();
            for (k, list) in temp_cn_dict_section.iter_mut().enumerate() {
                // let mut temp_cn_summ = 0;
                // for v1 in list.into_iter() {
                //     temp_cn_summ += v1.clone();
                // }
                section.insert(k as u8, *list as f64 / SECTION_SIZE as f64);
                // list.clear();
                *list = 0;
            }
            assert_eq!(temp_cn_dict_section, &mut [0_u32; 13]);
            self.cn_dict_sections.push(section.clone())
        }
        temp_energy_section_1000
    }

    pub fn save_lowest_energy(
        &mut self,
        iiter: &u64,
        lowest_energy_struct: &mut sim::LowestEnergy,
    ) {
        if *iiter as f64 >= self.niter as f64 * self.optimization_cut_off_perc {
            if &lowest_energy_struct.energy > &(self.total_energy_1000 as f64 / 1000.) {
                let mut empty_neighbor_cn: HashMap<u8, u32> = HashMap::new();
                let empty_set: HashSet<&u32> =
                    HashSet::from_iter(self.possible_moves.iter().map(|(_, empty)| empty));
                for empty in empty_set {
                    if self.cn[*empty as usize] > 3 {
                        empty_neighbor_cn
                            .entry(self.cn[*empty as usize] as u8)
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
        };
    }

    pub fn write_trajectorys(
        &self,
        iiter: &u64,
        // xyz: &mut Vec<[f64; 3]>,
        trajectory_option: &mut Option<Trajectory>,
        trajectory_last_frames_option: &mut Option<Trajectory>,
    ) {
        // let mut xyz: Vec<[f64; 3]> = Vec::new();
        if let Some(trajectory_last_frames) = trajectory_last_frames_option {
            if self.niter - iiter <= self.last_frames_trajectory.unwrap() {
                self.write_traj(trajectory_last_frames);
            }
        }
        // let mut xyz: Vec<[f64; 3]> = Vec::new();
        if let Some(trajectory) = trajectory_option {
            if self.niter - iiter > self.trajectory_frequency.unwrap_or(0) && iiter % 100 == 0 {
                self.write_traj(trajectory);
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

        trajectory.write(&mut frame).unwrap();
    }
    pub fn temp_energy_calculation(
        &self,
        move_from: u32,
        move_to: u32,
        total_temp_energy: &mut i64,
    ) {
        let lower_position = cmp::min(&move_from, &move_to).clone();
        let higher_position = cmp::max(&move_from, &move_to).clone();
        for o in self
            .nn_pair
            .get(&(lower_position as u64 + ((higher_position as u64) << 32)))
            .unwrap()
        {
            if self.occ[*o as usize] != 0 {
                *total_temp_energy += sim::energy_calculation(o, &self.cn);
                if o == &move_to {
                    continue;
                }
                *total_temp_energy -= self.former_energy_dict[o];
            }
        }
        *total_temp_energy -= self.former_energy_dict[&move_from];
    }

    pub fn perform_move(&mut self, move_from: u32, move_to: u32) -> u8 {
        self.occ[move_to as usize] = self.occ[move_from as usize]; // covers different alloys also
        self.occ[move_from as usize] = 0;

        self.onlyocc.remove(&move_from);
        self.onlyocc.insert(move_to);

        let mut cn_change: i32 = 0;

        self.cn_dict[self.cn[move_from as usize]] -= 1;
        for o in &self.nn[&move_from] {
            if self.occ[*o as usize] == 1 && o != &move_to {
                self.cn_dict[self.cn[*o as usize]] -= 1;
                self.cn_dict[self.cn[*o as usize] - 1] += 1;
                cn_change -= 1;
            }
            self.cn[*o as usize] -= 1;
        }
        for o in &self.nn[&move_to] {
            if self.occ[*o as usize] == 1 && o != &move_from {
                self.cn_dict[self.cn[*o as usize] as usize] -= 1;
                self.cn_dict[(self.cn[*o as usize] + 1) as usize] += 1;
                cn_change += 1
            }
            self.cn[*o as usize] += 1;
        }

        self.cn_dict[self.cn[move_to as usize]] += 1;

        cn_change.unsigned_abs() as u8
    }

    pub fn accept_move(&mut self, total_temp_energy: i64, move_from: u32, move_to: u32) {
        let lower_position = cmp::min(&move_from, &move_to).clone();
        let higher_position = cmp::max(&move_from, &move_to).clone();
        for o in &self.nn_pair[&(lower_position as u64 + ((higher_position as u64) << 32))] {
            self.former_energy_dict
                .insert(*o, sim::energy_calculation(o, &self.cn));
        }
        self.total_energy_1000 = total_temp_energy;

        self.update_possible_moves(move_from, move_to)
    }

    fn update_possible_moves(&mut self, move_from: u32, move_to: u32) {
        for neighbor_atom in self.nn[&move_from] {
            self.possible_moves.remove_item(move_from, neighbor_atom);
            if self.occ[neighbor_atom as usize] != 0 {
                // greater than one because of neighbor moving in this spot
                if self.cn[move_from as usize] > 1 {
                    self.possible_moves.add_item(neighbor_atom, move_from)
                }
            }
        }

        for empty_neighbor in self.nn[&move_to] {
            self.possible_moves.remove_item(empty_neighbor, move_to);
            if self.occ[empty_neighbor as usize] == 0 {
                // greater than one because of neighbor moving in this spot
                if self.cn[empty_neighbor as usize] > 1 {
                    self.possible_moves.add_item(move_to, empty_neighbor)
                }
            }
        }
    }
    pub fn calculate_current_temp(&self, iiter: u64) -> f64 {
        let heating_temp = 5500.;
        let cut_off = self.optimization_cut_off_perc;
        if self.start_temperature.is_some() {
            if (iiter + 1) as f64 <= self.niter as f64 * cut_off {
                heating_temp
                    - ((iiter + 1) as f64 / (self.niter as f64 * cut_off))
                        * (heating_temp - self.start_temperature.unwrap())
            } else {
                self.start_temperature.unwrap()
                    - ((iiter + 1) as f64 - self.niter as f64 * cut_off)
                        / (self.niter as f64 * (1. - cut_off))
                        * (self.start_temperature.unwrap() - self.temperature)
            }
        } else {
            if (iiter + 1) as f64 <= self.niter as f64 * cut_off {
                heating_temp
                    - (iiter as f64 / (self.niter as f64 * cut_off))
                        * (heating_temp - self.temperature)
            } else {
                self.temperature
            }
        }
    }
}
