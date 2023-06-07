use serde::{Deserialize, Serialize};
use std::collections::HashMap;

const M_BETA: i64 = -0330;
const M_ALPHA: i64 = 3960;

// const KB: f64 = 8.6173324e-5;

#[derive(Serialize, Deserialize)]
pub struct LowestEnergy {
    pub energy: f64,
    pub cn_total: HashMap<u8, u32>,
    pub empty_cn: HashMap<u8, u32>,
    pub iiter: u64,
}
#[derive(Serialize, Deserialize)]
pub struct Seed {
    pub rust: String,
    pub choose_seed: [u8; 32],
    pub e_number_seed: [u8; 32],
}
#[derive(Serialize, Deserialize)]
pub struct Start {
    pub start_energy: f64,
    pub start_cn: HashMap<u8, u32>,
}

#[derive(Serialize, Deserialize)]
pub struct Results {
    pub start: Start,
    pub lowest_energy_struct: LowestEnergy,
    pub number_all_atoms: u32,
    pub energy_section_list: Vec<f64>,
    pub cn_dict_sections: Vec<HashMap<u8, f64>>,
    pub seed: Seed,
    pub unique_levels: HashMap<i64, (HashMap<u8, u32>, u64)>,
}

pub fn energy_calculation(atom: &u32, cn: &Vec<usize>) -> i64 {
    let energy_1000 = cn[*atom as usize] as i64 * M_BETA + M_ALPHA;
    energy_1000
}
