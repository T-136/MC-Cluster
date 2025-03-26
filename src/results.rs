use serde::{Deserialize, Serialize, Serializer};
use serde_with::serde_as;
use std::collections::{BTreeMap, HashMap, HashSet};

use super::Simulation;

// const KB: f64 = 8.6173324e-5;

#[derive(Serialize, Deserialize, Default)]
pub struct LowestEnergy {
    pub energy: f64,
    #[serde(serialize_with = "ordered_map")]
    pub cn_total: HashMap<u8, u32>,
    // #[serde(serialize_with = "ordered_map")]
    pub empty_cn: HashMap<String, u32>,
    #[serde(serialize_with = "ordered_map")]
    pub cn_dict_at_supp: HashMap<u8, u32>,
    pub iiter: u64,
    #[serde(skip_serializing)]
    pub onlyocc: HashSet<u32, fnv::FnvBuildHasher>,
}

impl LowestEnergy {
    pub fn new() -> LowestEnergy {
        LowestEnergy {
            energy: f64::INFINITY,
            ..Default::default()
        }
    }

    pub fn update(&mut self, sim: &Simulation, iiter: &u64) -> bool {
        if self.energy > (sim.total_energy_1000 as f64 / 1000.) {
            let empty_neighbor_cn = sim.count_empty_sites(&sim.onlyocc);
            self.empty_cn = empty_neighbor_cn;
            self.energy = sim.total_energy_1000 as f64 / 1000.;
            self.iiter = *iiter;

            let mut cn_hash_map: HashMap<u8, u32> = HashMap::new();
            for (i, v) in sim.cn_dict.into_iter().enumerate() {
                cn_hash_map.insert(i as u8, v);
            }
            self.cn_total = cn_hash_map;

            let mut cn_hash_map_at_supp: HashMap<u8, u32> = HashMap::new();
            for (i, v) in sim.cn_dict_at_supp.into_iter().enumerate() {
                cn_hash_map_at_supp.insert(i as u8, v);
            }
            self.cn_dict_at_supp = cn_hash_map_at_supp;
            self.onlyocc = sim.onlyocc.clone();

            true
        } else {
            false
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct Start {
    pub start_energy: f64,
    #[serde(serialize_with = "ordered_map")]
    pub start_cn: HashMap<u8, u32>,
}

impl Start {
    pub fn new(total_energy_1000: i64, cn_dict: &[u32]) -> Start {
        let start_energy = total_energy_1000 as f64 / 1000.;
        let mut start_cn_dict = HashMap::new();
        for (k, v) in cn_dict.iter().enumerate() {
            start_cn_dict.insert(k as u8, *v);
        }
        Start {
            start_energy,
            start_cn: start_cn_dict,
        }
    }
}

#[serde_as]
#[derive(Serialize, Deserialize)]
pub struct Results {
    pub start: Start,
    pub lowest_energy_struct: LowestEnergy,
    pub number_all_atoms: u32,
    pub energy_section_list: Vec<f64>,
    pub cn_dict_sections: Vec<HashMap<u8, f64>>,
}

fn ordered_map<S>(value: &HashMap<u8, u32>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let ordered: BTreeMap<_, _> = value.iter().collect();
    ordered.serialize(serializer)
}
