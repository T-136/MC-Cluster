use serde::{Deserialize, Serialize, Serializer};
use serde_with::serde_as;
use std::collections::{BTreeMap, HashMap};

// const KB: f64 = 8.6173324e-5;

#[derive(Serialize, Deserialize)]
pub struct LowestEnergy {
    pub energy: f64,
    #[serde(serialize_with = "ordered_map")]
    pub cn_total: HashMap<u8, u32>,
    // #[serde(serialize_with = "ordered_map")]
    pub empty_cn: HashMap<String, u32>,
    #[serde(serialize_with = "ordered_map")]
    pub cn_dict_at_supp: HashMap<u8, u32>,
    pub iiter: u64,
}

impl Default for LowestEnergy {
    fn default() -> LowestEnergy {
        LowestEnergy {
            energy: f64::INFINITY,
            cn_total: HashMap::new(),
            empty_cn: HashMap::new(),
            cn_dict_at_supp: HashMap::new(),
            iiter: 0,
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
