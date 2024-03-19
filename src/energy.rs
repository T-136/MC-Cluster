use rayon::array;
use std::ops::Index;

// Pt
const M_BETA: i64 = -330;
const M_ALPHA: i64 = 3960;

// const support_e: i64 = -211;
#[derive(Clone, Debug)]
pub struct EnergyValues<T> {
    pub complet_energy: T,
    pub co_ads_energy: Option<T>,
}

#[derive(Clone, Debug)]
pub enum EnergyInput {
    LinearCn(EnergyValues<[i64; 2]>),
    Cn(EnergyValues<[i64; 13]>),
    LinearGcn(EnergyValues<[i64; 2]>),
    Gcn(EnergyValues<[i64; 145]>),
}

pub fn energy_1000_calculation(
    energy: &EnergyInput,
    cn: usize,
    at_support: u8,
    support_e: i64,
) -> i64 {
    match energy {
        EnergyInput::LinearCn(e) => {
            e.complet_energy[0] * cn as i64 + e.complet_energy[1] + support_e * at_support as i64
        }
        EnergyInput::Cn(e) => {
            if let Some(co_ads) = e.co_ads_energy {
                e.complet_energy[cn]
                    + co_ads[cn] * (at_support as i64)
                    + support_e * (at_support as i64)
            } else {
                e.complet_energy[cn] + support_e * (at_support as i64)
            }
        }
        EnergyInput::LinearGcn(e) => {
            todo!();
            // e[0] * cn as i64 + e[1]
        }
        EnergyInput::Gcn(e) => {
            if let Some(co_ads) = e.co_ads_energy {
                e.complet_energy[cn] - co_ads[cn] + support_e * at_support as i64
            } else {
                e.complet_energy[cn] + support_e * at_support as i64
            }
        }
    }
    // enrico_table(cn[*atom as usize])
    // cn[*atom as usize] as i64 * M_BETA + M_ALPHA
}

pub fn energy_diff_cn<'a, T, I, O>(
    energy: &EnergyValues<T>,
    cn_from_list: I,
    cn_to_list: O,
    move_from_cn: usize,
    move_to_cn: usize,
    is_from_at_support: u8,
    is_to_at_support: u8,
    support_e: i64,
) -> i64
where
    T: Index<usize, Output = i64>,
    I: Iterator<Item = (usize, u8)>,
    O: Iterator<Item = (usize, u8)>,
{
    if let Some(co_ads) = energy.co_ads_energy.as_ref() {
        let mut energy_diff_1000 = 0;
        for (cn_from, is_at_supp) in cn_from_list {
            energy_diff_1000 -= energy.complet_energy[cn_from];
            energy_diff_1000 += energy.complet_energy[cn_from - 1];
            if is_at_supp == 1 {
                energy_diff_1000 -= co_ads[cn_from];
                energy_diff_1000 += co_ads[cn_from - 1];
            }
        }
        for (cn_to, is_at_supp) in cn_to_list {
            energy_diff_1000 -= energy.complet_energy[cn_to];
            energy_diff_1000 += energy.complet_energy[cn_to + 1];
            if is_at_supp == 1 {
                energy_diff_1000 -= co_ads[cn_to];
                energy_diff_1000 += co_ads[cn_to + 1];
            }
        }

        energy_diff_1000 -= energy.complet_energy[move_from_cn];
        if is_from_at_support == 1 {
            energy_diff_1000 -= support_e;
            energy_diff_1000 -= co_ads[move_from_cn];
        }
        energy_diff_1000 += energy.complet_energy[move_to_cn - 1];
        if is_to_at_support == 1 {
            energy_diff_1000 += support_e;
            energy_diff_1000 += co_ads[move_to_cn - 1];
        }
        // println!("{},{}", is_from_at_support, is_to_at_support);
        // println!(
        //     "e_diff: {} at_supp_form: {}, at_supp_to: {}",
        //     energy_diff_1000, is_from_at_support, is_to_at_support
        // );
        energy_diff_1000
    } else {
        let mut energy_diff_1000 = 0;
        for (cn_from, is_at_supp) in cn_from_list {
            energy_diff_1000 -= energy.complet_energy[cn_from];
            energy_diff_1000 += energy.complet_energy[cn_from - 1];
        }
        for (cn_to, is_at_supp) in cn_to_list {
            energy_diff_1000 -= energy.complet_energy[cn_to];
            energy_diff_1000 += energy.complet_energy[cn_to + 1];
        }
        energy_diff_1000 -= energy.complet_energy[move_from_cn];
        if is_from_at_support == 1 {
            energy_diff_1000 -= support_e;
        }
        energy_diff_1000 += energy.complet_energy[move_to_cn - 1];
        if is_to_at_support == 1 {
            energy_diff_1000 += support_e;
        }

        energy_diff_1000
    }
}

#[inline(always)]
fn add_energy<T>(energy: &EnergyValues<T>, cn_old: usize, cn_new: usize) -> i64
where
    T: Index<usize, Output = i64>,
{
    energy.complet_energy[cn_new] - energy.complet_energy[cn_old]
}

#[inline(always)]
fn add_energy_supp<T>(energy: &EnergyValues<T>, cn_old: usize, cn_new: usize) -> i64
where
    T: Index<usize, Output = i64>,
{
    energy.complet_energy[cn_new] + energy.co_ads_energy.as_ref().unwrap()[cn_new]
        - (energy.complet_energy[cn_old] + energy.co_ads_energy.as_ref().unwrap()[cn_old])
}

pub fn energy_diff_gcn<T, I, O, P>(
    energy: &EnergyValues<T>,
    cn_from_list: I,
    cn_to_list: O,
    intersect_change: Option<P>,
    move_from_gcn: usize,
    move_to_gcn: usize,
    from_at_support: u8,
    to_at_support: u8,
    // support_e: i64,
) -> i64
where
    T: Index<usize, Output = i64>,
    // <T as Index<usize>>::Output: i64,
    // T: Index<usize>,
    I: Iterator<Item = (usize, usize)>,
    O: Iterator<Item = (usize, usize)>,
    P: Iterator<Item = (usize, usize)>,
{
    if energy.co_ads_energy.is_some() {
        let mut energy_diff_1000 = 0;
        for (gcn_from_old, gcn_from_new) in cn_from_list {
            // energy_diff_1000 -= energy.complet_energy[gcn_from_old];
            energy_diff_1000 += add_energy_supp(&energy, gcn_from_old, gcn_from_new);
        }
        for (gcn_to_old, gcn_to_new) in cn_to_list {
            // energy_diff_1000 -= energy.complet_energy[gcn_to_old];
            // energy_diff_1000 += energy.complet_energy[gcn_to_new];
            energy_diff_1000 += add_energy_supp(&energy, gcn_to_old, gcn_to_new);
        }
        if let Some(intersect_change) = intersect_change {
            for (gcn_inter_old, gcn_inter_new) in intersect_change {
                // energy_diff_1000 -= energy.complet_energy[gcn_inter_old];
                // energy_diff_1000 += energy.complet_energy[gcn_inter_new];
                energy_diff_1000 += add_energy_supp(&energy, gcn_inter_old, gcn_inter_new);
            }
        }
        energy_diff_1000 += add_energy_supp(&energy, move_from_gcn, move_to_gcn);
        // energy_diff_1000 += energy.complet_energy[move_to_gcn];
        // energy_diff_1000 -= energy.complet_energy[move_from_gcn];
        // energy_diff_1000 -= support_e * from_at_support as i64;
        // energy_diff_1000 += support_e * to_at_support as i64;

        energy_diff_1000
    } else {
        let mut energy_diff_1000 = 0;
        for (gcn_from_old, gcn_from_new) in cn_from_list {
            energy_diff_1000 += add_energy(&energy, gcn_from_old, gcn_from_new);
        }
        for (gcn_to_old, gcn_to_new) in cn_to_list {
            energy_diff_1000 += add_energy(&energy, gcn_to_old, gcn_to_new);
        }
        if let Some(intersect_change) = intersect_change {
            for (gcn_inter_old, gcn_inter_new) in intersect_change {
                energy_diff_1000 += add_energy(&energy, gcn_inter_old, gcn_inter_new);
            }
        }
        energy_diff_1000 += add_energy(&energy, move_to_gcn, move_from_gcn);
        // energy_diff_1000 += energy.complet_energy[move_to_gcn];
        // energy_diff_1000 -= energy.complet_energy[move_from_gcn];

        energy_diff_1000
    }
}
// }

pub fn energy_diff_l_cn(
    energy: [i64; 2],
    cn_from: usize,
    cn_to: usize,
    from_at_support: u8,
    to_at_support: u8,
    support_e: i64,
) -> i64 {
    let e = (2 * ((cn_to as i64) * energy[0] + energy[1])) + (support_e * to_at_support as i64)
        - (2 * ((cn_from as i64) * energy[0] + energy[1]))
        - (support_e * from_at_support as i64);
    e
}
pub fn energy_diff_l_gcn<I, O>(
    energy: [i64; 2],
    cn_from_list: I,
    cn_to_list: O,
    move_from_cn: usize,
    move_to_cn: usize,
) -> i64
where
    I: Iterator<Item = usize>,
    O: Iterator<Item = usize>,
{
    let mut energy_diff_1000 = 0;
    for cn_from in cn_from_list {
        energy_diff_1000 -= energy[0] * cn_from as i64 + energy[1];
        energy_diff_1000 += energy[0] * (cn_from as i64 - 1) + energy[1];
    }
    for cn_to in cn_to_list {
        energy_diff_1000 -= energy[0] * cn_to as i64 + energy[1];
        energy_diff_1000 += energy[0] * (cn_to as i64 + 1) + energy[1];
    }
    // energy_diff_1000 -= energy[move_from_cn];
    // energy_diff_1000 += energy[move_to_cn - 1];

    2 * energy_diff_1000
}
