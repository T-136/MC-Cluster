use std::ops::Index;

#[derive(Clone, Debug)]
pub struct EnergyValues<T> {
    pub complet_energy: T,
    pub co_ads_energy: Option<T>,
}

#[derive(Clone, Debug)]
pub enum EnergyInput {
    LinearCn(EnergyValues<[i64; 2]>),
    Cn(EnergyValues<[i64; 13]>),
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
    }
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
