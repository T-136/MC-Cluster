// Pt
const M_BETA: i64 = -0330;
const M_ALPHA: i64 = 3960;

use core::panic;

use crate::EnergyInput;

// Co
// const M_BETA: i64 = -0219;
// const M_ALPHA: i64 = 2628;

// enrico

const ENRICO_TABLE: [i64; 13] = [
    2690, 1879, 1660, 1441, 1242, 1043, 814, 725, 546, 427, 0, 219, 0,
];
fn enrico_table(cn: usize) -> i64 {
    match cn {
        0 => 2690,
        1 => 1879,
        2 => 1660,
        3 => 1441,
        4 => 1242,
        5 => 1043,
        6 => 814,
        7 => 725,
        8 => 546,
        9 => 427,
        10 => 0,
        11 => 219,
        12 => 0,
        _ => panic!("impossible cn found: {}", cn),
    }
}

pub fn energy_1000_calculation(energy: &super::EnergyInput, cn: usize) -> i64 {
    match energy {
        EnergyInput::LinearCn(e) => e[0] * cn as i64 + e[1],
        EnergyInput::Cn(e) => e[cn],
        EnergyInput::LinearGcn(e) => e[0] * cn as i64 + e[1],
        EnergyInput::Gcn(e) => e[cn],
    }
    // enrico_table(cn[*atom as usize])
    // cn[*atom as usize] as i64 * M_BETA + M_ALPHA
}

pub fn energy_diff_cn<I, O>(
    energy: [i64; 13],
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
        energy_diff_1000 -= energy[cn_from];
        energy_diff_1000 += energy[cn_from - 1];
    }
    for cn_to in cn_to_list {
        energy_diff_1000 -= energy[cn_to];
        energy_diff_1000 += energy[cn_to + 1];
    }
    energy_diff_1000 -= energy[move_from_cn];
    energy_diff_1000 += energy[move_to_cn - 1];

    energy_diff_1000
}
pub fn energy_diff_gcn<I, O>(
    energy: [i64; 60],
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
        energy_diff_1000 -= energy[cn_from];
        energy_diff_1000 += energy[cn_from - 1];
    }
    for cn_to in cn_to_list {
        energy_diff_1000 -= energy[cn_to];
        energy_diff_1000 += energy[cn_to + 1];
    }
    energy_diff_1000 -= energy[move_from_cn];
    energy_diff_1000 += energy[move_to_cn - 1];

    energy_diff_1000
}

pub fn energy_diff_l_cn(energy: [i64; 2], cn_from: usize, cn_to: usize) -> i64 {
    let e = (2 * ((cn_to as i64 - 1) * energy[0] + energy[1]))
        - (2 * ((cn_from as i64) * energy[0] + energy[1]));
    // if e < 0 {
    //     println!("{}", e);
    // }
    e
}
