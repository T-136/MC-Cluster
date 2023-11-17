// Pt
const M_BETA: i64 = -0330;
const M_ALPHA: i64 = 3960;

use core::panic;

// Co
// const M_BETA: i64 = -0219;
// const M_ALPHA: i64 = 2628;

// enrico
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
        _ => panic!("impossible cn found"),
    }
}

pub fn energy_1000_calculation(atom: &u32, cn: &Vec<usize>) -> i64 {
    // enrico_table(cn[*atom as usize])
    cn[*atom as usize] as i64 * M_BETA + M_ALPHA
}

pub fn energy_diff(cn_from: usize, cn_to: usize) -> i64 {
    // 2 * enrico_table(cn_to) - 2 * enrico_table(cn_from)
    (2 * ((cn_to as i64 - 1) * M_BETA + M_ALPHA)) - (2 * ((cn_from as i64) * M_BETA + M_ALPHA))
}
