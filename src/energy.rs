// Pt
// const M_BETA: i64 = -0330;
// const M_ALPHA: i64 = 3960;

use core::panic;

// Co
const M_BETA: i64 = -0219;
const M_ALPHA: i64 = 2628;

// ciobacia
fn ciobacia_table(cn: usize) -> i64 {
    match cn {
        0 => 4427,
        1 => 3977,
        2 => 3576,
        3 => 3191,
        4 => 2816,
        5 => 2449,
        6 => 2087,
        7 => 1731,
        8 => 1378,
        9 => 1029,
        10 => 684,
        11 => 341,
        12 => 0,
        _ => panic!("impossible cn found"),
    }
}

pub fn energy_1000_calculation(atom: &u32, cn: &Vec<usize>) -> i64 {
    ciobacia_table(cn[*atom as usize])
    // cn[*atom as usize] as i64 * M_BETA + M_ALPHA
}

pub fn energy_diff(cn_from: usize, cn_to: usize) -> i64 {
    2 * ciobacia_table(cn_to) - 2 * ciobacia_table(cn_from)
    // (2 * ((cn_to as i64 - 1) * M_BETA + M_ALPHA)) - (2 * ((cn_from as i64) * M_BETA + M_ALPHA))
}
