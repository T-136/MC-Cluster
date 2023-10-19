// Pt
// const M_BETA: i64 = -0330;
// const M_ALPHA: i64 = 3960;

// Co
const M_BETA: i64 = -0219;
const M_ALPHA: i64 = 2628;

pub fn energy_calculation(atom: &u32, cn: &Vec<usize>) -> i64 {
    let energy_1000 = cn[*atom as usize] as i64 * M_BETA + M_ALPHA;
    energy_1000
}

pub fn energy_diff(cn_from: usize, cn_to: usize) -> i64 {
    (2 * ((cn_to as i64 - 1) * M_BETA + M_ALPHA)) - (2 * ((cn_from as i64) * M_BETA + M_ALPHA))
}
