// Pt
const M_BETA: i64 = -330;
const M_ALPHA: i64 = 3960;

// const support_e: i64 = -211;

#[derive(Clone, Debug)]
pub enum EnergyInput {
    LinearCn([i64; 2]),
    Cn([i64; 13]),
    LinearGcn([i64; 2]),
    Gcn([i64; 145]),
}

pub fn energy_1000_calculation(
    energy: &EnergyInput,
    cn: usize,
    at_support: u8,
    support_e: i64,
) -> i64 {
    match energy {
        EnergyInput::LinearCn(e) => e[0] * cn as i64 + e[1] + support_e * at_support as i64,
        EnergyInput::Cn(e) => e[cn] + support_e * at_support as i64,
        EnergyInput::LinearGcn(e) => {
            println!("e1 {}", e[0]);
            println!("e2 {}", e[1]);
            println!("cn {}", cn);
            e[0] * cn as i64 + e[1]
        }
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
    from_at_support: u8,
    to_at_support: u8,
    support_e: i64,
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
    energy_diff_1000 -= support_e * from_at_support as i64;
    energy_diff_1000 += support_e * to_at_support as i64;

    energy_diff_1000
}
pub fn energy_diff_gcn<I, O, P>(
    energy: [i64; 145],
    cn_from_list: I,
    cn_to_list: O,
    intersect_change: P,
    move_from_gcn: usize,
    move_to_gcn: usize,
) -> i64
where
    I: Iterator<Item = (usize, usize)>,
    O: Iterator<Item = (usize, usize)>,
    P: Iterator<Item = (usize, usize)>,
{
    let mut energy_diff_1000 = 0;
    for (gcn_from_old, gcn_from_new) in cn_from_list {
        energy_diff_1000 -= energy[gcn_from_old];
        energy_diff_1000 += energy[gcn_from_new];
    }
    for (gcn_to_old, gcn_to_new) in cn_to_list {
        energy_diff_1000 -= energy[gcn_to_old];
        energy_diff_1000 += energy[gcn_to_new];
    }
    for (gcn_inter_old, gcn_inter_new) in intersect_change {
        energy_diff_1000 -= energy[gcn_inter_old];
        energy_diff_1000 += energy[gcn_inter_new];
    }
    energy_diff_1000 += energy[move_to_gcn];
    energy_diff_1000 -= energy[move_from_gcn];

    energy_diff_1000
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
