use chemfiles::{Atom, Frame, Trajectory, UnitCell};
use fnv::FnvBuildHasher;
use fnv::FnvHashMap;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs;
use std::io::{self, BufRead};

pub fn write_occ_as_xyz(
    // trajectory: &mut Trajectory,
    save_folder: String,
    onlyocc: HashSet<u32, fnv::FnvBuildHasher>,
    xsites_positions: &Vec<[f64; 3]>,
    unit_cell: &[f64; 3],
    atom_pos: &Vec<super::AtomPosition>,
) {
    let mut trajectory = Trajectory::open(save_folder.clone() + "/lowest_energy.xyz", 'w').unwrap();
    let mut xyz: Vec<[f64; 3]> = Vec::new();
    for (j, ii) in onlyocc.iter().enumerate() {
        xyz.insert(j, xsites_positions[*ii as usize]);
    }
    let mut frame = Frame::new();
    frame.set_cell(&UnitCell::new(unit_cell.clone()));

    for atom in xyz.into_iter() {
        frame.add_atom(&Atom::new("Pt"), [atom[0], atom[1], atom[2]], None);
    }
    for (i, atom) in atom_pos.iter().enumerate() {
        if atom.occ == 2 {
            frame.add_atom(&Atom::new("Al"), xsites_positions[i], None);
        }
    }

    trajectory
        .write(&frame)
        .unwrap_or_else(|x| eprintln!("{}", x));
}

fn fmt_scient(num: &str) -> f64 {
    let mut parts = num.split('e');

    let pre_num = parts.next().unwrap();
    let exp = parts.next().unwrap();

    let base: f64 = 10.;
    pre_num.parse::<f64>().unwrap() * base.powi(exp.parse::<i32>().unwrap())
}

pub fn read_atom_sites(input_file: &str, nsites: u32) -> Vec<[f64; 3]> {
    println!("reading atom_sites from: {}", input_file);
    let mut xsites_positions: Vec<[f64; 3]> = Vec::with_capacity(nsites as usize);
    let pairlist = fs::File::open(input_file).expect("Should have been able to read the file");
    let lines = io::BufReader::new(pairlist);

    for line in lines.lines() {
        let r = line.unwrap();
        let list: Vec<&str> = r.split_whitespace().clone().collect();
        let temp_str_vec: [&str; 3] = [list[0], list[1], list[2]];
        let temp_vec: [f64; 3] = temp_str_vec.map(|i| fmt_scient(i));
        xsites_positions.push(temp_vec);
    }
    xsites_positions
}

pub fn read_nn(pairlist_file: &str) -> HashMap<u32, [u32; super::CN], FnvBuildHasher> {
    println!("reading pairlists from: {}", pairlist_file);

    let pairlist = fs::File::open(pairlist_file).expect("Should have been able to read the file");

    let lines = io::BufReader::new(pairlist);
    let mut nn: HashMap<u32, [u32; super::CN], FnvBuildHasher> =
        FnvHashMap::with_capacity_and_hasher(5400, Default::default());

    for line in lines.lines() {
        let r = line.unwrap();
        let list: Vec<&str> = r.split_whitespace().clone().collect();
        let mut neighbors: [u32; super::CN] = [0; super::CN];
        let prime = list.first().clone();
        for (i, l) in list.iter().skip(1).enumerate() {
            neighbors[i] = l.parse::<u32>().unwrap()
        }
        nn.insert(prime.unwrap().parse::<u32>().unwrap(), neighbors);
    }
    nn
}

pub fn read_nn_pair_no_intersec(
    nn_pairlist_file: &str,
) -> HashMap<u64, [[u32; super::NN_PAIR_NO_INTERSEC_NUMBER]; 2], FnvBuildHasher> {
    let nn_pairlist =
        fs::File::open(nn_pairlist_file).expect("Should have been able to read the file");

    let lines = io::BufReader::new(nn_pairlist);

    let mut nn_pair: HashMap<u64, [[u32; super::NN_PAIR_NO_INTERSEC_NUMBER]; 2], FnvBuildHasher> =
        FnvHashMap::with_capacity_and_hasher(32000, Default::default());

    for line in lines.lines() {
        let r = line.unwrap();
        let test: Vec<&str> = r.split_whitespace().clone().collect();
        let site: u32 = std::cmp::min(
            test[0].parse::<u32>().unwrap(),
            test[1].parse::<u32>().unwrap(),
        );
        let j: u32 = std::cmp::max(
            test[0].parse::<u32>().unwrap(),
            test[1].parse::<u32>().unwrap(),
        );
        let mut neighbors: [[u32; super::NN_PAIR_NO_INTERSEC_NUMBER]; 2] =
            [[0; super::NN_PAIR_NO_INTERSEC_NUMBER]; 2];

        for (i, l) in test.iter().skip(2).enumerate() {
            if i < 7 {
                neighbors[0][i] = l.parse::<u32>().unwrap()
            } else if i < 14 {
                neighbors[1][i - 7] = l.parse::<u32>().unwrap()
            }
        }
        nn_pair
            // .entry()
            // .and_modify(|map| {
            //     map.insert(j, neighbors.clone());
            // })
            .insert(site as u64 + ((j as u64) << 32), neighbors);
        // println!("{:?}", line.unwrap());
    }

    return nn_pair;
}

pub fn read_nn_pairlists(
    nn_pairlist_file: &str,
) -> HashMap<u64, [u32; super::NN_PAIR_NUMBER], FnvBuildHasher> {
    let nn_pairlist =
        fs::File::open(nn_pairlist_file).expect("Should have been able to read the file");

    let lines = io::BufReader::new(nn_pairlist);

    let mut nn_pair: HashMap<u64, [u32; super::NN_PAIR_NUMBER], FnvBuildHasher> =
        FnvHashMap::with_capacity_and_hasher(32000, Default::default());

    for line in lines.lines() {
        let r = line.unwrap();
        let test: Vec<&str> = r.split_whitespace().clone().collect();
        let site: u32 = std::cmp::min(
            test[0].parse::<u32>().unwrap(),
            test[1].parse::<u32>().unwrap(),
        );
        let j: u32 = std::cmp::max(
            test[0].parse::<u32>().unwrap(),
            test[1].parse::<u32>().unwrap(),
        );
        let mut neighbors: [u32; super::NN_PAIR_NUMBER] = [0; super::NN_PAIR_NUMBER];

        for (i, l) in test.iter().skip(2).enumerate() {
            neighbors[i] = l.parse::<u32>().unwrap()
        }
        nn_pair
            // .entry()
            // .and_modify(|map| {
            //     map.insert(j, neighbors.clone());
            // })
            .insert(site as u64 + ((j as u64) << 32), neighbors);
        // println!("{:?}", line.unwrap());
    }

    return nn_pair;
}
