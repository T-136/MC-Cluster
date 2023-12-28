use chemfiles::Frame;
use chemfiles::Trajectory;
use fnv::FnvBuildHasher;
use fnv::FnvHashMap;
use std::collections::HashMap;
use std::fs;
use std::io::{self, BufRead};
use time_graph::instrument;
use vasp_poscar::Poscar;

#[instrument]
pub fn read_sample(input_file: &str) -> Vec<[f64; 3]> {
    if input_file.contains(".poscar") {
        let newatoms = Poscar::from_path(input_file).unwrap();
        let xyz = newatoms.scaled_cart_positions();
        xyz.into_owned()
    } else if input_file.contains(".xyz") {
        let mut trajectory = Trajectory::open(input_file, 'r').unwrap();
        let mut frame = Frame::new();
        trajectory.read(&mut frame).unwrap();
        frame.positions().to_owned()
    } else {
        panic!("no .poscar or .xyz, cant read file");
    }
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
pub fn read_nnn(pairlist_file: &str) -> HashMap<u32, [u32; super::GCN], FnvBuildHasher> {
    println!("reading pairlists from: {}", pairlist_file);

    let pairlist = fs::File::open(pairlist_file).expect("Should have been able to read the file");

    let lines = io::BufReader::new(pairlist);
    let mut nnn: HashMap<u32, [u32; super::GCN], FnvBuildHasher> =
        FnvHashMap::with_capacity_and_hasher(5400, Default::default());

    for line in lines.lines() {
        let r = line.unwrap();
        let list: Vec<&str> = r.split_whitespace().clone().collect();
        let mut neighbors: [u32; super::GCN] = [0; super::GCN];
        let prime = list.first().clone();
        for (i, l) in list.iter().skip(1).enumerate() {
            neighbors[i] = l.parse::<u32>().unwrap()
        }
        nnn.insert(prime.unwrap().parse::<u32>().unwrap(), neighbors);
    }
    nnn
}
pub fn read_nnn_pair_no_intersec(
    nnn_pairlist_file: &str,
) -> HashMap<
    u64,
    (
        Vec<Vec<u32>>,
        Vec<Vec<u32>>,
        Vec<(u32, Vec<u32>, Vec<u32>, Vec<u32>)>,
    ),
    fnv::FnvBuildHasher,
> {
    let nnn_pairlist =
        fs::File::open(nnn_pairlist_file).expect("Should have been able to read the file");

    let reader = io::BufReader::new(nnn_pairlist);

    let nnn_pair_no_bit_shifting: HashMap<
        u32,
        HashMap<
            u32,
            (
                Vec<Vec<u32>>,
                Vec<Vec<u32>>,
                Vec<(u32, Vec<u32>, Vec<u32>, Vec<u32>)>,
            ),
            fnv::FnvBuildHasher,
        >,
    > = serde_json::from_reader(reader).unwrap();

    let mut nnn_pair: HashMap<
        u64,
        (
            Vec<Vec<u32>>,
            Vec<Vec<u32>>,
            Vec<(u32, Vec<u32>, Vec<u32>, Vec<u32>)>,
        ),
        fnv::FnvBuildHasher,
    > = FnvHashMap::with_capacity_and_hasher(5400, Default::default());
    for (k1, d1) in nnn_pair_no_bit_shifting.into_iter() {
        for (k2, d2) in d1.into_iter() {
            nnn_pair.insert((k1 as u64 + ((k2 as u64) << 32)), d2);
        }
    }

    // let nnn_pair: HashMap<u64, [HashMap<u32, Vec<u32>, FnvBuildHasher>; 2], FnvBuildHasher> =
    //     serde_json::from_reader(reader).unwrap();

    return nnn_pair;
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
            } else {
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

pub fn read_nnn_pairlists(nn_pairlist_file: &str) -> HashMap<u64, [u32; 74], FnvBuildHasher> {
    let nn_pairlist =
        fs::File::open(nn_pairlist_file).expect("Should have been able to read the file");

    let lines = io::BufReader::new(nn_pairlist);

    let mut nn_pair: HashMap<u64, [u32; 74], FnvBuildHasher> =
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
        let mut neighbors: [u32; 74] = [0; 74];

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
