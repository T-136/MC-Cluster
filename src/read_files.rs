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

pub fn read_nn(pairlist_file: &str) -> HashMap<u32, [u32; 12], FnvBuildHasher> {
    println!("reading pairlists from: {}", pairlist_file);

    let pairlist = fs::File::open(pairlist_file).expect("Should have been able to read the file");

    let lines = io::BufReader::new(pairlist);
    let mut nn: HashMap<u32, [u32; 12], FnvBuildHasher> =
        HashMap::with_capacity_and_hasher(5400, Default::default());

    for line in lines.lines() {
        let r = line.unwrap();
        let list: Vec<&str> = r.split_whitespace().clone().collect();
        let mut neighbors: [u32; 12] = [0; 12];
        let prime = list.first().clone();
        for (i, l) in list.iter().skip(1).enumerate() {
            neighbors[i] = l.parse::<u32>().unwrap()
        }
        nn.insert(prime.unwrap().parse::<u32>().unwrap(), neighbors);
        // println!("{:?}", line.unwrap());
    }
    nn
}

pub fn read_nn_pairlists(nn_pairlist_file: &str) -> HashMap<u64, [u32; 20], FnvBuildHasher> {
    let nn_pairlist =
        fs::File::open(nn_pairlist_file).expect("Should have been able to read the file");

    let lines = io::BufReader::new(nn_pairlist);

    let mut nn_pair: HashMap<u64, [u32; 20], FnvBuildHasher> =
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
        let mut neighbors: [u32; 20] = [0; 20];

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
