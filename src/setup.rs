use fnv::FnvBuildHasher;
use std::collections::{HashMap, HashSet};

pub fn create_input_cluster(
    number_of_atoms: &u32,
    xsites_positions: &Vec<[f64; 3]>,
    nn: &HashMap<u32, [u32; 12], FnvBuildHasher>,
    nsites: u32,
) -> (Vec<u8>, HashSet<u32>) {
    let center_of_mass: [f64; 3] = {
        let mut d: [Vec<f64>; 3] = [Vec::new(), Vec::new(), Vec::new()];
        for coord in xsites_positions {
            d[0].push(coord[0]);
            d[1].push(coord[1]);
            d[2].push(coord[2]);
        }
        [
            d[0].iter().sum::<f64>() / d[0].len() as f64,
            d[1].iter().sum::<f64>() / d[1].len() as f64,
            d[2].iter().sum::<f64>() / d[2].len() as f64,
        ]
    };
    let iclose: u32 = {
        let mut minimum: f64 = f64::INFINITY;
        let mut index_of_center: u32 = 0;
        // let mut v = Vec::new();
        for (i, xsite_position) in xsites_positions.iter().enumerate() {
            let dist = (center_of_mass[0] - xsite_position[0]).powf(2.)
                + (center_of_mass[1] - xsite_position[1]).powf(2.)
                + (center_of_mass[2] - xsite_position[2]).powf(2.);
            if minimum >= dist {
                minimum = dist.clone();
                index_of_center = i as u32;
            }
        }
        index_of_center
    };
    println!(
        "{:?}, {:?}, {:?}",
        center_of_mass, iclose, xsites_positions[iclose as usize]
    );
    let mut occ: Vec<u8> = Vec::with_capacity(xsites_positions.len());
    let mut onlyocc: HashSet<u32> = HashSet::with_capacity(*number_of_atoms as usize);

    for _ in 0..nsites {
        occ.push(0);
    }

    // occ.entry(&iclose).and_modify(1) = 1;
    onlyocc.insert(iclose);

    let mut ks = 1;

    loop {
        let mut onlyocc_temp_storag: HashSet<u32> = HashSet::new();
        for site in onlyocc.iter() {
            for j in &nn[site] {
                if !onlyocc_temp_storag.contains(&j) && !onlyocc.contains(&j) {
                    onlyocc_temp_storag.insert(*j);
                    // onlyocc.insert(j);
                    ks += 1;
                    if &ks == number_of_atoms {
                        break;
                    }
                }
            }
            if &ks == number_of_atoms {
                break;
            }
        }
        for v in onlyocc_temp_storag {
            onlyocc.insert(v);
        }
        println!("{}", onlyocc.len());
        if &ks == number_of_atoms {
            break;
        }
    }
    println!("{}", onlyocc.len());

    for site in onlyocc.iter() {
        occ[*site as usize] = 1_u8;
    }

    (occ, onlyocc)
    // build_onlyocc = [iclose]
    // loc_occ =[]
    // number_atoms = int(number_atoms)
    // ks = len(build_onlyocc)
    // while True:
    //     loc_occ = copy.deepcopy(build_onlyocc)
    //     for site in loc_occ:
    //         for j in nn[site]:
    //             if j not in build_onlyocc:
    //                 build_onlyocc.append(j)
    //                 ks += 1
    //                 if ks == number_atoms:
    //                     break
    //         if ks == number_atoms:
    //             break
    //     # build_onlyocc = loc_occ.copy()
    //     if ks == number_atoms:
    //         break
    // todel = [i for i in sites if i not in build_onlyocc]
    // base_traj = atoms.copy()
    // del base_traj[todel]
    // base_traj.write('./../input_cluster/test.poscar')
}

pub fn occ_onlyocc_from_xyz(
    xyz: &Vec<[f64; 3]>,
    nsites: u32,
    xsites_positions: &Vec<[f64; 3]>,
) -> (Vec<u8>, HashSet<u32>) {
    let mut occ: Vec<u8> = Vec::with_capacity(nsites as usize);
    for _ in 0..nsites {
        occ.push(0 as u8);
    }
    let mut onlyocc: HashSet<u32> = HashSet::with_capacity(xyz.len());

    for x in xyz.iter() {
        for site in 0..nsites {
            let dist = (x[0] - xsites_positions[site as usize][0]).powf(2.)
                + (x[1] - xsites_positions[site as usize][1]).powf(2.)
                + (x[2] - xsites_positions[site as usize][2]).powf(2.);
            if dist < 0.15 {
                occ[site as usize] = 1;
                onlyocc.insert(site);
            }
        }
    }
    for (i, o) in occ.iter().enumerate() {
        if *o == 1_u8 {
            if !onlyocc.contains(&(i as u32)) {
                println!("{}", i);
            }
        }
    }
    (occ, onlyocc)
}
