use fnv::FnvBuildHasher;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

fn create_support(
    atom_pos: &mut Vec<super::AtomPosition>,
    xsites_positions: &Vec<[f64; 3]>,
    support_indices: &[i32],
    // support: &Support,
    nn: &HashMap<u32, [u32; 12], FnvBuildHasher>,
    iclose: u32,
) {
    let center_of_mass: &[f64; 3] = &xsites_positions[iclose as usize];
    let mut nn_support = vec![0_u8; xsites_positions.len()];
    let mut support_vec = Vec::new();
    let refpos = xsites_positions[0];
    for (i, xyz) in xsites_positions.iter().enumerate() {
        let new_vec = elementwise_subtraction(xyz, center_of_mass);
        if 1e-7 >= dot_product(&new_vec, support_indices).abs() {
            atom_pos[i].occ = 2;
            support_vec.push(i as u32)
        };
    }
    let mut second_layer_fixpoint = 0;
    for neighbor in nn[&iclose] {
        if atom_pos[neighbor as usize].occ == 0 {
            second_layer_fixpoint = neighbor;
        }
    }
    for sup in support_vec.iter() {
        for neighbor in nn[sup] {
            if atom_pos[neighbor as usize].occ != 2 {
                atom_pos[neighbor as usize].nn_support = 1;
            }
        }
    }
}

pub fn create_input_cluster(
    atom_pos: &mut Vec<super::AtomPosition>,
    number_of_atoms: &u32,
    xsites_positions: &Vec<[f64; 3]>,
    nn: &HashMap<u32, [u32; 12], FnvBuildHasher>,
    nsites: u32,
    // support: Option<&Support>,
    support_indices: Option<&Vec<i32>>,
) -> HashSet<u32, FnvBuildHasher> {
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
    println!("nsites: {}", nsites);
    assert_eq!(xsites_positions.len(), nsites as usize);
    let mut onlyocc: HashSet<u32, FnvBuildHasher> =
        fnv::FnvHashSet::with_capacity_and_hasher(*number_of_atoms as usize, Default::default());

    // if let Some(supp) = support {
    //     supp.support_indices
    //         .map(|sup| create_support(atom_pos, xsites_positions, sup, nn, iclose));
    // }
    if let Some(supp) = support_indices {
        create_support(atom_pos, xsites_positions, supp, nn, iclose);
    }

    // occ.entry(&iclose).and_modify(1) = 1;
    if atom_pos[iclose as usize].occ == 0 {
        onlyocc.insert(iclose);
    } else {
        for neighbor in nn[&iclose] {
            if atom_pos[neighbor as usize].occ == 0 {
                onlyocc.insert(neighbor);
                break;
            }
        }
    }

    let mut ks = 1;
    assert!(number_of_atoms < &nsites);

    loop {
        let mut onlyocc_temp_storag: HashSet<u32> = HashSet::new();
        for site in onlyocc.iter() {
            for j in &nn[site] {
                if !onlyocc_temp_storag.contains(&j)
                    && !onlyocc.contains(&j)
                    && atom_pos[*j as usize].occ == 0
                {
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
        if &ks == number_of_atoms {
            break;
        }
    }

    for site in onlyocc.iter() {
        atom_pos[*site as usize].occ = 1_u8;
    }

    onlyocc
}

pub fn occ_onlyocc_from_xyz(
    atom_pos: &mut Vec<super::AtomPosition>,
    xyz: Arc<Vec<(String, [f64; 3])>>,
    nsites: u32,
    xsites_positions: &Vec<[f64; 3]>,
    atom_names: super::AtomNames,
) -> HashSet<u32, FnvBuildHasher> {
    let mut occ: Vec<u8> = Vec::with_capacity(nsites as usize);
    for _ in 0..nsites {
        occ.push(0 as u8);
    }
    let mut onlyocc: HashSet<u32, FnvBuildHasher> =
        fnv::FnvHashSet::with_capacity_and_hasher(xyz.len(), Default::default());

    for x in xyz.iter() {
        for site in 0..nsites {
            let dist = (x.1[0] - xsites_positions[site as usize][0]).powf(2.)
                + (x.1[1] - xsites_positions[site as usize][1]).powf(2.)
                + (x.1[2] - xsites_positions[site as usize][2]).powf(2.);
            if dist < 0.15 {
                if &x.0 == atom_names.atom.as_ref().unwrap() {
                    atom_pos[site as usize].occ = 1;
                    onlyocc.insert(site);
                } else if let Some(supp) = atom_names.support.as_ref() {
                    if supp == &x.0 {
                        atom_pos[site as usize].occ = 2;
                    }
                }
            }
        }
    }
    // for (i, o) in occ.iter().enumerate() {
    //     if *o == 1_u8 {
    //         if !onlyocc.contains(&(i as u32)) {
    //             println!("occ: {}", i);
    //         }
    //     }
    // }
    onlyocc
}

fn elementwise_subtraction(vec_a: &[f64], vec_b: &[f64]) -> Vec<f64> {
    vec_a.iter().zip(vec_b).map(|(a, b)| a - b).collect()
}

fn dot_product(vec_a: &[f64], vec_b: &[i32]) -> f64 {
    vec_a
        .iter()
        .zip(vec_b.iter())
        .map(|(x, y)| x * *y as f64)
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_pos_in_pane() {
        let center = vec![15., 15., 15.];
        let vec_a = vec![16., 14., 15.];
        let plane_vec = vec![1, 1, 1];
        let new_vec = elementwise_subtraction(&vec_a, &center);
        assert!(dot_product(&new_vec, &plane_vec) <= 1e-7);

        let vec_a = vec![17., 14., 14.];
        let plane_vec = vec![1, 1, 1];
        let new_vec = elementwise_subtraction(&vec_a, &center);
        assert!(dot_product(&new_vec, &plane_vec) <= 1e-7);
    }
}
