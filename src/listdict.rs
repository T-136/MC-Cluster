use fnv::FnvBuildHasher;
use fnv::FnvHashMap;
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use std::collections::HashMap;

fn pairing_function(a: u64, b: u64) -> usize {
    // ((a + b) * (a + b + 1) / 2 + a) as usize
    // 2_u64.pow(a) * (2 * b + 1) - 1
    if a >= b {
        (a * a + a + b) as usize
    } else {
        (a + b * b) as usize
    }
}

#[derive(Clone)]
pub struct ListDict {
    item_to_position: Vec<Option<usize>>,
    items: Vec<(u32, u32)>,
}

impl ListDict {
    pub fn new() -> ListDict {
        let largest_atom_position = 15_u64.pow(3) * 4;
        // println!("max size: {}", largest_atom_position);
        let size_vec = pairing_function(largest_atom_position, largest_atom_position);
        // println!("max size: {}", size_vec);
        let item_to_position = vec![None; size_vec];
        ListDict {
            item_to_position,
            items: Vec::with_capacity(32000),
        }
    }

    pub fn add_item(&mut self, move_from: u32, move_to: u32) {
        match self
            .item_to_position
            .get_mut(pairing_function(move_from as u64, move_to as u64))
            .unwrap()
        {
            Some(_) => {
                return;
            }
            None => {
                self.items.push((move_from, move_to));
                self.item_to_position[pairing_function(move_from as u64, move_to as u64)] =
                    Some(self.items.len() - 1);
            }
        }
    }
    pub fn remove_item(&mut self, move_from: u32, move_to: u32) {
        if let Some(position) =
            self.item_to_position[pairing_function(move_from as u64, move_to as u64)]
        {
            let (move_from_new, move_to_new) = self.items.pop().unwrap();
            if position != self.items.len() {
                self.items[position] = (move_from_new, move_to_new);
                if let Some(e) = self
                    .item_to_position
                    .get_mut(pairing_function(move_from_new as u64, move_to_new as u64))
                {
                    *e = Some(position)
                };
            }
            self.item_to_position[pairing_function(move_from as u64, move_to as u64)] = None;
        }
    }

    pub fn choose_random_item(&self, rng_choose: &mut SmallRng) -> (u32, u32) {
        self.items.choose(rng_choose).unwrap().clone()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, (u32, u32)> {
        self.items.iter()
    }

    // pub fn contains(&self, move_from: u32, move_to: u32) -> bool {
    //     self.item_to_position
    //         .contains_key(&(move_from as u64 + ((move_to as u64) << 32)))
    // }

    // pub fn remove_by_index(&mut self, index: usize) {
    //     self.item_to_position.remove(&self.items.swap_remove(index));
    // }

    // pub fn drain_filter(&mut self, cn: &Vec<usize>, move_from: &u32, move_to: &u32) {
    //     let mut i = 0;
    //     while i < self.items.len() {
    //         let (o, u) = self.items[i];
    //         if (cn[u as usize] == 0 || &o == move_from || &u == move_to) {
    //             let (move_from, move_to) = self.items.remove(i);
    //             self.item_to_position
    //                 .remove(&(move_from as u64 + ((move_to as u64) << 32)));
    //         } else {
    //             i += 1;
    //         }
    //     }
    // }
    // pub fn filter(&self) {
    //     self.items.iter().filter()
    // }

    // pub fn iter_mut(self) -> std::vec::IntoIter<(u64, u64)> {
    //     self.items.into_iter()
    // }

    pub fn len(&self) -> usize {
        self.items.len()
    }
}
