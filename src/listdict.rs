use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use std::collections::HashMap;

fn pairing_function(a: u64, b: u64) -> u64 {
    // ((a + b) * (a + b + 1) / 2 + a) as usize
    // 2_u64.pow(a) * (2 * b + 1) - 1
    if a >= b {
        (a * a + a + b)
    } else {
        (a + b * b)
    }
}

#[derive(Clone)]
pub struct ListDict {
    move_to_position: HashMap<u64, usize, ahash::RandomState>,
    pub moves: Vec<(u32, u32)>,
}

impl ListDict {
    pub fn new(grid_size: [u32; 3]) -> ListDict {
        let largest_atom_position = grid_size[0] * grid_size[1] * grid_size[2] * 4;
        let item_to_position: HashMap<u64, usize, ahash::RandomState> = HashMap::default();
        // let item_to_position: HashMap<u64, usize, fnv::FnvBuildHasher> =
        //     fnv::FnvHashMap::with_capacity_and_hasher(32000, Default::default());
        ListDict {
            move_to_position: item_to_position,
            moves: Vec::with_capacity((largest_atom_position * 3) as usize),
        }
    }

    pub fn add_item(&mut self, move_from: u32, move_to: u32) {
        match self
            .move_to_position
            .entry((move_from as u64 + ((move_to as u64) << 32)))
        {
            std::collections::hash_map::Entry::Vacant(e) => {
                self.moves.push((move_from, move_to));
                e.insert(self.moves.len() - 1);
            }
            _ => return,
        }
    }
    pub fn remove_item(&mut self, move_from: u32, move_to: u32) {
        if let Some(position) = self
            .move_to_position
            .remove(&(move_from as u64 + ((move_to as u64) << 32)))
        {
            let (move_from, move_to) = self.moves.pop().unwrap();
            if position != self.moves.len() {
                self.moves[position] = (move_from, move_to);
                self.move_to_position
                    .insert((move_from as u64 + ((move_to as u64) << 32)), position);
            }
        }
    }

    pub fn choose_random_item(&self, rng_choose: &mut SmallRng) -> (u32, u32) {
        self.moves.choose(rng_choose).unwrap().clone()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, (u32, u32)> {
        self.moves.iter()
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

    pub fn _len(&self) -> usize {
        self.moves.len()
    }
}
