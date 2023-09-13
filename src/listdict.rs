use rand::rngs::SmallRng;
use rand::seq::SliceRandom;

fn pairing_function(a: u32, b: u32) -> usize {
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
    move_to_position: Vec<Option<usize>>,
    moves_list: Vec<(u32, u32)>,
}

impl ListDict {
    pub fn new(grid_size: [u32; 3]) -> ListDict {
        let largest_atom_position = grid_size[0] * grid_size[1] * grid_size[2] * 4;
        let size_vec = pairing_function(largest_atom_position, largest_atom_position);
        let item_to_position = vec![None; size_vec];
        ListDict {
            move_to_position: item_to_position,
            moves_list: Vec::with_capacity((largest_atom_position * 3) as usize),
        }
    }

    pub fn add_item(&mut self, move_from: u32, move_to: u32) {
        if self.move_to_position[pairing_function(move_from, move_to)].is_none() {
            self.moves_list.push((move_from, move_to));
            self.move_to_position[pairing_function(move_from, move_to)] =
                Some(self.moves_list.len() - 1);
        }
    }
    pub fn remove_item(&mut self, move_from: u32, move_to: u32) {
        if let Some(position) = self.move_to_position[pairing_function(move_from, move_to)] {
            let (move_from_new, move_to_new) = self.moves_list.pop().unwrap();
            if position != self.moves_list.len() {
                self.moves_list[position] = (move_from_new, move_to_new);
                if let Some(move_index) = self
                    .move_to_position
                    .get_mut(pairing_function(move_from_new, move_to_new))
                {
                    *move_index = Some(position)
                };
            }
            self.move_to_position[pairing_function(move_from, move_to)] = None;
        }
    }

    pub fn choose_random_item(&self, rng_choose: &mut SmallRng) -> (u32, u32) {
        self.moves_list.choose(rng_choose).unwrap().clone()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, (u32, u32)> {
        self.moves_list.iter()
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
        self.moves_list.len()
    }
}
