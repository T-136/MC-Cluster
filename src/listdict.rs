use boomphf::hashmap::BoomHashMap;
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
    move_to_position: BoomHashMap<u64, Option<usize>>,
    moves: Vec<(u32, u32)>,
}

impl ListDict {
    pub fn new(grid_size: [u32; 3], nn_pairlist: Vec<u64>) -> ListDict {
        let largest_atom_position = grid_size[0] * grid_size[1] * grid_size[2] * 4;
        let nn_pairlist_balues: Vec<Option<usize>> = vec![None; nn_pairlist.len()];

        // let phf = Mphf::new(1.7, &nn_pairlist);
        let item_to_position: BoomHashMap<u64, Option<usize>> =
            BoomHashMap::new(nn_pairlist, nn_pairlist_balues);

        // let item_to_position: HashMap<u64, usize, fnv::FnvBuildHasher> =
        //     fnv::FnvHashMap::with_capacity_and_hasher(32000, Default::default());
        ListDict {
            move_to_position: item_to_position,
            moves: Vec::with_capacity((largest_atom_position * 3) as usize),
        }
    }

    pub fn add_item(&mut self, move_from: u32, move_to: u32) {
        let x = self
            .move_to_position
            .get_mut(&(move_from as u64 + ((move_to as u64) << 32)))
            .unwrap();
        if x.is_none() {
            self.moves.push((move_from, move_to));
            *x = Some(self.moves.len() - 1);
        }
    }

    pub fn remove_item(&mut self, move_from: u32, move_to: u32) {
        let new_position = match self
            .move_to_position
            .get(&(move_from as u64 + ((move_to as u64) << 32)))
            .unwrap()
        {
            Some(position) => {
                let (new_move_from, new_move_to) = self.moves.pop().unwrap();
                if *position != self.moves.len() {
                    if *position > self.moves.len() {
                        println!("{}", self.moves.len());
                        println!(
                            "{:?}",
                            self.moves.iter().position(|&x| x == ((move_from, move_to)))
                        );
                    }
                    self.moves[*position] = (new_move_from, new_move_to);
                    Some((*position, new_move_from, new_move_to))
                    // self.move_to_position
                    //     .get_mut(&(move_from as u64 + ((move_to as u64) << 32)))
                    //     .map(|val| *val = Some(*position));
                } else {
                    None
                }
            }
            _ => return,
        };
        let val = self
            .move_to_position
            .get_mut(&(move_from as u64 + ((move_to as u64) << 32)))
            .unwrap();
        {
            *val = None
        }
        if let Some((pos, new_move_from, new_move_to)) = new_position {
            if let Some(val) = self
                .move_to_position
                .get_mut(&(new_move_from as u64 + ((new_move_to as u64) << 32)))
                .unwrap()
            {
                *val = pos
            }
        }

        // println!(
        //     "{:?}",
        //     self.move_to_position
        //         .get(&(move_from as u64 + ((move_to as u64) << 32)))
        // );
        // println!(
        //     "{:?}",
        //     self.moves.iter().position(|&x| x == ((move_from, move_to)))
        // );
    }

    pub fn count_hashmap(&self) -> u32 {
        let mut count = 0;
        // println!("{:?}", self.move_to_position);
        for (x, v) in self.move_to_position.iter() {
            // println!("n1 {:?}", v);
            if v.is_some() {
                count += 1;
            }
        }
        count
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
