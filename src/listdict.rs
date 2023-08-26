use fnv::FnvBuildHasher;
use fnv::FnvHashMap;
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use std::collections::HashMap;

#[derive(Clone)]
pub struct ListDict {
    item_to_position: HashMap<u64, usize, FnvBuildHasher>,
    items: Vec<(u32, u32)>,
}

impl ListDict {
    pub fn new() -> ListDict {
        let item_to_position: HashMap<u64, usize, FnvBuildHasher> =
            FnvHashMap::with_capacity_and_hasher(32000, Default::default());
        ListDict {
            item_to_position,
            items: Vec::with_capacity(32000),
        }
    }

    pub fn add_item(&mut self, move_from: u32, move_to: u32) {
        if let std::collections::hash_map::Entry::Vacant(e) = self
            .item_to_position
            .entry(move_from as u64 + ((move_to as u64) << 32))
        {
            self.items.push((move_from, move_to));
            e.insert(self.items.len() - 1);
        } else {
            return;
        }
    }
    pub fn remove_item(&mut self, move_from: u32, move_to: u32) {
        if let Some(position) = self
            .item_to_position
            .remove(&(move_from as u64 + ((move_to as u64) << 32)))
        {
            let (move_from, move_to) = self.items.pop().unwrap();
            if position != self.items.len() {
                self.items[position] = (move_from, move_to);
                self.item_to_position
                    .insert(move_from as u64 + ((move_to as u64) << 32), position);
            }
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

    // pub fn len(&self) -> usize {
    //     self.items.len()
    // }
}
