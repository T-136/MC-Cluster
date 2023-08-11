use fnv::FnvBuildHasher;
use fnv::FnvHashMap;
use std::collections::HashMap;
use std::collections::HashSet;

#[derive(Clone)]
pub struct EnergyChange {
    move_from_map: HashMap<u32, HashSet<u32>, FnvBuildHasher>,
    move_to_map: HashMap<u32, HashSet<u32>, FnvBuildHasher>,
    energy1000_change_map: HashMap<(u32, u32), i64, FnvBuildHasher>,
}

impl EnergyChange {
    pub fn new() -> EnergyChange {
        EnergyChange {
            move_from_map: FnvHashMap::default(),
            move_to_map: FnvHashMap::default(),
            energy1000_change_map: FnvHashMap::default(),
        }
    }

    pub fn add(&mut self, move_from: u32, move_to: u32, energy1000_diff: i64) {
        self.move_from_map
            .entry(move_from)
            .or_insert(HashSet::new())
            .insert(move_to);
        self.move_to_map
            .entry(move_to)
            .or_insert(HashSet::new())
            .insert(move_from);

        self.energy1000_change_map
            .insert((move_from, move_to), energy1000_diff);
    }

    pub fn get(&self, move_from: u32, move_to: u32) -> Option<&i64> {
        self.energy1000_change_map.get(&(move_from, move_to))
    }

    pub fn contains_move(&self, move_from: u32, move_to: u32) -> bool {
        self.energy1000_change_map
            .contains_key(&(move_from, move_to))
    }

    pub fn contains_position(&self, position: u32) -> Option<&HashSet<u32>> {
        if let Some(x) = self.move_from_map.get(&position) {
            return Some(x);
        } else {
            self.move_to_map.get(&position)
        }
    }

    pub fn remove(&mut self, position: u32) {
        if self.move_from_map.contains_key(&position) {
            let move_to_vec = self.move_from_map.remove(&position).unwrap();
            for move_to in move_to_vec.into_iter() {
                self.delete(position, move_to)
            }
        }
        if self.move_to_map.contains_key(&position) {
            let move_from_vec = self.move_to_map.remove(&position).unwrap();
            for move_from in move_from_vec.into_iter() {
                self.delete(position, move_from)
            }
        }
    }

    fn delete(&mut self, move_from: u32, move_to: u32) {
        self.energy1000_change_map.remove(&(move_from, move_to));
        if let Some(x) = self.move_from_map.get_mut(&move_from) {
            x.remove(&move_to);
        }
        if let Some(x) = self.move_to_map.get_mut(&move_to) {
            x.remove(&move_from);
        }
    }
}
