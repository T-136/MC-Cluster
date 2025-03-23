use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use std::collections::HashMap;

#[derive(Clone)]
pub struct ListDict {
    move_to_position: HashMap<u64, usize, ahash::RandomState>,
    pub moves: Vec<(u32, u32, Option<i64>)>, // [(from, to, energy_change)]
    total_energy_change: Option<i64>,
}

impl ListDict {
    pub fn new(nsites: u32) -> ListDict {
        let item_to_position: HashMap<u64, usize, ahash::RandomState> = HashMap::default();
        ListDict {
            move_to_position: item_to_position,
            moves: Vec::with_capacity((nsites * 3) as usize),
            total_energy_change: None,
        }
    }

    pub fn add_item(&mut self, move_from: u32, move_to: u32, energy_change: Option<i64>) {
        if let std::collections::hash_map::Entry::Vacant(e) = self
            .move_to_position
            .entry(move_from as u64 + ((move_to as u64) << 32))
        {
            self.moves.push((move_from, move_to, energy_change));
            e.insert(self.moves.len() - 1);
        }
    }

    pub fn remove_item(&mut self, move_from: u32, move_to: u32) {
        if let Some(position) = self
            .move_to_position
            .remove(&(move_from as u64 + ((move_to as u64) << 32)))
        {
            let (move_from, move_to, energy_change) = self.moves.pop().unwrap();
            if position != self.moves.len() {
                self.moves[position] = (move_from, move_to, energy_change);
                self.move_to_position
                    .insert(move_from as u64 + ((move_to as u64) << 32), position);
            }
        }
    }

    pub fn choose_random_item_mc(&self, rng_choose: &mut SmallRng) -> (u32, u32, Option<i64>) {
        *self.moves.choose(rng_choose).unwrap()
    }

    pub fn _iter(&self) -> std::slice::Iter<'_, (u32, u32, Option<i64>)> {
        self.moves.iter()
    }

    pub fn _len(&self) -> usize {
        self.moves.len()
    }
}
