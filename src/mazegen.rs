use rand::Rng;
use rand::rngs::ThreadRng;
use rand::seq::IteratorRandom;
use std::collections::HashSet;

fn calculate_tiles_distance(tilei: [usize; 2], tilef: [usize; 2]) -> f32 {
    let dx = if tilef[0] > tilei[0] {
        tilef[0] - tilei[0]
    } else {
        tilei[0] - tilef[0]
    };
    let dy = if tilef[1] > tilei[1] {
        tilef[1] - tilei[1]
    } else {
        tilei[1] - tilef[1]
    };
    f32::sqrt((dx * dx + dy * dy) as f32)
}

fn select_random_tile_from_set(set: &Vec<[usize; 2]>, rng: &mut ThreadRng) -> [usize; 2] {
    set[rng.gen_range(0..set.len())]
}

fn select_random_tiles_from_set(set: &Vec<[usize; 2]>, weight: f32, rng: &mut ThreadRng) -> Vec<[usize; 2]> {
    let amount = (weight * set.len() as f32) as usize;
    set.iter().copied().choose_multiple(rng, amount)
}

enum Iter<I: Iterator> {
    Standard(I),
    Reverse(std::iter::Rev<I>),
}

impl<I: Iterator + DoubleEndedIterator> Iterator for Iter<I> {
    type Item = I::Item;
    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Iter::Standard(iter) => iter.next(),
            Iter::Reverse(iter) => iter.next(),
        }
    }
}

pub struct MazeManager {
    pub grid: Vec<Vec<bool>>,
    pub true_tiles: HashSet<[usize; 2]>,
    pub false_tiles: HashSet<[usize; 2]>,
    pub width: usize,
    pub height: usize,
    center_x: usize,
    center_y: usize,
    area: f32
}

impl MazeManager {
    pub fn instantiate_map(width: usize, height: usize, fill: bool) -> MazeManager {
        let mut true_tiles = HashSet::new();
        let mut false_tiles = HashSet::new();
        if fill {
            for x in 0..width {
                for y in 0..height {
                    true_tiles.insert([x, y]);
                }
            }
        } else {
            for x in 0..width {
                for y in 0..height {
                    false_tiles.insert([x, y]);
                }
            }
        }
        MazeManager {
            grid: vec![vec![fill; width]; height],
            width,
            height,
            center_x: f32::ceil(width as f32 * 0.5) as usize,
            center_y: f32::ceil(height as f32 * 0.5) as usize,
            area: (width * height) as f32,
            true_tiles,
            false_tiles
        }
    }

    fn calculate_directions(&self, x: usize, y: usize, direction_type: bool) -> Vec<[usize; 2]> {
        let directions = if direction_type {
            vec![[-1, -1], [-1, 0], [-1, 1], [0, -1], [0, 1], [1, -1], [1, 0], [1, 1]]
        } else {
            vec![[-1, 0], [0, -1], [0, 1], [1, 0]]
        };
        let mut directions_result = Vec::new();
        for [dx, dy] in directions {
            let new_x = x as isize + dx;
            let new_y = y as isize + dy;
            if new_x >= 0 && new_x < self.width as isize && new_y >= 0 && new_y < self.height as isize {
                directions_result.push([new_x as usize, new_y as usize]);
            }
        }
        directions_result
    }

    pub fn set_tile(&mut self, x: usize, y: usize, set_to: bool) {
        if self.grid[y][x] != set_to {
            if !set_to {
                self.true_tiles.remove(&[x, y]);
                self.false_tiles.insert([x, y]);
                self.grid[y][x] = false;
            } else {
                self.false_tiles.remove(&[x, y]);
                self.true_tiles.insert([x, y]);
                self.grid[y][x] = true;
            }
        }
    }
    
    fn list_set_tile(&mut self, x: usize, y: usize, set_to: bool, completed_coords: &mut Vec<[usize; 2]>) {
        self.set_tile(x, y, set_to);
        completed_coords.push([x, y]);
    }

    fn check_surrounding_tiles(&self, direction_type: bool, compare: bool, x: usize, y: usize) -> usize {
        let mut flag = 0;
        let directions = self.calculate_directions(x, y, direction_type);
        let max_directions_len = if direction_type {
            8
        } else { 
            4
        };
        if directions.len() != max_directions_len {
            flag += max_directions_len - directions.len();
        }
        for direction in directions {
            if self.grid[direction[1]][direction[0]] == compare {
                flag += 1;
            }
        }
        flag
    }

    fn select_random_tile(&self, state: bool, rng: &mut ThreadRng) -> [usize; 2] {
        if state {
            return self.true_tiles.iter().copied().choose(rng).unwrap();
        } else {
            return self.false_tiles.iter().copied().choose(rng).unwrap();
        }
    }

    fn select_random_tile_from_distance(&self, state: bool, coord: [usize; 2], min_distance: f32, rng: &mut ThreadRng) -> [usize; 2] {
        let mut list = if state {
            self.true_tiles.iter().copied().collect::<Vec<[usize; 2]>>()
        } else {
            self.false_tiles.iter().copied().collect::<Vec<[usize; 2]>>()
        };
        list.retain(|&coordinate| {
            let distance = calculate_tiles_distance(coord, coordinate);
            distance >= min_distance
        });
        if list.len() != 0 {
            return list[rng.gen_range(0..list.len())];
        } else {
            return self.select_random_tile_from_distance(state, coord, min_distance - 1.0, rng);
        }
    }

    fn select_random_tile_pair(&self, state: bool, rng: &mut ThreadRng) -> [[usize; 2]; 2] {
        if state {
            let tiles = self.true_tiles.iter().copied().choose_multiple(rng, 2);
            return [tiles[0], tiles[1]];
        } else {
            let tiles = self.false_tiles.iter().copied().choose_multiple(rng, 2);
            return [tiles[0], tiles[1]];
        }
    }

    fn select_random_tile_pair_from_distance(&self, state: bool, min_distance: f32, rng: &mut ThreadRng) -> [[usize; 2]; 2] {
        let list = if state {
            self.true_tiles.iter().copied().collect::<Vec<[usize; 2]>>()
        } else {
            self.false_tiles.iter().copied().collect::<Vec<[usize; 2]>>()
        };
        let mut compare_list = Vec::new();
        let mut result_list = Vec::new();
        for entry in list {
            for compare in &compare_list {
                if calculate_tiles_distance(entry, *compare) >= min_distance {
                    result_list.push([entry, *compare]);
                }
            }
            compare_list.push(entry);
        }
        if result_list.len() != 0 {
            return result_list[rng.gen_range(0..result_list.len())];
        } else {
            return self.select_random_tile_pair_from_distance(state, min_distance - 1.0, rng);
        }
    }

    fn select_random_tiles(&self, state: bool, weight: f32, rng: &mut ThreadRng) -> Vec<[usize; 2]> {
        if state {
            let amount = (weight * self.true_tiles.len() as f32) as usize;
            return self.true_tiles.iter().copied().choose_multiple(rng, amount);
        } else {
            let amount = (weight * self.false_tiles.len() as f32) as usize;
            return self.false_tiles.iter().copied().choose_multiple(rng, amount);
        }
    }

    fn select_random_tile_pairs(&self, state: bool, weight: f32, rng: &mut ThreadRng) -> Vec<[[usize; 2]; 2]> {
        let tiles = if state {
            let mut amount = (weight * self.true_tiles.len() as f32) as usize;
            if amount % 2 == 1 {
                amount -= 1;
            }
            self.true_tiles.iter().copied().choose_multiple(rng, amount)
        } else {
            let mut amount = (weight * self.false_tiles.len() as f32) as usize;
            if amount % 2 == 1 {
                amount -= 1;
            }
            self.false_tiles.iter().copied().choose_multiple(rng, amount)
        };
        let mut pairs = Vec::new();
        for tile in (0..tiles.len()).step_by(2) {
            pairs.push([tiles[tile], tiles[tile + 1]]);
        }
        pairs
    }
    
    fn create_horizontal_branch(&mut self, fill_first_tile: bool, set_to: bool, direction: bool, x: usize, y: usize, length: usize, sway_up_weight: f32, sway_down_weight: f32, rng: &mut ThreadRng) -> Vec<[usize; 2]> {
        let straight_path_weight = sway_up_weight + sway_down_weight;
        let mut y_coord = y;
        let mut completed_coords = Vec::new();
        let iter = if direction {
            let x_start = if fill_first_tile {
                x
            } else {
                completed_coords.push([x, y]);
                x + 1
            };
            Iter::Standard(x_start..usize::min(x + length, self.width))
        } else {
            let x_end = if fill_first_tile {
                x + 1
            } else {
                completed_coords.push([x, y]);
                x
            };
            if x < length + 1 {
                Iter::Reverse((0..x_end).rev())
            } else {
                Iter::Reverse((x - length - 1..x_end).rev())
            }
        };
        for x_coord in iter {
            self.list_set_tile(x_coord, y_coord, set_to, &mut completed_coords);
            let sway_chance = rng.gen_range(0.0..1.0);
            if sway_chance < sway_up_weight {
                if y_coord != 0 {
                    y_coord -= 1;
                    self.list_set_tile(x_coord, y_coord, set_to, &mut completed_coords);
                }
            } else if sway_chance < straight_path_weight {
                if y_coord != self.height - 1 {
                    y_coord += 1;
                    self.list_set_tile(x_coord, y_coord, set_to, &mut completed_coords);
                }
            }
        }
        completed_coords
    }
    
    fn create_vertical_branch(&mut self, fill_first_tile: bool, set_to: bool, direction: bool, x: usize, y: usize, length: usize, sway_left_weight: f32, sway_right_weight: f32, rng: &mut ThreadRng) -> Vec<[usize; 2]> {
        let straight_path_weight = sway_left_weight + sway_right_weight;
        let mut x_coord = x;
        let mut completed_coords = Vec::new();
        let iter = if direction {
            let y_start = if fill_first_tile {
                y
            } else {
                completed_coords.push([x, y]);
                y + 1
            };
            Iter::Standard(y_start..usize::min(y + length, self.height))
        } else {
            let y_end = if fill_first_tile {
                y + 1
            } else {
                completed_coords.push([x, y]);
                y
            };
            if y < length + 1 {
                Iter::Reverse((0..y_end).rev())
            } else {
                Iter::Reverse((y - length - 1..y_end).rev())
            }
        };
        for y_coord in iter {
            self.list_set_tile(x_coord, y_coord, set_to, &mut completed_coords);
            let sway_chance = rng.gen_range(0.0..1.0);
            if sway_chance < sway_left_weight {
                if x_coord != 0 {
                    x_coord -= 1;
                    self.list_set_tile(x_coord, y_coord, set_to, &mut completed_coords);
                }
            } else if sway_chance < straight_path_weight {
                if x_coord != self.width - 1 {
                    x_coord += 1;
                    self.list_set_tile(x_coord, y_coord, set_to, &mut completed_coords);
                }
            }
        }
        completed_coords
    }
    
    fn create_tunnel(&mut self, fill_first_tile: bool, set_to: bool, coordi: [usize; 2], coordf: [usize; 2], rng: &mut ThreadRng) -> Vec<[usize; 2]> {
        let mut completed_coords = Vec::new();
        let mut cur_coord = coordi;
        self.list_set_tile(coordf[0], coordf[1], set_to, &mut completed_coords);
        if fill_first_tile {
            self.list_set_tile(cur_coord[0], cur_coord[1], set_to, &mut completed_coords);
        }
        while cur_coord[0] != coordf[0] || cur_coord[1] != coordf[1] {
            if cur_coord[0] != coordf[0] && cur_coord[1] != coordf[1] {
                if rng.gen_range(0.0..1.0) < 0.5 {
                    if cur_coord[0] < coordf[0] {
                        cur_coord[0] += 1;
                    } else {
                        cur_coord[0] -= 1;
                    }
                } else {
                    if cur_coord[1] < coordf[1] {
                        cur_coord[1] += 1;
                    } else {
                        cur_coord[1] -= 1;
                    }
                }
            } else {
                if cur_coord[0] < coordf[0] {
                    cur_coord[0] += 1;
                } else if cur_coord[0] > coordf[0] {
                    cur_coord[0] -= 1;
                } else if cur_coord[1] < coordf[1] {
                    cur_coord[1] += 1;
                } else {
                    cur_coord[1] -= 1;
                }
            }
            self.list_set_tile(cur_coord[0], cur_coord[1], set_to, &mut completed_coords);
        }
        completed_coords
    }

    fn create_random_tunnel(&mut self, set_to: bool, coord: [usize; 2], weight: f32, up_factor: usize, down_factor: usize, left_factor: usize, right_factor: usize, rng: &mut ThreadRng) -> Vec<[usize; 2]> {
        let mut completed_coords = Vec::new();
        let mut cur_coord = coord;
        let mut coordinate_set = HashSet::new();
        while if set_to { self.true_tiles.len() as f32 } else { self.false_tiles.len() as f32 } / self.area < weight {
            self.list_set_tile(cur_coord[0], cur_coord[1], set_to, &mut completed_coords);
            coordinate_set.insert([cur_coord[0], cur_coord[1]]);
            let mut options = Vec::new();
            let mut backup_options = Vec::new();
            if cur_coord[0] != 0 {
                if !coordinate_set.contains(&[cur_coord[0] - 1, cur_coord[1]]) {
                    for _count in 0..left_factor {
                        options.push([cur_coord[0] - 1, cur_coord[1]]);
                    }
                } else {
                    for _count in 0..right_factor {
                        backup_options.push([cur_coord[0] - 1, cur_coord[1]]);
                    }
                }
            }
            if cur_coord[0] != self.width - 1 {
                if !coordinate_set.contains(&[cur_coord[0] + 1, cur_coord[1]]) {
                    for _count in 0..right_factor {
                        options.push([cur_coord[0] + 1, cur_coord[1]]);
                    }
                } else {
                    for _count in 0..left_factor {
                        backup_options.push([cur_coord[0] + 1, cur_coord[1]]);
                    }
                }
            }  
            if cur_coord[1] != 0 {
                if !coordinate_set.contains(&[cur_coord[0], cur_coord[1] - 1]) {
                    for _count in 0..up_factor {
                        options.push([cur_coord[0], cur_coord[1] - 1]);
                    }
                } else {
                    for _count in 0..down_factor {
                        backup_options.push([cur_coord[0], cur_coord[1] - 1]);
                    }
                }
            }
            if cur_coord[1] != self.height - 1 {
                if !coordinate_set.contains(&[cur_coord[0], cur_coord[1] + 1]) {
                    for _count in 0..down_factor {
                        options.push([cur_coord[0], cur_coord[1] + 1]);
                    }
                } else {
                    for _count in 0..up_factor {
                        backup_options.push([cur_coord[0], cur_coord[1] + 1]);
                    }
                }
            }  
            if options.len() != 0 {
                let pick = options[rng.gen_range(0..options.len())];
                cur_coord = pick;
                coordinate_set.insert(pick);
            } else {
                cur_coord = backup_options[rng.gen_range(0..backup_options.len())];
            }
        }
        completed_coords
    }

    fn erode_tiles(&mut self, exact_flag: bool, repeat_until_unmodified: bool, direction_type: bool, state: bool, compare_state: bool, weight: f32, flag: usize, rng: &mut ThreadRng) {
        let tiles = self.select_random_tiles(state, weight, rng);
        let mut modified = false;
        for tile in tiles {
            let surrounding_tiles_flag = self.check_surrounding_tiles(direction_type, compare_state, tile[0], tile[1]);
            let check = if !exact_flag {
                surrounding_tiles_flag >= flag
            } else {
                surrounding_tiles_flag == flag
            };
            if check {
                self.set_tile(tile[0], tile[1], !state);
                if !modified {
                    modified = true;
                }
            }
        }
        if repeat_until_unmodified && modified {
            self.erode_tiles(false, true, direction_type, state, compare_state, weight, flag, rng);
        }
    }

    pub fn instantiate_basic_erosion_maze(width: usize, height: usize, weight: f32, rng: &mut ThreadRng) -> MazeManager {
        let mut map = MazeManager::instantiate_map(width, height, false);
        map.erode_tiles(false, false, true, false, false, weight, 7, rng);
        map
    }

    pub fn instantiate_random_tunnel_maze(width: usize, height: usize, weight: f32, directions_weighted_by_coordinates: bool, max_direction_factor: usize, rng: &mut ThreadRng) -> MazeManager {
        let mut map = MazeManager::instantiate_map(width, height, true);
        let coord = map.select_random_tile(true, rng);
        let mut directions = Vec::new();
        if directions_weighted_by_coordinates {
            directions = vec![1; 4];
            let regions_x = map.center_x / max_direction_factor;
            let regions_y = map.center_y / max_direction_factor;
            if coord[0] > map.center_x {
                let dx = coord[0] - map.center_x;
                directions[2] = 1 + dx / regions_x;
            } else {
                let dx = map.center_x - coord[0];
                directions[3] = 1 + dx / regions_x;
            }
            if coord[1] > map.center_y {
                let dy = coord[1] - map.center_y;
                directions[0] = 1 + dy / regions_y;
            } else {
                let dy = map.center_y - coord[1];
                directions[1] = 1 + dy / regions_y;
            }
        } else {
            for _count in 0..4 {
                directions.push(rng.gen_range(1..=max_direction_factor));
            }
        }
        map.create_random_tunnel(false, coord, weight, directions[0], directions[1], directions[2], directions[3], rng);
        map
    }

    pub fn instantiate_tunnel_maze(width: usize, height: usize, fill_holes: bool, weight: f32, original_parent_tunnel_distance: f32, min_min_tunnel_distance: f32, max_min_tunnel_distance: f32, rng: &mut ThreadRng) -> MazeManager {
        let mut map = MazeManager::instantiate_map(width, height, true);
        let tile_pair = map.select_random_tile_pair_from_distance(true, original_parent_tunnel_distance, rng);
        let amount = (map.area * weight) as usize;
        let mut parent_tunnel = map.create_tunnel(true, false, tile_pair[0], tile_pair[1], rng);
        while map.false_tiles.len() < amount {
            let child_coord = parent_tunnel[0];
            let child_final_coord = map.select_random_tile_from_distance(true, child_coord, rng.gen_range(min_min_tunnel_distance..=max_min_tunnel_distance), rng);
            parent_tunnel = map.create_tunnel(false, false, child_coord, child_final_coord, rng);
        }
        if fill_holes {
            map.erode_tiles(false, true, false, true, false, 1.0, 4, rng);
        }
        map
    }

    pub fn instantiate_branch_maze(width: usize, height: usize, fill_holes: bool, min_branch_length: usize, max_branch_length: usize, weight: f32, branch_direction_sway_factor: f32, rng: &mut ThreadRng) -> MazeManager {
        let mut map = MazeManager::instantiate_map(width, height, true);
        let tile = map.select_random_tile(true, rng);
        let mut branch_count = 0;
        let mut direction_x_counters = [0, 0];
        let mut direction_y_counters = [0, 0];
        let direction = if rng.gen_range(0.0..1.0) < 0.5 {
            direction_x_counters[0] += 1;
            false
        } else {
            direction_x_counters[1] += 1;
            true
        };
        let y_weight = if tile[1] > map.center_y {
            [0.33 + ((tile[1] - map.center_y) as f32 / map.center_y as f32) * branch_direction_sway_factor, 0.33]
        } else {
            [0.33, 0.33 + ((map.center_y - tile[1]) as f32 / map.center_y as f32) * branch_direction_sway_factor]
        };
        let mut parent_branch = map.create_horizontal_branch(true, false, direction, tile[0], tile[1], rng.gen_range(min_branch_length..=max_branch_length), y_weight[0], y_weight[1], rng);
        let amount = (map.area * weight) as usize;
        while map.false_tiles.len() < amount {
            let child_coord = select_random_tile_from_set(&parent_branch, rng);
            if branch_count % 2 == 0 {
                let direction = if direction_y_counters[0] < direction_y_counters[1] {
                    direction_y_counters[0] += 1;
                    false
                } else {
                    direction_y_counters[1] += 1;
                    true
                };
                let y_weight = if child_coord[1] > map.center_y {
                    [0.33 + ((child_coord[1] - map.center_y) as f32 / map.center_y as f32) * branch_direction_sway_factor, 0.33]
                } else {
                    [0.33, 0.33 + ((map.center_y - child_coord[1]) as f32 / map.center_y as f32) * branch_direction_sway_factor]
                };
                parent_branch = map.create_vertical_branch(false, false, direction, child_coord[0], child_coord[1], rng.gen_range(min_branch_length..=max_branch_length), y_weight[0], y_weight[1], rng);
            } else {
                let direction = if direction_x_counters[0] < direction_x_counters[1] {
                    direction_x_counters[0] += 1;
                    false
                } else {
                    direction_x_counters[1] += 1;
                    true
                };
                let x_weight = if child_coord[0] > map.center_x {
                    [0.33 + ((child_coord[0] - map.center_x) as f32 / map.center_y as f32) * branch_direction_sway_factor, 0.33]
                } else {
                    [0.33, 0.33 + ((map.center_x - child_coord[0]) as f32 / map.center_y as f32) * branch_direction_sway_factor]
                };
                parent_branch = map.create_horizontal_branch(false, false, direction, child_coord[0], child_coord[1], rng.gen_range(min_branch_length..=max_branch_length), x_weight[0], x_weight[1], rng);
            }
            branch_count += 1;
        }
        if fill_holes {
            map.erode_tiles(false, true, false, true, false, 1.0, 4, rng);
        }
        map
    }

    pub fn instantiate_branched_erosion_maze(width: usize, height: usize, weight: f32, horizontal_step: usize, vertical_step: usize, rng: &mut ThreadRng) -> MazeManager {
        let mut map = MazeManager::instantiate_map(width, height, true);
        for y in (0..height).step_by(horizontal_step) {
            map.create_horizontal_branch(true, false, true, 0, y, width, 0.33, 0.33, rng);
        }
        for x in (0..width).step_by(vertical_step) {
            map.create_vertical_branch(true, false, true, x, 0, height, 0.33, 0.33, rng);
        }
        map.erode_tiles(false, false, true, true, false, weight, 7, rng);
        map
    }

    fn square_quantize_check_next_square_magnitude(&self, coord: [usize; 2], depth: usize, occupied: &Vec<Vec<bool>>) -> bool {
        if coord[0] + depth >= self.width || coord[1] + depth >= self.height {
            return false;
        }
        for y in coord[1]..=coord[1] + depth {
            for x in coord[0]..=coord[0] + depth {
                if occupied[y][x] || !self.grid[y][x] {
                    return false;
                }
            }
        }
        true
    }
    
    pub fn square_quantize_maze(&self) -> Vec<[usize; 3]> {
        let mut y = 0;
        let mut x = 0;
        let mut completed_squares = Vec::new();
        let mut occupied = vec![vec![false; self.width]; self.height];
        while y < self.height {
            if occupied[y][x] {
                x += 1;
                if x == self.width {
                    x = 0;
                    y += 1;
                }
                continue;
            }
            if self.grid[y][x] {
                let mut square_magnitude = 1;
                while self.square_quantize_check_next_square_magnitude([x, y], square_magnitude, &occupied) {
                    square_magnitude += 1;
                }
                completed_squares.push([x, y, square_magnitude]);
                for yy in y..y + square_magnitude {
                    for xx in x..x + square_magnitude {
                        occupied[yy][xx] = true;
                    }
                }
                x += square_magnitude;
                if x == self.width {
                    x = 0;
                    y += 1;
                }
            } else {
                x += 1;
                if x == self.width {
                    x = 0;
                    y += 1;
                }
            }
        }
        completed_squares
    }
}
