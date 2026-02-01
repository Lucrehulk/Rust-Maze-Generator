mod mazegen;
use crate::mazegen::MazeManager;

fn main() {
    /* Some example generations */
    // let world = MazeManager::instantiate_basic_erosion_maze(100, 100, 0.25, &mut rand::rng());
    // let world = MazeManager::instantiate_random_tunnel_maze(100, 100, 0.75, true, 3, &mut rand::rng());
    // let world = MazeManager::instantiate_tunnel_maze(100, 100, true, 0.2, 140.0, 50.0, 100.0, &mut rand::rng());
    // let world = MazeManager::instantiate_branch_maze(100, 100, true, 10, 100, 0.4, 0.3, &mut rand::rng());
    let world = MazeManager::instantiate_branched_erosion_maze(100, 100, 0.5, 3, 3, &mut rand::rng());
    println!("{:?}", world.grid);
    println!("{:?}", world.square_quantize_maze());
}
