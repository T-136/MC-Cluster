use clap::Parser;
use mc::Simulation;
use std::fs;
use std::thread;
// use time_graph::enable_data_collection;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short, long)]
    start_cluster: Option<String>,

    #[arg(short, long)]
    atoms_input: Option<u32>,

    /// folder to store results
    #[arg(short, long, default_value_t = String::from("./sim/"))]
    folder: String,

    /// iterations
    #[arg(short, long, default_value_t = 50000)]
    iterations: u64,

    #[arg(short, long, default_value_t = 300.)]
    temperature: f64,

    #[arg(short, long)]
    begin_temperature: Option<f64>,

    #[arg(short, long, default_value_t = 1)]
    repetition: usize,

    #[arg(short, long, default_value_t = String::from("../555-pair"))]
    grid_folder: String,

    #[arg(short, long, default_value_t = String::from("../input_cluster/bulk.poscar"))]
    core_file: String,

    #[arg(short, long)]
    last_frames_trajectory: Option<u64>,

    #[arg(short, long)]
    write_trajectory_frequency: Option<u64>,

    #[arg(short, long, default_value_t = 0.50)]
    optimization_cut_off_perc: f64,
}

fn file_paths(grid_folder: String) -> (String, String, String, String) {
    (
        format!("{}/pairlist", grid_folder),
        format!("{}/nn_pairlist", grid_folder),
        format!("{}/nnn_pairlist", grid_folder),
        format!("{}/atom_sites", grid_folder),
    )
}

fn main() {
    // enable_data_collection(true);
    println!("determined next-nearest neighbor list");

    let args = Args::parse();
    let save_folder: String = args.folder;
    let temperature: f64 = args.temperature;
    let start_temperature: Option<f64> = args.begin_temperature;
    let nsites: u32 = 15 * 15 * 15 * 4;
    let atoms_input = args.atoms_input;
    if !std::path::Path::new(&save_folder).exists() {
        fs::create_dir_all(&save_folder).unwrap();
    }

    let input_file: Option<String> = args.start_cluster;

    #[allow(unused_variables)]
    let (pairlist_file, nn_pairlist_file, nnn_pairlist_file, atom_sites) =
        file_paths(args.grid_folder);

    let niter = args.iterations;
    let trajectory_frequency: Option<u64> = args.write_trajectory_frequency;
    let last_frames_trajectory: Option<u64> = args.last_frames_trajectory;
    let bulk_file_name: String = args.core_file;
    let optimization_cut_off_perc: f64 = args.optimization_cut_off_perc;

    let repetition = args.repetition;

    let mut handle_vec = Vec::new();
    for rep in 0..repetition {
        let input_file = input_file.clone();
        let save_folder = save_folder.clone();
        let pairlist_file = pairlist_file.clone();
        let nn_pairlist_file = nn_pairlist_file.clone();
        // let nnn_pairlist_file = nnn_pairlist_file.clone();
        let atom_sites = atom_sites.clone();
        let bulk_file_name = bulk_file_name.clone();
        // let optimization_cut_off_perc = optimization_cut_off_perc.clone();
        handle_vec.push(thread::spawn(move || {
            let mut sim = Simulation::new(
                niter,
                nsites,
                input_file,
                atoms_input,
                temperature,
                start_temperature,
                save_folder,
                pairlist_file,
                nn_pairlist_file,
                // nnn_pairlist_file,
                atom_sites,
                trajectory_frequency,
                last_frames_trajectory,
                bulk_file_name,
                rep,
                optimization_cut_off_perc,
            );
            let exp = sim.run();
            sim.write_exp_file(&exp);
        }));
    }
    for handle in handle_vec {
        handle.join().unwrap();
    }
    // mc::find_simulation_with_lowest_energy(save_folder).unwrap_or_else(|err| {
    //     println!(
    //         "{:?}",
    //         format!("deleting folders with heigh energy not successful {err}")
    //     )
    // });
}
