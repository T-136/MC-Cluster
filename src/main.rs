use clap::Parser;
use mc::Simulation;
use std::fs;
use std::panic;
use std::thread;

fn fmt_scient(num: &str) -> u64 {
    let mut parts = num.split(['e', 'E']);

    let pre_num = parts.next().unwrap();
    let exp = parts.next().unwrap_or(&"1");
    if parts.next().is_some() {
        panic!("wrong iterations input");
    }

    let base: u64 = 10;
    pre_num.parse::<u64>().expect("wrong iterations input")
        * base.pow(exp.parse::<u32>().expect("wrong iterations input"))
}

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
    #[arg(short, long)]
    iterations: String,

    #[arg(short, long, default_value_t = 300.)]
    temperature: f64,

    #[arg(short, long)]
    begin_temperature: Option<f64>,

    #[arg(short, long, value_delimiter = '-', default_values_t = vec!(0,1))]
    repetition: Vec<usize>,

    #[arg(short, long, default_value_t = String::from("../171717-pair"))]
    grid_folder: String,

    #[arg(short, long, default_value_t = String::from("../input_cluster/bulk.poscar"))]
    core_file: String,

    #[arg(long)]
    last_frames_xyz: Option<u64>,

    #[arg(short, long, default_value_t = 1)]
    write_trajectory_frequency: u64,

    #[arg(short, long, default_value_t = 0.50)]
    optimization_cut_off_perc: f64,

    #[arg(short, long, allow_hyphen_values = true)]
    unique_levels: i32,
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
    let atoms_input = args.atoms_input;
    let unique_levels = args.unique_levels;
    if !std::path::Path::new(&save_folder).exists() {
        fs::create_dir_all(&save_folder).unwrap();
    }

    let input_file: Option<String> = args.start_cluster;

    #[allow(unused_variables)]
    let (pairlist_file, nn_pairlist_file, nnn_pairlist_file, atom_sites) =
        file_paths(args.grid_folder);

    let niter_str = args.iterations;
    let niter = fmt_scient(&niter_str);
    let last_traj_frequency: u64 = args.write_trajectory_frequency;
    let last_frames_trajectory: Option<u64> = args.last_frames_xyz;
    let bulk_file_name: String = args.core_file;
    let optimization_cut_off_perc: f64 = args.optimization_cut_off_perc;

    let repetition = args.repetition;

    let mut handle_vec = Vec::new();
    for rep in repetition[0]..repetition[1] {
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
                input_file,
                atoms_input,
                temperature,
                start_temperature,
                save_folder,
                pairlist_file,
                nn_pairlist_file,
                // nnn_pairlist_file,
                atom_sites,
                last_traj_frequency,
                last_frames_trajectory,
                bulk_file_name,
                rep,
                optimization_cut_off_perc,
            );
            let exp = sim.run(unique_levels);
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
