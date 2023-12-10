use clap::ArgGroup;
use clap::Parser;
use mc::EnergyInput;
use mc::Simulation;
use std::fs;
use std::panic;
use std::thread;

fn fmt_scient(num: &str) -> u64 {
    let mut parts = num.split(['e', 'E']);

    let pre_num = parts.next().unwrap();
    let exp = parts.next().unwrap_or("0");
    if parts.next().is_some() {
        panic!("wrong iterations input");
    }

    let base: u64 = 10;
    pre_num.parse::<u64>().expect("wrong iterations input")
        * base.pow(exp.parse::<u32>().expect("wrong iterations input"))
}

fn file_or_energy(inp: &str, energy: mc::EnergyInput) -> mc::EnergyInput {
    if inp.starts_with("./") {
        let contents = fs::read_to_string(inp).expect("can't find energy file");

        let mut string_iter = contents.split([';', ',']);
        match energy {
            EnergyInput::LinearCn(energy_vec) => {
                energy_vec.map(|_| string_iter.next().unwrap().parse::<i64>().unwrap());
                EnergyInput::LinearCn(energy_vec)
            }
            EnergyInput::Cn(energy_vec) => {
                energy_vec.map(|_| string_iter.next().unwrap().parse::<i64>().unwrap());
                EnergyInput::Cn(energy_vec)
            }
            EnergyInput::LinearGcn(energy_vec) => {
                energy_vec.map(|_| string_iter.next().unwrap().parse::<i64>().unwrap());
                EnergyInput::LinearGcn(energy_vec)
            }
            EnergyInput::Gcn(energy_vec) => {
                energy_vec.map(|_| string_iter.next().unwrap().parse::<i64>().unwrap());
                EnergyInput::Gcn(energy_vec)
            }
        }
    } else {
        let mut string_iter = inp.split(',');
        // println!("string iter {:?}", string_iter.collect());
        println!("string iter {:?}", string_iter);

        match energy {
            EnergyInput::LinearCn(mut energy_vec) => {
                for x in energy_vec.iter_mut() {
                    *x = string_iter.next().unwrap().parse::<i64>().unwrap()
                }
                EnergyInput::LinearCn(energy_vec)
            }
            EnergyInput::Cn(energy_vec) => {
                energy_vec.map(|_| string_iter.next().unwrap().parse::<i64>().unwrap());
                EnergyInput::Cn(energy_vec)
            }
            EnergyInput::LinearGcn(energy_vec) => {
                energy_vec.map(|_| string_iter.next().unwrap().parse::<i64>().unwrap());
                EnergyInput::LinearGcn(energy_vec)
            }
            EnergyInput::Gcn(energy_vec) => {
                energy_vec.map(|_| string_iter.next().unwrap().parse::<i64>().unwrap());
                EnergyInput::Gcn(energy_vec)
            }
        }
    }
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[clap(group(
        ArgGroup::new("my-group")
            .required(true)
            .args(&["e_l_cn", "e_cn", "e_l_gcn", "e_gcn"]),
    ))]
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

    #[arg(long, allow_hyphen_values(true))]
    e_l_cn: Option<String>,

    #[arg(long, allow_hyphen_values(true))]
    e_cn: Option<String>,

    #[arg(long, allow_hyphen_values(true))]
    e_l_gcn: Option<String>,

    #[arg(long, allow_hyphen_values(true))]
    e_gcn: Option<String>,

    #[arg(short, long, value_delimiter = '-', default_values_t = vec!(0,1))]
    repetition: Vec<usize>,

    #[arg(short, long, default_value_t = String::from("../999-pair"))]
    grid_folder: String,

    #[arg(short, long, default_value_t = String::from("../input_cluster/bulk.poscar"))]
    core_file: String,

    #[arg(short, long, default_value_t = false)]
    write_snap_shots: bool,

    #[arg(long, default_value_t = false)]
    heat_map: bool,

    #[arg(short, long, value_delimiter = '/', default_values_t = vec!(1,2))]
    optimization_cut_off_fraction: Vec<u64>,

    #[arg(short, long, allow_hyphen_values = true)]
    unique_levels: i32,
}

fn file_paths(grid_folder: String) -> (String, String, String, String, String, String) {
    (
        format!("{}/pairlist", grid_folder),
        format!("{}/nn_pairlist", grid_folder),
        format!("{}/nnn_pairlist", grid_folder),
        format!("{}/atom_sites", grid_folder),
        format!("{}/nn_pair_no_intersec", grid_folder),
        format!("{}/nnn_pair_no_intersec", grid_folder),
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
    let (
        pairlist_file,
        nn_pairlist_file,
        nnn_pairlist_file,
        atom_sites,
        nn_pair_no_int_file,
        nnn_pair_no_int_file,
    ) = file_paths(args.grid_folder);

    let niter_str = args.iterations;
    let niter = fmt_scient(&niter_str);
    let mut write_snap_shots: bool = args.write_snap_shots;
    let heat_map: bool = args.heat_map;
    if heat_map {
        write_snap_shots = true;
    }
    let bulk_file_name: String = args.core_file;
    let optimization_cut_off_fraction: Vec<u64> = args.optimization_cut_off_fraction;

    let repetition = args.repetition;

    let energy = if args.e_l_cn.is_some() {
        file_or_energy(&args.e_l_cn.unwrap(), EnergyInput::LinearCn([0; 2]))
    } else if args.e_l_cn.is_some() {
        file_or_energy(&args.e_cn.unwrap(), EnergyInput::Cn([0; 13]))
    } else if args.e_l_cn.is_some() {
        file_or_energy(&args.e_l_gcn.unwrap(), EnergyInput::LinearGcn([0; 2]))
    } else if args.e_l_cn.is_some() {
        file_or_energy(&args.e_gcn.unwrap(), EnergyInput::Gcn([0; 60]))
    } else {
        panic!("no energy")
    };

    println!("energy: {:?}", energy);
    println!("{:?}", repetition);

    let mut handle_vec = Vec::new();
    for rep in repetition[0]..repetition[1] {
        let input_file = input_file.clone();
        let save_folder = save_folder.clone();
        let pairlist_file = pairlist_file.clone();
        let nn_pair_no_int_file = nn_pair_no_int_file.clone();
        let nnn_pair_no_int_file = nnn_pair_no_int_file.clone();
        // let nn_pairlist_file = nn_pairlist_file.clone();
        // let nnn_pairlist_file = nnn_pairlist_file.clone();
        let atom_sites = atom_sites.clone();
        let bulk_file_name = bulk_file_name.clone();
        let optimization_cut_off_fraction = optimization_cut_off_fraction.clone();
        let energy = energy.clone();
        handle_vec.push(thread::spawn(move || {
            let mut sim = Simulation::new(
                niter,
                input_file,
                atoms_input,
                temperature,
                start_temperature,
                save_folder,
                pairlist_file,
                nn_pair_no_int_file,
                nnn_pair_no_int_file,
                // nn_pairlist_file,
                // nnn_pairlist_file,
                atom_sites,
                write_snap_shots,
                heat_map,
                bulk_file_name,
                rep,
                optimization_cut_off_fraction,
                energy,
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
