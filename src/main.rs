use clap::ArgGroup;
use clap::Parser;
use core::panic;
use lazy_static::lazy_static;
use mc::energy;
use mc::energy::EnergyInput;
use mc::energy::EnergyValues;
use mc::GridStructure;
use mc::Simulation;
use std::collections::HashMap;
use std::fs;
use std::io::BufReader;
use std::sync::Arc;
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

fn collect_energy_values<const N: usize>(
    mut energy_vec: [i64; N],
    inp: String,
) -> EnergyValues<[i64; N]> {
    let contents = if inp.chars().next().unwrap().is_numeric() || inp.starts_with('-') {
        inp
    } else if inp.ends_with(".json") {
        let file = fs::File::open(inp).expect("can't find energy file");
        let reader = BufReader::new(file);
        // let res: Result<Results, serde_json::Error> = serde_json::from_reader(reader);
        let res: Result<HashMap<String, Vec<i64>, fnv::FnvBuildHasher>, serde_json::Error> =
            serde_json::from_reader(reader);
        let mut energy: [i64; N] = [0; N];
        let mut CO_ads: [i64; N] = [0; N];
        for (i, val) in res
            .as_ref()
            .unwrap()
            .get("energy")
            .unwrap()
            .iter()
            .enumerate()
        {
            energy[i] = *val;
        }
        let CO_ads_opt = if let Some(it) = res.unwrap().get("CO_ads") {
            for (i, val) in it.iter().enumerate() {
                CO_ads[i] = *val;
            }
            Some(CO_ads)
        } else {
            None
        };
        return EnergyValues {
            complet_energy: energy,
            co_ads_energy: CO_ads_opt,
        };
    } else {
        fs::read_to_string(inp).expect("can't find energy file")
    };
    let mut string_iter = contents.trim().split(',');
    for x in energy_vec.iter_mut() {
        *x = string_iter
            .next()
            .unwrap()
            .trim()
            .parse::<i64>()
            .unwrap_or_else(|err| {
                panic!(
                    "iter received from input file: {:?}, err: {}",
                    string_iter, err
                )
            });
    }
    return EnergyValues {
        complet_energy: energy_vec,
        co_ads_energy: None,
    };
}

#[derive(Parser, Debug)]
#[clap(group(
        ArgGroup::new("startstructure")
            .required(true)
            .args(&["start_cluster", "atoms_input"]),
    ))]
struct StartStructure {
    #[arg(short, long)]
    start_cluster: Option<String>,

    #[arg(short, long, value_delimiter = ',')]
    atoms_input: Option<Vec<u32>>,
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[clap(group(
        ArgGroup::new("energy")
            .multiple(true)
            .required(true)
            .args(&["e_l_cn", "e_cn"]),
    ))]
struct Args {
    #[clap(flatten)]
    start_structure: StartStructure,

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
    support_e: Option<i64>,

    #[arg(long, allow_hyphen_values(true))]
    e_l_cn: Option<String>,

    #[arg(long, allow_hyphen_values(true))]
    e_cn: Option<String>,

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

fn file_paths(grid_folder: String) -> (String, String, String, String, String) {
    (
        format!("{}/nearest_neighbor", grid_folder),
        format!("{}/next_nearest_neighbor", grid_folder),
        format!("{}/nn_pairlist", grid_folder),
        format!("{}/atom_sites", grid_folder),
        format!("{}/nn_pair_no_intersec", grid_folder),
    )
}

fn unpack_atoms_input(atoms: Vec<u32>) -> (Option<u32>, Option<Vec<u32>>) {
    if atoms.len() == 1 {
        return (Some(atoms[0]), None);
    } else if atoms.len() == 4 {
        return (Some(atoms[0]), Some(vec![atoms[1], atoms[2], atoms[3]]));
    } else {
        panic!("wrong atoms input, user one of the two input options: \n number of atoms: '-a x' \n or number of atoms with miller indices: '-a x h k l' ")
    }
}

fn main() {
    // enable_data_collection(true);
    println!("determined next-nearest neighbor list");

    let args = Args::parse();
    let save_folder: String = args.folder;
    let temperature: f64 = args.temperature;
    let start_temperature: Option<f64> = args.begin_temperature;
    let unique_levels = args.unique_levels;
    if !std::path::Path::new(&save_folder).exists() {
        fs::create_dir_all(&save_folder).unwrap();
    }

    // let (atoms_input, sup) =
    let (atoms_input, support_indices) = if let Some(atoms_input) = args.start_structure.atoms_input
    {
        unpack_atoms_input(atoms_input)
    } else {
        (None, None)
    };

    let input_file: Option<String> = args.start_structure.start_cluster;

    lazy_static! {

    static ref grid_folder: String = Args::parse().grid_folder;

    // static ref pairlist_file: String = format!("{}/nearest_neighbor", grid_folder);
    static ref pairlist_file: String = grid_folder.clone() + "/nearest_neighbor";
    static ref n_pairlist_file: String = grid_folder.clone() + "/next_nearest_neighbor";
    static ref nn_pairlist_file: String = grid_folder.clone() + "/nn_pairlist";
    static ref atom_sites: String = grid_folder.clone() + "/atom_sites";
    static ref nn_pair_no_int_file: String = grid_folder.clone() + "/nn_pair_no_intersec";
    static ref bulk_file_name: String = Args::parse().core_file;
    }

    let niter_str = args.iterations;
    let niter = fmt_scient(&niter_str);
    let mut write_snap_shots: bool = args.write_snap_shots;
    let heat_map: bool = args.heat_map;
    if heat_map {
        write_snap_shots = true;
    }
    let optimization_cut_off_fraction: Vec<u64> = args.optimization_cut_off_fraction;
    let repetition = args.repetition;
    let support_e = args.support_e.unwrap_or(0);

    let energy = if args.e_l_cn.is_some() {
        EnergyInput::LinearCn(collect_energy_values([0; 2], args.e_l_cn.unwrap()))
    } else if args.e_cn.is_some() {
        EnergyInput::Cn(collect_energy_values([0; 13], args.e_cn.unwrap()))
    } else {
        panic!("no energy")
    };

    println!("energy: {:?}", energy);
    println!("{:?}", repetition);

    let mut handle_vec = Vec::new();
    lazy_static! {
        static ref gridstructure: GridStructure = GridStructure::new(
            &pairlist_file,
            // &n_pairlist_file,
            &nn_pair_no_int_file,
            &atom_sites,
            &bulk_file_name,
        );
    }

    // lazy_static! {
    //     static ref gridstructure_stat: GridStructure = gridstructure;
    // }
    // let gridstructure = Arc::new(gridstructure);

    for rep in repetition[0]..repetition[1] {
        let input_file = input_file.clone();
        let save_folder = save_folder.clone();
        let optimization_cut_off_fraction = optimization_cut_off_fraction.clone();
        let energy = energy.clone();
        // let gridstructure_arc = Arc::clone(&gridstructure);
        let support_indices = support_indices.clone();

        handle_vec.push(thread::spawn(move || {
            let mut sim = Simulation::new(
                niter,
                input_file,
                atoms_input,
                temperature,
                start_temperature,
                save_folder,
                write_snap_shots,
                heat_map,
                rep,
                optimization_cut_off_fraction,
                energy,
                support_indices,
                &gridstructure,
                support_e,
            );
            let exp = sim.run(unique_levels);
            sim.write_exp_file(&exp);
            sim.count_empty_sites()
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
