use clap::ArgGroup;
use clap::Parser;
use clap::ValueEnum;
use core::panic;
use mc::energy;
use mc::energy::EnergyInput;
use mc::energy::EnergyValues;
use mc::CreateStructure;
use mc::GridStructure;
use mc::Simulation;
use mc::Structure;
use std::collections::HashMap;
use std::fs;
use std::io::BufReader;
use std::sync::Arc;
use std::thread;

fn prepend<T>(v: Vec<T>, s: &[T]) -> Vec<T>
where
    T: Clone,
{
    let mut tmp: Vec<_> = s.to_owned();
    tmp.extend(v);
    tmp
}

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

#[derive(Parser, Debug, Clone)]
#[clap(group(
        ArgGroup::new("startstructure")
            .required(true)
            .multiple(true)
            .args(&["start_cluster", "atoms", "support"]),
    ))]
struct StartStructure {
    #[arg(short, long, group = "mode", conflicts_with_all = ["atoms", "support"])]
    start_cluster: Option<String>,

    #[command(flatten)]
    create_cluster: Option<Createcluster>,
}

#[derive(Parser, Debug, Clone)]
#[clap(group(
        ArgGroup::new("createcluster")
            .multiple(true)
            .args(&["atoms", "support"]),
    ))]
struct Createcluster {
    /// Atom name and number of that atom seperated by a comma. "Pt,4000"
    #[arg(short, long, value_delimiter = ',', allow_hyphen_values(true))]
    atoms: Vec<String>,

    #[arg(short, long, value_delimiter = ',', allow_hyphen_values(true))]
    support: Option<Vec<String>>,
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
    #[command(flatten)]
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
}

fn file_paths(grid_folder: String) -> (String, String, String, String) {
    (
        format!("{}/nearest_neighbor", grid_folder),
        format!("{}/nn_pair_no_intersec", grid_folder),
        // format!("{}/next_nearest_neighbor", grid_folder),
        // format!("{}/nn_pairlist", grid_folder),
        format!("{}/atom_sites", grid_folder),
        format!("{}/grid_file.xyz", grid_folder),
    )
}

fn unpack_atoms_input(atoms: Vec<String>) -> (String, u32) {
    return (
        atoms[0].clone(),
        atoms[1].parse().expect("bad number of atoms"),
    );
}

fn unpack_support_input(atoms_opt: Option<Vec<String>>) -> (Option<String>, Option<Vec<i32>>) {
    if let Some(atoms) = atoms_opt {
        if atoms.len() == 1 {
            return (Some(atoms[0].clone()), None);
        } else if atoms.len() == 4 {
            return (
                Some(atoms[0].clone()),
                Some(vec![
                    atoms[1].parse::<i32>().expect("wrong support vector"),
                    atoms[2].parse::<i32>().expect("wrong support vector"),
                    atoms[3].parse::<i32>().expect("wrong support vector"),
                ]),
            );
        } else {
            panic!("wrong support input")
            // panic!("wrong support input, use one of the two input options: \n number of atoms: '-a x' \n or number of atoms with miller indices: '-a x h k l' ")
        }
    } else {
        (None, None)
    }
}

fn main() {
    // enable_data_collection(true);
    println!("determined next-nearest neighbor list");

    let args = Args::parse();
    let save_folder: String = args.folder;
    let temperature: f64 = args.temperature;
    let start_temperature: Option<f64> = args.begin_temperature;
    if !std::path::Path::new(&save_folder).exists() {
        fs::create_dir_all(&save_folder).unwrap();
    }

    let start_structure = if let Some(start_cluster) = args.start_structure.start_cluster {
        Structure::StartStructure(start_cluster)
    } else if let Some(create_cluster) = args.start_structure.create_cluster {
        let (supp_atom_name, support_indices) = unpack_support_input(create_cluster.support);
        let (atom_name, atom_count) = unpack_atoms_input(create_cluster.atoms);

        Structure::CreateCluster(CreateStructure {
            atom_name,
            atom_count,
            support_vector: support_indices,
            support_atom_name: supp_atom_name,
        })
    } else {
        panic!("either file path to start structure or vlaues for creating a cluster are required");
    };

    let support_e = args.support_e.unwrap_or(0);

    let grid_folder: String = Args::parse().grid_folder;

    let niter_str = args.iterations;
    let niter = fmt_scient(&niter_str);
    let mut write_snap_shots: bool = args.write_snap_shots;
    let heat_map: bool = args.heat_map;
    if heat_map {
        write_snap_shots = true;
    }
    let optimization_cut_off_fraction: Vec<u64> = args.optimization_cut_off_fraction;
    let repetition = args.repetition;

    let repetition = if repetition.len() == 1 {
        prepend(repetition, &[0])
    } else {
        repetition
    };

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
    let gridstructure: GridStructure = GridStructure::new(file_paths(grid_folder));

    let gridstructure = Arc::new(gridstructure);

    for rep in repetition[0]..repetition[1] {
        // let input_file = input_file.clone();
        let save_folder = save_folder.clone();
        let optimization_cut_off_fraction = optimization_cut_off_fraction.clone();
        let energy = energy.clone();
        let gridstructure_arc = Arc::clone(&gridstructure);
        let start_structure = start_structure.clone();

        handle_vec.push(thread::spawn(move || {
            let mut sim = Simulation::new(
                niter,
                start_structure,
                temperature,
                start_temperature,
                save_folder,
                write_snap_shots,
                heat_map,
                rep,
                optimization_cut_off_fraction,
                energy,
                gridstructure_arc,
                support_e,
            );
            let exp = sim.run();
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
