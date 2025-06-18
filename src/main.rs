use chemfiles::{Frame, Trajectory};
use clap::{ArgGroup, Parser};
use core::panic;
use std::collections::HashMap;
use std::io::BufReader;
use std::sync::Arc;
use std::{fs, thread};
use MC_Cluster::energy::{EnergyInput, EnergyValues};
use MC_Cluster::{CreateStructure, GridStructure, Simulation, Structure};

fn atoms_input(atom_name: &str, atom_names: &mut MC_Cluster::AtomNames) {
    if let Some(supp) = atom_names.support.as_ref() {
        if atom_name == supp {
            return;
        }
    }
    if atom_names.atom.is_none() {
        atom_names.atom = Some(atom_name.to_string());
    } else if let Some(atom) = atom_names.atom.as_ref() {
        if atom != atom_name {
            panic!("to many atoms in input file");
        }
    }
}

fn read_sample(
    input_file: &str,
    atom_names: &mut MC_Cluster::AtomNames,
) -> Vec<(String, [f64; 3])> {
    if input_file.contains(".poscar") {
        todo!();
        // let newatoms = Poscar::from_path(input_file).unwrap();
        // println!("poscar naem {:?}", newatoms.group_symbols().unwrap().next());
        // let xyz = newatoms.scaled_cart_positions();
        // xyz.into_owned()
    } else if input_file.contains(".xyz") {
        let mut trajectory = Trajectory::open(input_file, 'r').unwrap();
        let mut frame = Frame::new();
        trajectory.read(&mut frame).unwrap();
        let mut atom_vec: Vec<(String, [f64; 3])> = Vec::new();
        let positions = frame.positions().to_owned();
        for (i, atom) in frame.iter_atoms().enumerate() {
            atoms_input(&atom.name(), atom_names);
            atom_vec.push((atom.name().clone(), positions[i]));
        }
        atom_vec
    } else {
        panic!("no .poscar or .xyz, cant read file");
    }
}

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

fn collect_energy_values<const N: usize>(inp: String) -> EnergyValues<[i64; N]> {
    let json =
        if inp.chars().next().unwrap().is_numeric() || inp.starts_with('-') || inp.starts_with('{')
        {
            let res: Result<HashMap<String, Vec<i64>, fnv::FnvBuildHasher>, serde_json::Error> =
                serde_json::from_str(&inp);
            res
        } else if inp.ends_with(".json") {
            let file = fs::File::open(inp).expect("can't find energy file");
            let reader = BufReader::new(file);
            let res: Result<HashMap<String, Vec<i64>, fnv::FnvBuildHasher>, serde_json::Error> =
                serde_json::from_reader(reader);
            res
        } else {
            panic!("energy input is neither JSON nor a file path");
            // fs::read_to_string(inp).expect("can't find energy file")
        };
    let mut energy: [i64; N] = [0; N];
    #[allow(non_snake_case)]
    let mut CO_ads: [i64; N] = [0; N];
    for (i, val) in json
        .as_ref()
        .unwrap()
        .get("CN_energy")
        .unwrap()
        .iter()
        .enumerate()
    {
        energy[i] = *val;
    }
    #[allow(non_snake_case)]
    let CO_ads_opt = if let Some(it) = json.unwrap().get("ads_e_CO") {
        for (i, val) in it.iter().enumerate() {
            CO_ads[i] = *val;
        }
        Some(CO_ads)
    } else {
        None
    };
    EnergyValues {
        complet_energy: energy,
        co_ads_energy: CO_ads_opt,
    }
}

#[derive(Parser, Debug, Clone)]
#[clap(group(
        ArgGroup::new("startstructure")
            .required(true)
            .multiple(true)
            .args(&["start_cluster", "atoms", "support"]),
    ))]
struct StartStructure {
    #[arg(short, long, group = "mode", conflicts_with_all = ["atoms"])]
    start_cluster: Option<String>,

    /// Atom name and number of that atom seperated by a comma. "Pt,4000"
    #[arg(short, long, value_delimiter = ',', allow_hyphen_values(true))]
    atoms: Option<Vec<String>>,

    /// When creating a new particle using the atoms flag, write the Atom name and a vector orthogonal to the support
    /// surface Al,1,1,1. When starting from a xyz file only the support atom is required.
    #[arg(short, long, value_delimiter = ',', allow_hyphen_values(true))]
    support: Option<Vec<String>>,
}

// #[derive(Parser, Debug, Clone)]
// #[clap(group(
//         ArgGroup::new("createcluster")
//             .multiple(true)
//             .args(&["atoms", "support"]),
//     ))]
// struct Createcluster {
//     /// Atom name and number of that atom seperated by a comma. "Pt,4000"
//     #[arg(short, long, value_delimiter = ',', allow_hyphen_values(true))]
//     atoms: Vec<String>,
//
//     #[arg(short, long, value_delimiter = ',', allow_hyphen_values(true))]
//     support: Option<Vec<String>>,
// }

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

    /// Output folder
    #[arg(short, long, default_value_t = String::from("./sim/"))]
    folder: String,

    /// Determines the length of the simulation
    #[arg(short, long)]
    iterations: String,

    /// Temperature where the annealing process starts
    #[arg(short, long, default_value_t = 5000.)]
    begin_temperature: f64,

    /// Temperature at which the annealing process stops
    #[arg(short, long, default_value_t = 300.)]
    temperature: f64,

    /// Fraction of the simulation after which the annealing process is completed. After that, the
    /// temperature remains constant.
    #[arg(short, long, value_delimiter = '/', default_values_t = vec!(1,2))]
    optimization_cut_off_fraction: Vec<u64>,

    /// Support energy
    #[arg(long, allow_hyphen_values(true))]
    support_e: Option<i64>,

    /// File path or string containing a JSON-formatted list of energies
    #[arg(long, allow_hyphen_values(true))]
    e_l_cn: Option<String>,

    /// File path or string containing JSON-formatted energies
    #[arg(long, allow_hyphen_values(true))]
    e_cn: Option<String>,

    /// How many times the same simulation is run. Multiple runs allow for convergence tests.
    /// The number will be part of the simulation folder name. After running `-r 0-1`, you can run `-r 1-2`
    /// and the previous simulation will not be overwritten.
    #[arg(short, long, value_delimiter = '-', default_values_t = vec!(0,1))]
    repetition: Vec<usize>,

    /// Folder containing the setup files like neighbor sites. It can be created using the Python
    /// script `create_sites.py`.
    #[arg(short, long, default_value_t = String::from("../303030-pair"))]
    grid_folder: String,

    /// Track the simulation by taking snapshots of the cluster throughout the simulation
    #[arg(short, long, default_value_t = false)]
    xyz_trajectory: bool,

    /// Generate a heat map
    #[arg(long, default_value_t = false)]
    heat_map: bool,
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

fn unpack_support_input(atoms_opt: Option<Vec<String>>) -> (Option<String>, Option<Vec<i32>>) {
    if let Some(atoms) = atoms_opt {
        if atoms.len() == 1 {
            (Some(atoms[0].clone()), None)
        } else if atoms.len() == 4 {
            (
                Some(atoms[0].clone()),
                Some(vec![
                    atoms[1].parse::<i32>().expect("wrong support vector"),
                    atoms[2].parse::<i32>().expect("wrong support vector"),
                    atoms[3].parse::<i32>().expect("wrong support vector"),
                ]),
            )
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

    let mut atom_names = MC_Cluster::AtomNames::default();
    let args = Args::parse();
    let save_folder: String = args.folder;
    let temperature: f64 = args.temperature;
    let start_temperature: f64 = args.begin_temperature;
    if !std::path::Path::new(&save_folder).exists() {
        fs::create_dir_all(&save_folder).unwrap();
    }

    if let Some(supp_atom) = args.start_structure.support.clone() {
        atom_names.support = Some(supp_atom[0].clone());
    }
    let start_structure = if let Some(start_cluster) = args.start_structure.start_cluster {
        let xyz = read_sample(&start_cluster, &mut atom_names);
        Structure::StartStructure(Arc::new(xyz))
    } else if let Some(atoms) = args.start_structure.atoms {
        let (supp_atom_name, support_indices) = unpack_support_input(args.start_structure.support);
        let (atom_name, atom_count) = (
            atoms[0].clone(),
            atoms[1].parse().expect("bad number of atoms"),
        );
        atom_names.atom = Some(atom_name.clone());
        atom_names.support = supp_atom_name.clone();

        Structure::CreateCluster(CreateStructure {
            atom_name,
            atom_count,
            support_vector: support_indices,
            support_atom_name: supp_atom_name,
        })
    } else {
        panic!("either file path to start structure or vlaues for creating a cluster are required");
    };
    println!("atom_names {:?}", atom_names);

    let support_e = args.support_e.unwrap_or(0);

    let grid_folder: String = Args::parse().grid_folder;

    let niter_str = args.iterations;
    let niter = fmt_scient(&niter_str);
    let mut write_snap_shots: bool = args.xyz_trajectory;
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
        EnergyInput::LinearCn(collect_energy_values(args.e_l_cn.unwrap()))
    } else if args.e_cn.is_some() {
        EnergyInput::Cn(collect_energy_values(args.e_cn.unwrap()))
    } else {
        panic!("no energy")
    };

    println!("energy: {:?}", energy);
    println!("{:?}", repetition);

    let mut handle_vec = Vec::new();
    let gridstructure: GridStructure = GridStructure::new(file_paths(grid_folder));

    let gridstructure = Arc::new(gridstructure);

    for rep in repetition[0]..repetition[1] {
        let save_folder = save_folder.clone();
        let optimization_cut_off_fraction = optimization_cut_off_fraction.clone();
        let energy = energy.clone();
        let gridstructure_arc = Arc::clone(&gridstructure);
        let start_structure = start_structure.clone();
        let atom_names = atom_names.clone();

        handle_vec.push(thread::spawn(move || {
            let mut sim = Simulation::new(
                atom_names,
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
        }));
    }
    for handle in handle_vec {
        handle.join().unwrap();
    }
    // MC_Cluster::find_simulation_with_lowest_energy(save_folder).unwrap_or_else(|err| {
    //     println!(
    //         "{:?}",
    //         format!("deleting folders with heigh energy not successful {err}")
    //     )
    // });
}
