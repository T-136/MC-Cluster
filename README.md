<!-- --- -->
<!-- geometry: top=2.5cm, bottom=2.5cm, left=2.5cm, right=2.5cm -->
<!-- --- -->

# MC-Cluster a Monte Carlo simulation for nanoparticle

This Monte-Carlo simulation simulates nanoparticles at different temperatures. 
By applying a simulated annealing protocol, the user can find the most stable, entropically realistic structure of a nanoparticle.

![Image](https://github.com/user-attachments/assets/1aaf0711-65a0-447b-8a4c-2bdc36ef0774)
<!-- ![Image pdf](../particle_supp.png "particle with support") -->

## Overview 

### Simple energy model
The energy input is a simple mapping of coordination number to energy in eV. The input can be a direct mapping of every coordination number to a distinct energy. Alternatively, a linear energy increase with coordination number can be used, which speeds up the simulation.

### High performance
By using a simple energy model, the simulation can perform up to 10^11 iterations in one day. With a linear coordination number to energy mapping it can reach up to 10^12 iterations in one day. This allows for the simulation of larger nanoparticles. 

### Multithreaded
The code allows for as many parallel simulations as cores available on the computing node. This enables the user to obtain a statistically relevant ensemble of low-energy structures. Adding parallel simulations has negligible increase in required RAM.

### Track the Simulation using the .xyz-format
Throughout the simulation it is possible to track the energy and save snapshots of the particle.
Furthermore the particle with the lowest energy will be saved as an .xyz-file.

### Support
The simulation can simulate the effect a support can have on the particle.  
A monometallic support can be created by giving a vector that is orthogonal to the support pane.


## Usage

### Requirements

- [Rust](https://www.rust-lang.org/tools/install)
- [Python](https://www.python.org/) used version 3.11.10
- [ASE](https://wiki.fysik.dtu.dk/ase/) *used version: ase-3.22.1*

### Install from git 
```bash
git clone git@github.com:T-136/MC-Cluster.git
```
### Build the binary

Build the program with:
```bash
cargo build -r
```
After the build step is complete the compiled program can be found in your project folder under "./target/release/mc".

### Run the simulation

```bash
./target/release/MC-Cluster -a Pt,1000  --support Al,1,1,1  -t 1000 -i 1e7 -r 0-1  --e-cn ./example_data/cn_input_example.json -o 9/10 -g ./example_data/303030-grid --support-e 0 --xyz-trajectory
```

Use "-h" or "--help" to see the available flags and how to use them. 
```
  -s, --start-cluster <START_CLUSTER>

  -a, --atoms <ATOMS>
          Atom name and number of that atom seperated by a comma. "Pt,4000"
  -s, --support <SUPPORT>
          When creating a new particle using the atoms flag, write the Atom name and a vector orthogonal to the support surface Al,1,1,1. When starting from a xyz file only the support atom is required
  -f, --folder <FOLDER>
          Output folder [default: ./sim/]
  -i, --iterations <ITERATIONS>
          Determines the length of the simulation
  -b, --begin-temperature <BEGIN_TEMPERATURE>
          Temperature where the annealing process starts [default: 5000]
  -t, --temperature <TEMPERATURE>
          Temperature at which the annealing process stops [default: 300]
  -o, --optimization-cut-off-fraction <OPTIMIZATION_CUT_OFF_FRACTION>
          Fraction of the simulation after which the annealing process is completed. After that, the temperature remains constant [default: 1 2]
      --support-e <SUPPORT_E>
          Support energy
      --e-l-cn <E_L_CN>
          File path or string containing a JSON-formatted list of energies
      --e-cn <E_CN>
          File path or string containing JSON-formatted energies
  -r, --repetition <REPETITION>
          How many times the same simulation is run. Multiple runs allow for convergence tests. The number will be part of the simulation folder name. After running `-r 0-1`, you can run `-r 1-2` and the previous simulation will not be overwritten [default: 0 1]
  -g, --grid-folder <GRID_FOLDER>
          Folder containing the setup files like neighbor sites. It can be created using the Python script `create_sites.py` [default: ../303030-pair]
  -x, --xyz-trajectory <XYZ_TRAJECTORY>
          Set how many snapshots are saved in each simulation. Snapshots are spread out equally throughout the simulation
      --heat-map
          Generate a heat map
  -h, --help
          Print help
  -V, --version
          Print version
```

## License

Licensed under GNU Affero General Public License. 

Cite This: https://pubs.acs.org/doi/10.1021/acs.jpcc.5c02777


