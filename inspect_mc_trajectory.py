#!/opt/homebrew/bin/python3
import sys
from ase.visualize import view
import argparse
from ase import Atoms, Atom
import csv


def read_atom_sites(atoms_file):
    index_xyz = []
    with open(atoms_file) as f:
        for line in f:
            xyz = []
            for coordinate in line.split():
                xyz.append(float(coordinate))
            index_xyz.append(xyz)
    return index_xyz


def read_snapshots(snap_file, index, index_xyz):
    snap_shots = []
    snap_shots.append([])

    with open(snap_file, newline="") as csvfile:
        csv.field_size_limit(sys.maxsize)
        csv_reader = csv.reader(csvfile, delimiter=",", quotechar="|")
        for row in csv_reader:
            snap_shots.append(row)

    atoms = []
    if len(index) == 0:
        index = [x for x in range(1, len(snap_shots))]
    for i in index:
        atom_list = []
        for i, x in enumerate(snap_shots[i]):
            if x == "0":
                continue
            atom_list.append(Atom(x, position=index_xyz[i]))

        atoms.append(Atoms(atom_list))
    return atoms


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "dir",
        help="folder of the simualtion which contains the snap_shot_sections.csv-file e.g. './sim/900K_15000000I_711A_0/'",
        default=None,
    )  # default set path later
    parser.add_argument(
        "-i",
        "--index_sections",
        help="-1 is last, can be multiple seperated by blank",
        nargs="+",
        default=[],
        type=int,
    )
    parser.add_argument(
        "-g",
        "--grid",
        help="basis for reading data, has to be the same file used in the simualtion",
        default="/Users/tilman/Documents/Vertiefer/cluster_simulations/202020-pair/",
    )

    exp_folder = parser.parse_args().dir
    index = parser.parse_args().index_sections
    grid_folder = parser.parse_args().grid

    atoms_file = grid_folder + "atom_sites"
    atom_positions = read_atom_sites(atoms_file)

    snap_file = exp_folder + "snapshot_sections.csv"
    atoms = read_snapshots(snap_file, index, atom_positions)

    print(snap_file)

    view(atoms)
