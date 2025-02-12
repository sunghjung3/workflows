import numpy as np
import copy

from ase.build import fcc111
from ase.constraints import FixAtoms
from ase.ga.data import PrepareDB
from ase.ga.startgenerator import StartGenerator
from ase.ga.utilities import closest_distances_generator, get_all_atom_types
from ase.io import read
from ase.atom import Atom
import ase.data
from ase.geometry import get_distances


def main(n_atoms, element):
    db_file = 'gadb.db'
    atom_comp = n_atoms * [element]  # composition of atoms to add
    atom_numbers = [ase.data.atomic_numbers[e] for e in atom_comp]

    # import the skeleton
    skeleton = read('../../POSCAR')

    # define the volume in which atoms are placed
    # the volume is defined by a corner position (p0)
    # and three spanning vectors (v1, v2, v3)
    cell = skeleton.get_cell()
    p0 = np.array([0., 0., 0.])
    v1 = cell[0, :]
    v2 = cell[1, :]
    v3 = cell[2, :]

    # define the closest distance two atoms of a given species can be to each other
    unique_atom_types = get_all_atom_types(skeleton, atom_numbers)
    blmin = closest_distances_generator(atom_numbers=unique_atom_types,
                                        ratio_of_covalent_radii=0.7)

    # create the database to store information in
    d = PrepareDB(db_file_name=db_file,
                  simulation_cell=skeleton,
                  stoichiometry=atom_numbers)

    # generate the starting population
    population_size = 20
    for i in range(population_size):  # each new structure
        a = copy.deepcopy(skeleton)
        for j, atom_number in enumerate(atom_numbers):  # each atom to add
            acceptable = False
            while not acceptable:
                pos = p0 + np.random.rand()*v1 + np.random.rand()*v2 + np.random.rand()*v3
                pos = pos[np.newaxis, :]  # make 2D
                for atom in a:
                    pos_atom = atom.position[np.newaxis, :]
                    _, dist = get_distances(pos, pos_atom, cell=a.cell, pbc=a.pbc)  # with PBC
                    ituple = (atom_number, atom.number)
                    if dist < blmin[ituple]:
                        acceptable = False
                        break
                    else:
                        acceptable = True
            a.append(Atom(atom_number, position=pos[0]))
        d.add_unrelaxed_candidate(a)


if __name__ == '__main__':
    import os

    element = 'Mg'
    for i in range(1, 13):
        dir_name = f"{element}-{i}"
        print(dir_name)
        os.mkdir(dir_name)  # Create the directory
        
        os.chdir(dir_name)  # Change into the new directory
        main(i, element)    # Run main() inside the new directory
        
        os.chdir('..')  # Return to the original directory
