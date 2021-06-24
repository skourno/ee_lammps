from data_lmp import data_lammps, data_to_molecule, Mol_System_to_data
from molSys   import Mol_System, Create_config_on_lattice
from lattice  import lattice
from molecule import molecule
from atom     import atom

import numpy as np
import sys

WaterDataFile = 'data.singleTIP4P-2005'
NWaters       =  2134
NCations      =  29
NAnions       =  29
DataFileOut   = 'data.test'


# read a datafile that contains the water molecule
DataWat  = data_lammps(WaterDataFile)

water = data_to_molecule(DataWat) # create a molecule object from the Data input
Cat   = atom(idx=0, iSpc=1, Type=2, charge=+1.0, xyz=np.zeros(3,np.double))
An    = atom(idx=0, iSpc=2, Type=3, charge=-1.0, xyz=np.zeros(3,np.double))

# what to build is an array of tuples that contains the molecule/atom 
# and the corresponding num of the species to be created
what_to_build    = [None] * 3
what_to_build[0] = (water, NWaters)
what_to_build[1] = (Cat,   NCations)
what_to_build[2] = (An,    NAnions)

NSpecies         = NWaters + NCations + NAnions

edgeUC  = [3.1]*3 
latt    = lattice('cub', edgeUC, NSpecies) # type of latt, edge lengths of the UC, minimum # of Nodes in the latt

# create a molecular system with the what_to_build directives, on the the lattice latt
molSys  = Create_config_on_lattice(latt,what_to_build)

# num of random Swaps, one species of the swaped pair is always in the IonIndices slice
IonIndices = slice(NWaters,-1,1)
molSys.rand_swaps(NSwaps=10000,TargetSlice=IonIndices) 


DataOut = Mol_System_to_data(molSys)
DataOut.write(DataFileOut)
