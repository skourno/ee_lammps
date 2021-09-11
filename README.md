ee_lammps v1.1

This is a code designed to interface lammps and conduct expanded ensemble simulations. 
Example scripts can be found in the ee_lammps/examples folder.


Running expanded ensemble simulations
-------------------------------------------------------------------------------------------------

You can find an example script (ee_lammps.py) and an annotated sample input in examples/ee_run/.

The code supports expanded ensemble simulations to compute:
1) The electrostatic contribution to the chemical potential of ion pairs (salt molecules).
2) The electrostatic contribution to the chemical potential of charged particles.

More cases can be added by creating the corresponding child class of the test_object template class.
Specifically you need to:
1) Create the class that inherits the test_object class.
2) Create the constructor. Make sure you add all the info you need for the rest of the methods.
3) Define the .subEns_change() method. It should direct LAMMPS on how to change your sub-ensemble.
4) Define the .ee_coord() method. This should return the current value of your reaction coordinate/collective variable.
5) Define the .print_idx() method. This should return the index of the current sub-ensemble.
6) ... Add whatever optimization or specialized method you need.

Creating initial configurations
-------------------------------------------------------------------------------------------------

The code can be used to generate initial configurations on a lattice. There are suitable classes to declare 
molecules (Molecule), and atoms (Atom). There is an annotated example script in examples/init_config.



