from lammps         import lammps
from sim_lmp        import Initialize_tip4p_with_Ions_simulation, run_lammps_sim
from mpi4py         import MPI

comm = MPI.COMM_WORLD

lmp = Initialize_tip4p_with_Ions_simulation('tip4p05.data',298.15,98362)

lmp.command("dump myDump all xyz 500 dump.atom")
#lmp.command("dump_modify myDump scale no")

lmp.command("compute myRDF all rdf 50 3 4")
lmp.command("fix fmyRDF all ave/time 1000 1000 1000000 c_myRDF[*] file NaCl.rdf mode vector")

run_lammps_sim(lmp,2000000)