from lammps         import lammps
from sim_lmp        import start_lammps, setup_tip4p_with_Ions, set_lammps_dump, run_lammps_sim
from mpi4py         import MPI

comm = MPI.COMM_WORLD

lmp = start_lammps(ISeed=752663)

lmp.command("read_data data.equil.final")

setup_tip4p_with_Ions(lmp,Temp=298.15)



set_lammps_dump(lmp,NStepDump=250,DumpFileName="rdf.lammpstrj")


lmp.command("compute myRDF all rdf 50 3 4")
lmp.command("fix fmyRDF  all ave/time 200 5000 1000000 c_myRDF[*] file rdf_NaCl.xvg         mode vector")
lmp.command("fix fmyRDF2 all ave/time 200 500  100000  c_myRDF[*] file rdf_NaCl_rolling.xvg mode vector")

run_lammps_sim(lmp,1000000)