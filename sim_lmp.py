from lammps import lammps

import sys
import numpy as np
	
# ----------------------------------------------------------------
# Start lammps
# ----------------------------------------------------------------
def start_lammps(ISeed):
	lmp = lammps(cmdargs=["-screen","none"])
	lmp.command("units          real")
	lmp.command("atom_style     full")

	lmp.command("variable seed  string 11111")
	flag = lmp.set_variable("seed",ISeed)

	return lmp

# ----------------------------------------------------------------
# Initialize a simulation of tip4p/05 water with ions in lammps
#
# Temp         : simulation temperature
# ----------------------------------------------------------------
def setup_tip4p_with_Ions(lmp, Temp):
	# translate the arguments to lammps variables
	lmp.command("variable Text equal  298.0")
	flag = lmp.set_variable("Text" ,Temp)

	lmp.command("include        tip4p05.ff               ")  # the forcefield goes here

	lmp.command("velocity       all create ${Text} 1234  ")

	lmp.command("neighbor       2.0 bin                  ")
	lmp.command("neigh_modify   every 1 delay 0 check yes")
	
	lmp.command("thermo_style   custom step temp press epair evdwl ecoul elong")
	lmp.command("thermo_modify  flush yes")
	lmp.command("thermo 500") 
	
	lmp.command("timestep 2.0")
	
	lmp.command("fix constrain all shake 1.0e-4 100 0 b 1 a 1    ")
	lmp.command("fix integrate all nvt temp ${Text} ${Text} 100.0")
	lmp.command("fix removeMomentum all momentum 1 linear 1 1 1  ")


#------------------------------------------------------------
# run a simulation with length 'time' in number of timesteps
#------------------------------------------------------------
def run_lammps_sim(lmp,time):
	lmp.command("variable time string -1")
	flag = lmp.set_variable("time",time)
	lmp.command("run ${time}")