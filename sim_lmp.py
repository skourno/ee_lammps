from lammps import lammps

import sys
import numpy as np
	
def start_lammps(ISeed):
	"""
	Initialize a lammps simulation. 

	Returns a lammps simulation object that can be 
	passed to other routines.
	"""
	lmp = lammps(cmdargs=["-screen","none"])
	lmp.command("log            log.lammps.ee")
	lmp.command("units          real")
	lmp.command("atom_style     full")

	lmp.command("variable seed  string 11111")
	flag = lmp.set_variable("seed",ISeed)

	return lmp

# ----------------------------------------------------------------

# ----------------------------------------------------------------
def setup_tip4p_with_Ions(lmp, Temp):
	"""
	Initialize a simulation of tip4p/05 water with ions in lammps

	Temp : simulation temperature in K
	"""

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


# ----------------------------------------------------------------

# ----------------------------------------------------------------
def setup_LJ_with_Ions(lmp, Temp):
	"""
	Initialize a simulation of LJ water with ions in lammps

	Temp : simulation temperature in K
	"""

	# translate the arguments to lammps variables
	lmp.command("variable Text equal  300.0")
	flag = lmp.set_variable("Text" ,Temp)

	lmp.command("include        CGLJ.ff               ")  # the forcefield goes here

	lmp.command("velocity       all create ${Text} 1234  ")
	
	lmp.command("thermo_style   custom step temp press epair evdwl ecoul elong")
	lmp.command("thermo_modify  flush yes")
	lmp.command("thermo 500") 
	
	lmp.command("timestep 2.0")
	
	lmp.command("fix integrate all nvt temp ${Text} ${Text} 100.0")



#------------------------------------------------------------
def run_lammps_sim(lmp,NStepsRun):
	"""
	Run a simulation with length 'NStepsRun' in number of timesteps
	"""
	lmp.command("variable NStepsRun string 100")
	flag = lmp.set_variable("NStepsRun",NStepsRun)
	lmp.command("run ${NStepsRun}")

#------------------------------------------------------------
def set_lammps_dump(lmp,NStepDump=200,DumpFileName="expEns.lammpstrj"):
	"""
	Forces a lammps simulation to output dump configurations
	
	Argument list:
	lmp          -  the lammps instance
	NStepDump    -  write a config every this many steps (OPTIONAL)
	DumpFileName -  The name of the dump file            (OPTIONAL)
	"""
	lmp.command("variable DumpFileName  string dump.atom")
	lmp.command("variable NStepDump     string 200      ")

	flag = lmp.set_variable("DumpFileName",DumpFileName)
	flag = lmp.set_variable("NStepDump",   NStepDump)

	lmp.command("dump trj all atom ${NStepDump} ${DumpFileName}")


