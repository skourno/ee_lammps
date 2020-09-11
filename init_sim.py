from lammps import lammps

# ----------------------------------------------------------------
# Initialize a simulation of tip4p/05 water with ions in lammps
#
# DatafileName : the name of the .data file that contains the initial
#                configuration
#
# ----------------------------------------------------------------
def Initialize_tip4p_with_Ions_simulation(DatafileName, Temp, ISeed):
	
	lmp = lammps(cmdargs=["-screen","none"])

	# translate the arguments to lammps variables
	lmp.command("variable DatafileName string data.run")
	lmp.command("variable Text         equal  298.0"   )
	lmp.command("variable seed         string 11111"   )
	flag = lmp.set_variable("DatafileName",DatafileName)
	flag = lmp.set_variable("Text"        ,Temp        )
	flag = lmp.set_variable("seed"        ,ISeed       )

	lmp.command("units          real                     ")
	lmp.command("atom_style     full                     ")
	lmp.command("read_data      ${DatafileName}          ")
	lmp.command("include        tip4p05.ff               ")  # the forcefield goes here

	lmp.command("velocity       all create ${Text} 1234  ")

	lmp.command("neighbor       2.0 bin                  ")
	lmp.command("neigh_modify   every 1 delay 0 check yes")
	
	lmp.command("thermo_style   custom step temp press epair evdwl ecoul elong")
	lmp.command("thermo_modify  flush yes")
	lmp.command("thermo 100") 
	
	lmp.command("timestep 2.0")
	
	lmp.command("fix constrain all shake 1.0e-4 100 0 b 1 a 1    ")
	lmp.command("fix integrate all nvt temp ${Text} ${Text} 100.0")
	lmp.command("fix removeMomentum all momentum 1 linear 1 1 1  ")

	return lmp
