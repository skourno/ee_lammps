from lammps import lammps
from input  import input_data
from typing import Dict       # used to enforce types of variables in the argument list
from Hist   import Histogram

import sys
import numpy as np

class simulation():
	def __init__(self, inData: input_data):
	
		# General EE histogram info
		if (not inData.ee_histo_parsed):
			sys.exit('simulation.__init__ : ERROR - Exp Ens histogram data not found')
		global EEHist  = Histogram(inData.boundary_min, inData.boundary_max, inData.desired_inc)
		
		if (inData.use_wl):
			lnf_init      = inData.lnf
			lnf_scaler    = inData.lnf_scaler
			ratio_crit    = inData.ratio_crit
			lnf_crit      = inData.lnf_crit
			global WLHist = WL_histogram(self.EEHist,lnf_init,lnf_scaler,ratio_crit,lnf_crit)
		
		if (inData.use_tmmc):
			global TMHist = TMMC_histogram(self.EEHist)
		
		if (inData.write_wl_parsed):
			global write_wl     = True
			global WL_FILE_out  = open(inData.outFile_wl,"w") # File to write DOS from WL
	
		if (inData.write_tmmc_parsed):
			global write_tmmc     = True
			global TM_FILE_out  = open(inData.outFile_tmmc,"w") # File to write DOS from TMMC	

		DatafileName = inData.DataFile
		# Needs to be passed as arguments through an input file (TODO)
		n_atoms_water = 6402
		i_atom_Na_i   = 6403
		i_atom_Na_f   = 6460
		i_atom_Cl_i   = 6461
		i_atom_Cl_f   = 6518
		natoms        = 6518
		NIonsPairs    = 58


		if (not inData.sim_time_parsed):
			sys.exit('simulation.__init__ : ERROR - Simulationt times are not set')
		self.time_equil   = inData.time_equil 
		self.time_prod    = inData.time_prod
		self.time_sub_sim = inData.time_sub_sim

		if (inData.init_dos_parsed):
			wts_init = inData.dos_init
		else:
			# initialize to zero array
			wts_init = np.zeros(EEHist.NBins,np.double)

		global Temp  = inData.Temp
		global Beta  = 1.0 / (Temp * 8.314 / 4184.0)
		ISeed        = inData.iseed
		global lmp   = Initialize_tip4p_with_Ions_simulation(DatafileName, Temp, ISeed)
	


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

#------------------------------------------------------------
# run a simulation with length 'time' in fs
#------------------------------------------------------------
def run_lammps_sim(lmp,time):
	lmp.command("variable time string -1")
	flag = lmp.set_variable("time",time)
	lmp.command("run ${time}")