import sys

from Hist    import Histogram
from input   import input_data
from typing  import Dict       # used to enforce the parent class typing
from sim_lmp import start_lammps
from lammps  import lammps
from tmmc    import TMMC_histogram
from WL      import WL_histogram

import numpy as np

# ----------------------------------------------
# ----------------------------------------------
class simData():
	"""
	This class stores all the information on the 
	submitted simulation run. 
	"""

	Temp           = 298.0 # Simulation temperature in K
	Beta           = 0.0   # Thermodynamic beta (1/kT)
	EEHist         = 0     # Exp Ens histogram
	WLHist         = 0     # Wang-Landau histogram
	TMHist         = 0     # Transition Matrix Monte Carlo histogram
	NStepsUpdateTM = 0     # Update the TM every this many steps
	write_wl       = False # output the WL histogram if true
	WL_FILE_out    = ''    # file path and name to the WL output
	write_tmmc     = False # output the TMMC histogram if true
	TM_FILE_out    = ''    # file path and name to the TMMC output
	ISeed          = 0     # seed to be used for lammps initialization
	NSteps_equil   = 0     # equilibration run length in timesteps
	NSteps_prod    = 0     # production run length in timesteps
	NSteps_subEns  = 0     # simulation length in timesteps for simulation runs between sub-change attempts
	use_wl_bool    = False # True if wl is used
	use_tmmc_bool  = False # True if tmmc is used 
	DataFileName   = ''    # The data file name and path. Contains the init config and topology
	NWStepWL       = 0     # wl write step
	NWStepTM       = 0     # tmmc write step
	wts            = 0.0   # weights used to roam the expanded ensemble

	simInit_bool   = False # True if lammps has been initialized
	sysLoaded_bool = False # True if a configuration and topology has been read

	# configuration info (imported using import_config() )
	NAtoms         = 0   # Number of atoms in the simulation

	# ----------------------------------------------
	# ----------------------------------------------
	def __init__(self,inData: input_data,write_info_bool=False):
		"""
		You can initialize with an input_data object.
		
		OPTIONAL arg : write_info_bool, if True then
		at the end of execution the imported info will 
		be printed
		"""
		self.Temp          = inData.Temp                           # Simulation temperature in K
		self.Beta          = 1.0 / (self.Temp * 8.314 / 4184.0)    # Thermodynamic beta (1/kT)
	
		# General EE histogram info
		if (not inData.ee_histo_parsed):
			sys.exit('simulation.__init__ : ERROR - Exp Ens histogram data not found')
		self.EEHist  = Histogram(inData.boundary_min, inData.boundary_max, inData.desired_inc)

		if (inData.init_dos_parsed):
			wts_init = inData.dos_init
		else:
			# initialize to zero array
			wts_init = np.zeros(self.EEHist.NBins,np.double)

		# set the initial weights of the simulation
		self.EEHist.binValue = np.copy(wts_init)

		if (inData.use_wl_bool):
			lnf_init      = inData.lnf
			lnf_scaler    = inData.lnf_scaler
			ratio_crit    = inData.ratio_crit
			lnf_crit      = inData.lnf_crit
			self.WLHist   = WL_histogram(self.EEHist,lnf_init,lnf_scaler,ratio_crit,lnf_crit)
			self.use_wl_bool = True
		
		if (inData.use_tmmc_bool):
			self.TMHist         = TMMC_histogram(self.EEHist)
			self.use_tmmc_bool  = True
			self.NStepsUpdateTM = inData.NStepsUpdateTM
		
		if (inData.write_wl_parsed):
			self.write_wl     = True
			self.WL_FILE_out  = open(inData.outFile_wl,"w") # File to write DOS from WL
			self.NWStepWL     = inData.wstep_wl
	
		if (inData.write_tmmc_parsed):
			self.write_tmmc    = True
			self.TM_FILE_out   = open(inData.outFile_tmmc,"w") # File to write DOS from TMMC
			self.NWStepTM      = inData.wstep_tmmc

		if (not inData.iseed_parsed):
			self.ISeed = random.randint(1000,1000000) # get a random seed between 1000 and 1000000
		else:
			self.ISeed = inData.iseed

		if (not inData.sim_steps_parsed):
			sys.exit('simulation.__init__ : ERROR - Simulation number of steps variables are not set')
		self.NSteps_equil   = inData.NSteps_equil 
		self.NSteps_prod    = inData.NSteps_prod
		self.NSteps_subEns  = inData.NSteps_subEns 

		if (inData.read_data_parsed):
			self.DataFileName = inData.DataFile
		else:
			sys.exit('simulation.__init__ : ERROR - Use read_data in the input to import a configuration and topology')

		self.wts = np.zeros(self.EEHist.NBins,np.double)
		if  (inData.ee_method == 'wl' or inData.ee_method == 'wl_and_tmmc'):
			self.update_roaming_weights('wl')
		elif(inData.ee_method == 'tmmc'):
			self.update_roaming_weights('tmmc')
		else:
			sys.exit('simulation.__init__ : ERROR - Undefined roaming method. Use ee_method to define one')

		# print the initialized Sys info
		if(write_info_bool):
			#self.print_info
			pass

	# ----------------------------------------------
	# ----------------------------------------------
	def init_lammps_Sim(self):
		"""
		Initializes a lammps simulation using the data 
		stored in the class
		"""
		lmp          = start_lammps(self.ISeed)
		simInit_bool = True
		return lmp 

	# -----------------------------------------------
	# -----------------------------------------------
	def import_config(self,lmp):
		"""
		Reads the data file and initializes the associated
		class variables.
		"""
		lmp.command("variable DataFileName string data.run")
		flag = lmp.set_variable("DataFileName",self.DataFileName)
		lmp.command("read_data ${DataFileName} ")

		lmp.command("variable NAtoms   equal count(all)")

		self.NAtoms         = int(lmp.extract_variable("NAtoms","group",0))
		self.sysLoaded_bool = True


	# ----------------------------------------------
	# ----------------------------------------------
	def update_roaming_weights(self,ee_method):
		"""
		Updates the simulation weights acoording to the ee_method
		specified
		"""
		if  (ee_method == 'wl'):
			self.wts = np.copy(self.WLHist.wts)
		elif(ee_method == 'tmmc'):
			self.wts = np.copy(self.TMHist.wts)
		elif(ee_method == 'wl_and_tmmc'):
			if (self.WLHist.isItDone and self.TMHist.allSubsLogged):
				# then it is time to roam with the TMMC weights
				self.wts              = np.copy(self.TMHist.wts)
				self.use_wl_bool      = False
				self.TMHist.activated = True
			else:
				# if wl is not done or tmmc is not ready to produce weights
				# then we use wl weights
				self.wts = np.copy(self.WLHist.wts)




