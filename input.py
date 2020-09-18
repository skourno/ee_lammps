import sys
import numpy as np

from Hist import Histogram
# -------------------------------------------------------------
# A class to read the input file of the program and store the
# relevant data.
#
# All the routines used to read input entries return an errorMessage.
# If there is no error while executing the routine then the 
# message will simply be a 'No error' string
#
# -------------------------------------------------------------
class input_data:

	# boolean flags that are changed to true if the corresponding entry was found in the input
	ee_histo_parsed     = False
	use_wl_parsed       = False
	use_tmmc_parsed     = False
	iseed_parsed        = False
	read_data_parsed    = False
	sim_time_parsed     = False
	tp_IonPair_parsed   = False
	set_temp_parsed     = False
	init_dos_parsed     = False
	write_wl_parsed     = False
	write_tmmc_parsed   = False
	write_dump_parsed   = False
	roam_ee_with_parsed = False

	# info that can be read from the input file 
	# (in comment the corresponding keyword of the input)
	boundary_min        = 0.0    # ee_histo
	boundary_max        = 0.0    # ee_histo
	desired_inc         = 0.0    # ee_histo
	lnf                 = 0.0    # use_wl
	lnf_scaler          = 0.0    # use_wl 
	ratio_crit          = 0.0    # use_wl    
	lnf_crit            = 0.0    # use_wl  
	use_wl_bool         = False  # use_wl
	use_tmmc_bool       = False  # use_tmmc
	iseed               = 0      # iseed
	DataFile            = ''     # read_data
	time_equil          = 0.0    # sim_time    
	time_prod           = 0.0    # sim_time
	time_sub_sim        = 0.0    # sim_time
	cationName          = ''     # tp_IonPair
	anionName           = ''     # tp_IonPair
	iTypeTestCat        = 0      # tp_IonPair
	iTypeTestAn         = 0      # tp_IonPair
	Temp                = 298.15 # set_temp
	outFile_wl          = ''     # write_Wl
	wstep_wl            = 0      # write_Wl
	outFile_tmmc        = ''     # write_tmmc
	wstep_tmmc          = 0      # write_tmmc
	outFile_dump        = ''     # write_dump
	wstep_dump          = 0      # write_dump
	ee_method           = ''     # roam_ee_with

	NoErrorMessage      = 'No error' # no error flag

	# ----------------------------------------------------------------
	# Initialize function
	# ----------------------------------------------------------------
	def __init__(self,InputFile):

		NLinesRead = 0

		with open(InputFile) as input:
			for line in input:
				NLinesRead += 1
				if (line == "\n"):
					continue # empty line


				# ignore anything that follows after a '#' as a comment 
				line        = line.partition('#')[0]
				if (len(line) == 0):
					continue

				line        = line.rstrip()
				lineArgs    = line.split()

				if  (lineArgs[0] == 'ee_histo' ):
					errorMessage           = self.ee_histo(lineArgs)
					self.ee_histo_parsed   = True
				elif(lineArgs[0] == 'use_wl'   ):
					errorMessage           = self.use_wl(lineArgs)
					self.use_wl_parsed     = True
				elif(lineArgs[0] == 'use_tmmc' ):
					self.use_tmmc_parsed   = True
				elif(lineArgs[0] == 'iseed'    ):
					errorMessage           = self.read_iseed(lineArgs)
					self.iseed_parsed      = True
				elif(lineArgs[0] == 'read_data'):
					errorMessage           = self.read_data(lineArgs)
					self.read_data_parsed  = True
				elif(lineArgs[0] == 'sim_time' ):
					errorMessage           = self.sim_time(lineArgs)
					self.sim_time_parsed   = True				
				elif(lineArgs[0] == 'tp_IonPair'):
					errorMessage           = self.tp_IonPair(lineArgs)
					self.tp_IonPair_parsed = True		
				elif(lineArgs[0] == 'set_temp'):
					errorMessage           = self.set_temp(lineArgs)
					self.set_temp_parsed   = True		
				elif(lineArgs[0] == 'init_dos'):
					errorMessage           = self.init_dos(lineArgs)
					self.init_dos_parsed   = True		
				elif(lineArgs[0] == 'write_wl'    ):
					errorMessage           = self.write_wl(lineArgs)
					self.write_wl_parsed   = True		
				elif(lineArgs[0] == 'write_tmmc'  ):
					errorMessage           = self.write_tmmc(lineArgs)
					self.write_tmmc_parsed = True
				elif(lineArgs[0] == 'write_dump'  ):
					errorMessage           = self.write_dump(lineArgs)
					self.write_dump_parsed = True
				elif(lineArgs[0] == 'roam_ee_with'):
					errorMessage             = self.roam_ee_with(lineArgs)
					self.roam_ee_with_parsed = True
				else:
					print('input_data.__init__ : ERROR - Unrecognised input file entry at line %d\n' %NLinesRead)
					print('                      %s is not recognised\n' %lineArgs[0]    )
					sys.exit()     

				if (errorMessage != self.NoErrorMessage):
					print('input_data.__init__ : ERROR - Erroneous entry at line %d\n' %NLinesRead)
					print(' > Printing ERROR info ...\n')
					print(' > %s\n' %errorMessage)

					sys.exit()

  # ----------------------------------
	def ee_histo(self,lineArgs):
		if (len(lineArgs) != 4):
			return 'ee_histo: Argument mismatch. Specify: min,max,desired_inc' # error
		self.boundary_min = float(lineArgs[1])
		self.boundary_max = float(lineArgs[2])
		self.desired_inc  = float(lineArgs[3])

		return self.NoErrorMessage

	# ----------------------------------
	def use_wl(self,lineArgs):
		if   (len(lineArgs) < 2):
			return 'use_wl: Called without arguments. Please choose yes/no' # error
		elif (not self.ee_histo_parsed):
			return 'use_wl: Cant find an EE histogram. Define one using ee_histo' # error

		if  (lineArgs[1] == 'yes'):
			if (len(lineArgs) != 6):
				return 'use_wl: Argument mismatch. Specify(after yes): lnf, lnf_scaler, ratio_crit, lnf_crit' # error

			self.use_wl      = True
			self.lnf         = float(lineArgs[2])
			self.lnf_scaler  = float(lineArgs[3])  
			self.ratio_crit  = float(lineArgs[4])  
			self.lnf_crit    = float(lineArgs[5]) 
		elif(lineArgs[1] == 'no' ):
			self.use_wl = False

		return self.NoErrorMessage 

	# ----------------------------------
	def use_tmmc(self,lineArgs):
		if   (len(lineArgs) != 2):
			return 'use_tmmc: Argument mismatch. Please only choose yes/no' # error
		elif (not self.ee_histo_parsed):
			return 'use_tmmc: Cant find an EE histogram. Define one using ee_histo' # error	
		
		self.use_tmmc = False
		if (lineArgs[1] == 'yes'):
			self.use_tmmc = True

		return self.NoErrorMessage


	# ----------------------------------
	def read_iseed(self,lineArgs):
		if (len(lineArgs) != 2):
			return 'iseed: Argument mismatch. Only specify the desired seed' # error
		self.iseed      = int(lineArgs[1])
		return self.NoErrorMessage

	def read_data(self,lineArgs):
		if (len(lineArgs) != 2):
			return 'read_data: Argument mismatch. Only specify the data file name'

		self.DataFile = str(lineArgs[1])
		return self.NoErrorMessage

	# ----------------------------------
	def sim_time(self,lineArgs):
		if  (len(lineArgs) != 4):
			return 'sim_time: Argument mismatch. Specify: time_equil, time_prod, time_sub_sim' # error
		elif(float(lineArgs[1]) < 0 or float(lineArgs[2]) < 0 or float(lineArgs[3]) < 0):
			return 'sim time: Simulation time cannot be negative' # error
		self.time_equil   = int(lineArgs[1])
		self.time_prod    = int(lineArgs[2])
		self.time_sub_sim = int(lineArgs[3])
		return self.NoErrorMessage

	# ----------------------------------
	def tp_IonPair(self,lineArgs):
		if  (len(lineArgs) != 5):
			return 'tp_IonPair: Argument mismatch. Specify: cationName, iTypeCat, anionName, iTypeAn' # error
		self.iTypeTestCat = int(lineArgs[1])
		self.cationName   = str(lineArgs[2])
		self.iTypeTestAn  = int(lineArgs[3])
		self.anionName    = str(lineArgs[4])
		return self.NoErrorMessage

	# ----------------------------------
	def set_temp(self,lineArgs):
		if  (len(lineArgs) != 2):
			return 'set_temp: Argument mismatch. Only specify the simulation temperature' # error
		elif(float(lineArgs[1]) < 0):
			return 'set_temp: Negative temperature'  # error
		self.Temp = float(lineArgs[1])
		return self.NoErrorMessage

	# ----------------------------------
	def init_dos(self,lineArgs):
		if   (len(lineArgs) != 2):
			return 'init_dos: Argument mismatch. Specify just the file path and file name of the initial DOS' # error
		elif (not self.ee_histo_parsed):
			return 'init_dos: Cant find an EE histogram. Define one using ee_histo' # error

		EEHisto_aux    = Histogram(self.boundary_min, self.boundary_max, self.desired_inc)
		NSubs          = EEHisto_aux.NBins
		self.dos_init  = np.zeros(NSubs, np.double)
		dos_init_check = np.zeros(NSubs, dtype=bool)
		iEntry         = -1

		with open(str(lineArgs[1])) as dosFile:
			for line in dosFile:
				if (line == "\n"):
					continue # empty line
		
				# ignore anything that follows after a '#' as a comment 
				line        = line.partition('#')[0]
				if (len(line) == 0):
					continue

				line        = line.rstrip()
				lineArgsDos = line.split()

				if (iEntry == EEHisto_aux.NBins-1):
					dosFile.close()
					return 'init_dos: The dos_init provided has more sub-ensembles that ee_histo' # error

				# Dos should be located in the third column
				iEntry                += 1
				self.dos_init[iEntry]	 = lineArgsDos[2]		
				dos_init_check[iEntry] = True

		if (not all(dos_init_check)):
			return 'init_dos: The dos_init provided has less sub-ensembles that ee_histo' # error

		return self.NoErrorMessage

	# ----------------------------------
	def write_wl(self,lineArgs):
		if   (len(lineArgs) != 3):
			return 'write_wl: Argument mismatch. Specify: write_step(fs), fileName' # error
		elif (not self.use_wl_parsed):
			return 'write_wl: Cant find a WL histogram. Define one using use_wl' # error

		self.wstep_wl   = lineArgs[1]
		self.outFile_wl = lineArgs[2]

		return self.NoErrorMessage


	# ----------------------------------
	def write_tmmc(self,lineArgs):
		if   (len(lineArgs) != 3):
			return 'write_tmmc: Argument mismatch. Specify: write_step(fs), fileName' # error
		elif (not self.use_tmmc_parsed):
			return 'write_tmmc: Cant find a TMMC histogram. Define one using use_tmmc' # error

		self.wstep_tmmc   = lineArgs[1]
		self.outFile_tmmc = lineArgs[2]

		return self.NoErrorMessage

	# ----------------------------------
	def write_dump(self,lineArgs):
		if   (len(lineArgs) != 3):
			return 'write_dump: Argument mismatch. Specify: write_step(fs), fileName' # error

		self.wstep_dump   = lineArgs[1]
		self.outFile_dump = lineArgs[2]

		return self.NoErrorMessage

	# ----------------------------------
	def roam_ee_with(self,lineArgs):
		if   (len(lineArgs) != 2):
			return 'roam_ee_with: Argument mismatch. Specify how to roam the ee' # error

		use_wl_parsed   = self.use_wl_parsed
		use_tmmc_parsed	= self.use_tmmc_parsed

		ee_method       = str(lineArgs[1])

		if (not ee_method == 'wl'  ) and (not ee_method == 'tmmc') and (not ee_method == 'wl_and_tmmc'):
			return 'roam_ee_with: Supported options wl/tmmc/wl_and_tmmc'

		if  (ee_method == 'wl'  ) and (not use_wl_parsed):
			return 'roam_ee_with: Cant find a WL histogram. Define one using use_wl' # error
		elif(ee_method == 'tmmc') and (not use_tmmc_parsed):
			return 'roam_ee_with: Cant find a TMMC histogram. Define one using use_tmmc' # error
		elif(ee_method == 'wl_and_tmmc'):
			if (not use_wl_parsed) or (not use_tmmc_parsed):
				return 'roam_ee_with: Cant find WL or TMMC histogram. Define both using use_wl and use_tmmc' # error

		self.ee_method = ee_method

		return self.NoErrorMessage









                


				
