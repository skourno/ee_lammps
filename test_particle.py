from typing  import Dict       # used to enforce the parent class typing
from lammps  import lammps
from simData import simData
from input   import input_data
from mpi4py  import MPI

import sys
import numpy as np

# ----------------------------------------------------------
class test_object():
	"""
	Parent class of all test objects that should be inherited 
	and used as a template for all test particles. The list
	of functions displayed here are used by the main script 
	to execute operations on the test particles
	"""
	Type = 'test particle template'

	def __init__(self):
		sys.exit('test_object.__init__ : ERROR - Illegal call. This class is a template for other child classes')
	def subEns_change(self):
		sys.exit('test_object.subEns_change : ERROR - Illegal call. subEns_change() is undefined within the class')
	def ee_coord(self):
		sys.exit('test_object.subEns_change : ERROR - Illegal call. ee_coord() is undefined within the class')
	def print_idx(self):
		sys.exit('test_object.print_idx : ERROR - Illegal call. print_idx() is undefined within the class')



# ----------------------------------------------------------
class test_IonPair(test_object):
	"""
	An ion pair is chosen as the test object. The charge is 
	gradually turned on/off when a subensemble change is
	attempted.
	When using a test ion pair the ions of the pair should
	appear last in the list of atom types. 
	"""

	idxCat     = -1     # global index of the test cation
	idxAn      = -1     # global index of the test anion
	Dcharge    =  0.0   # charge parturbation between sub-ensembles
	charge     =  0.0   # currect charge of the test particles
	NIonPairs  =  0     # number of ion pairs (when charge is full some other ions are chosen instead of the current pair)
	chargeMin  =  0.0   # is the minimum charge of the ee simulation
	chargeMax  =  1.0   # is the maximum charge of the ee simulation
	fullCharge =  1.0   # set by default to be 1.0
	idxCat1    =  0     # the index of the first Cation (of the test part type)
	idxAn1     =  0     # the index of the first Anion  (of the test part type)
	Type       = 'Ion Pair' 

	# ----------------------------------------------------------
	def __init__(self,   nameCation: str,\
		                 iCatType:   int,\
		                 nameAnion:  str,\
		                 iAnType:    int,\
		                 sim:    simData,\
		                 lmp:     lammps):
		"""
		Initialize the a test Ion Pair using:
			nameCation : string containing the Cation name
			iCatType   : index of the atom type that is the Cation
			nameAnion  : string containing the Anion name
			iAnType    : index of the atom type that is the Anion 
			sim        : an initialized simulation object
			lmp        : an initialized lammps simulation
		"""

		# Initializing auxiliary variables within lammps
		lmp.command("variable idx_testCat  string -1")
		lmp.command("variable idx_testAn   string -1")
		lmp.command("variable q_testCat    string  0")
		lmp.command("variable q_testAn     string  0")
		lmp.command("variable nameCat      string  cation")
		lmp.command("variable nameAn       string  anion" )
		lmp.command("variable iCat         string  0")
		lmp.command("variable iAn          string  0")
	
		flag          = lmp.set_variable("nameCat" , nameCation)
		flag          = lmp.set_variable("nameAn"  , nameAnion )
		flag          = lmp.set_variable("iCat"    , iCatType  )
		flag          = lmp.set_variable("iAn"     , iAnType   )

		# make auxiliary ion groups
		lmp.command("group ${nameCat} type ${iCat}")
		lmp.command("group ${nameAn}  type ${iAn}" )

		lmp.command("variable NCations equal count(${nameCat})")

		print("test_IonPair.__init__ : WARNING - Assuming Cation is the next to last type (NAtTypes-1) and anion last (NAtTypes)")
		self.NIonPairs  = int(lmp.extract_variable("NCations","group",0))

		if (sim.NAtoms <= 0):
			print("test_IonPair.__init__ : ERROR - Found illegal num of Atoms in simData %d\n" %sim.NAtoms )
			print('                                Load a configuration before Initializing a test_IonPair')	

		# calculate the first Anion and Cation indices. Assuming they appear last after all atoms
		self.idxCat1 = sim.NAtoms - 2*self.NIonPairs + 1
		self.idxAn1  = sim.NAtoms -   self.NIonPairs + 1

		# Choose the initial test ion pair to be idxCat1 and idxAn1
		self.idxCat   = self.idxCat1
		self.idxAn    = self.idxAn1
		flag1          = lmp.set_variable("idx_testCat", self.idxCat)
		flag2          = lmp.set_variable("idx_testAn" , self.idxAn )

		self.Dcharge   = sim.EEHist.width_bin # inherited from the ee histo
		self.charge    = sim.EEHist.max       # start with full charge
		self.chargeMin = sim.EEHist.min
		self.chargeMax = sim.EEHist.max
		self.fullCharge = sim.EEHist.max
	
		# set a string flag that denotes the type the test particle
		self.Type     = 'Ion Pair' 

	#-----------------------------------------------
	def ee_coord(self):
		"""
		The exp ens coordinate of the test Ion pair 
		is the abs(charge)
		"""

		return self.charge

	#-----------------------------------------------
	def subEns_change(self,lmp,iDir):
		"""
		Changes the sub Ens to the specified direction iDir
		lmp  : is an initialized lammps simulation
		iDir : is the direction of the change
		       -1 decreases the sub index
		        0 is a remain move
		       +1 increases the sub index
		"""
		delta_q      = iDir * self.Dcharge
		q_testCat    = self.charge + delta_q
		q_testAn     = - q_testCat
		self.charge  = + q_testCat

		if (q_testCat == 0.0):
			q_testCat = 1.e-10
			q_testAn  = 1.e-10

		flag = lmp.set_variable("q_testCat", q_testCat)
		flag = lmp.set_variable("q_testAn" , q_testAn)

		lmp.command("set atom ${idx_testCat} charge ${q_testCat}")
		lmp.command("set atom ${idx_testAn}  charge ${q_testAn}" )

	#-----------------------------------------------
	def print_idx(self):
		"""Prints the indices of the test ion pair"""
		return "C%d A%d" %(self.idxCat, self.idxAn)


	#-----------------------------------------------
	def shuffle_testPart(self,lmp,comm):
		"""
		If all ion pairs have full fractional charges, 
		randomly select ions to be the test pair.

		lmp  : an initialized lammps simulation
		comm : is the mpi instance. Used for parallel
		execution
		"""

		if (comm.Get_rank() == 0):
			NIonPairs     = self.NIonPairs
			idx_shift_cat = np.random.randint(0,NIonPairs)
			idx_shift_an  = np.random.randint(0,NIonPairs)

			self.idxCat   = self.idxCat1 + idx_shift_cat 
			self.idxAn    = self.idxAn1  + idx_shift_an
	
		self.idxCat = comm.bcast(self.idxCat,root=0)
		self.idxAn  = comm.bcast(self.idxAn ,root=0)
		flag        = lmp.set_variable("idx_testCat",self.idxCat)
		flag        = lmp.set_variable("idx_testAn", self.idxAn)

# ----------------------------------------------------------
class test_Ion(test_object):
	"""
	An ion is chosen as the test object. The charge is 
	gradually turned on/off when a subensemble change is
	attempted.
	When using a test ion, it should
	appear last in the list of atom types. 
	"""

	idx        = -1     # global index of the test cation
	Dcharge    =  0.0   # charge parturbation between sub-ensembles
	charge     =  0.0   # currect charge of the test particle
	chargeMin  =  0.0   # is the minimum charge of the ee simulation (absolute)
	chargeMax  =  1.0   # is the maximum charge of the ee simulation (absolute)
	fullCharge =  1.0   # set by default to be 1.0                   (absolute)
	Type       = 'Ion' 
	sign       =  0     # +1 for cation and -1 for anion

	# ----------------------------------------------------------
	def __init__(self,   name:    str,\
		                 iType:   int,\
		                 charge:  float,\
		                 sim:     simData,\
		                 lmp:     lammps):
		"""
		Initialize the a test Ion Pair using:
			name       : string containing the ion name
			iType      : index of the atom type that is the ion
			sim        : an initialized simulation object
			lmp        : an initialized lammps simulation
		"""

		# Initializing auxiliary variables within lammps
		lmp.command("variable idx_testIon  string -1")
		lmp.command("variable q_testIon    string  0")
		lmp.command("variable nameIon      string  ion")
		lmp.command("variable iIon         string  0")
	
		flag          = lmp.set_variable("nameIon"    , name )
		flag          = lmp.set_variable("iIon"    , iType)

		# make auxiliary ion groups
		lmp.command("group ${nameIon} type ${iIon}")

		if (sim.NAtoms <= 0):
			print("test_Ion.__init__ : ERROR - Found illegal num of Atoms in simData %d\n" %sim.NAtoms )
			print('                            Load a configuration before Initializing a test_IonPair')	

		# the ion should always be the last atom
		self.idx = sim.NAtoms 

		# Choose the initial test ion pair to be idxCat1 and idxAn1
		flag1          = lmp.set_variable("idx_testIon", self.idx)

		self.Dcharge   = sim.EEHist.width_bin # inherited from the ee histo
		self.chargeMin = sim.EEHist.min
		self.chargeMax = sim.EEHist.max
		self.fullCharge = sim.EEHist.max
	
		# set a string flag that denotes the type of the test particle
		self.Type     = 'Ion'


		self.charge     = charge
		self.sign       = np.sign(charge)

	#-----------------------------------------------
	def ee_coord(self):
		"""
		The exp ens coordinate of the test Ion pair 
		is the abs(charge)
		"""

		return abs(self.charge)

	#-----------------------------------------------
	def subEns_change(self,lmp,iDir):
		"""
		Changes the sub Ens to the specified direction iDir
		lmp  : is an initialized lammps simulation
		iDir : is the direction of the change
		       -1 decreases the sub index
		        0 is a remain move
		       +1 increases the sub index
		"""


		delta_q        = iDir * self.Dcharge * self.sign
		self.charge    = self.charge + delta_q
		q_testIon      = self.charge

		if (self.charge  == 0.0):
			q_testIon  = 1.e-10 * self.sign

		flag = lmp.set_variable("q_testIon", q_testIon)

		lmp.command("set atom ${idx_testIon} charge ${q_testIon}")

	#-----------------------------------------------
	def print_idx(self):
		"""Prints the indices of the test ion pair"""
		return "Ion%d" %(self.idx)