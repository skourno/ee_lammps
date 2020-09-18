from typing  import Dict       # used to enforce the parent class typing
from lammps  import lammps
from simData import simData
from input   import input_data

import sys
# ----------------------------------------------------------
# Parent class of all test objects that should be inherited 
# and used as a template for all test particles. The list
# of functions displayed here are used by the main script 
# to execute operations on the test particles
# ----------------------------------------------------------
class test_object():
	testPartType = ''

	def __init__(self):
		sys.exit('test_object.__init__ : ERROR - Illegal call. This class is a template for other child classes')
	def subEns_change(self):
		sys.exit('test_object.subEns_change : ERROR - Illegal call. subEns_change() is undefined within the class')
	def print_info(self):
		sys.exit('test_object.print_info : ERROR - Illegal call. print_info() is undefined within the class')

# ----------------------------------------------------------
# An ion pair is chosen as the test object. The charge is 
# gradually turned on/off when a subensemble change is
# attempted.
#
# When using a test ion pair the ions of the pair should
# appear last in the list of atom types.
# ----------------------------------------------------------
class test_IonPair(test_object):
	idxCat    = -1     # global index of the test cation
	idxAn     = -1     # global index of the test anion
	Dcharge   =  0.0   # charge parturbation between sub-ensembles
	charge    =  0.0   # currect charge of the test particles
	NIonPairs =  0     # number of ion pairs (when charge is full some other ions are chosen instead of the current pair)

	def __init__(self, nameCation: str,\
		                 iCatType:   int,\
		                 nameAnion:  str,\
		                 iAnType:    int,\
		                 sim:    simData,\
		                 lmp:     lammps):

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

		# Choose the initial test ion pair
		self.idxCat   = sim.NAtoms - 2*self.NIonPairs + 1
		self.idxAn    = sim.NAtoms -   self.NIonPairs + 1
		flag          = lmp.set_variable("idx_testCat", self.idxCat)
		flag          = lmp.set_variable("idx_testAn" , self.idxAn )

		self.Dcharge  = sim.EEHist.width_bin # inherited from the ee histo
		self.charge   = sim.EEHist.max       # start with full charge
	
		# set a string flag that denotes the type the test particle
		self.testPartType  = 'Ion pair' 