from readUt    import clean_and_split
from atom      import atom
from molecule  import bond, angle, molecule

import numpy as np
import sys

class data_lammps:
	"""This class reads and stores data files from lammps"""
	NAtoms      = 0
	NBonds      = 0
	NAngles     = 0
	NAtomTypes  = 0
	NBondTypes  = 0
	NAngleTypes = 0

	Box_lo      = np.zeros(3,np.double)
	Box_hi      = np.zeros(3,np.double)

	At          = []
	Bnd         = []
	Ang         = []

	def __init__(self,DataFilePath):
		"""Initialize with just the file path to the lammps data file"""
		DataFile = open(DataFilePath,"r")

		self.__read_generalInfo(DataFile)

		while True:
			line = DataFile.readline()

			if not line:
				DataFile.close()
				break

			if (line == "\n"):
				continue # empty line

			# ignore anything that follows after a '#' as a comment 
			line        = line.partition('#')[0]

			# ignore empty lines
			if (len(line) == 0):
				continue
		
			line        = line.rstrip()
			lineArgs    = line.split()

			if   (lineArgs[0] == 'Atoms'):
				self.__read_Atoms(DataFile)
			elif (lineArgs[0] == 'Bonds'):
				self.__read_Bonds(DataFile)
			elif (lineArgs[0] == 'Angles'):
				self.__read_Angles(DataFile)
			else:
				sys.exit('data_lammps.__init__ : ERROR - Unexpected keyword %s in DataFile' %lineArgs[0])

				
	def __read_generalInfo(self,DataFile):
		line = DataFile.readline() # empty line
		line = DataFile.readline() # empty line

		line         = DataFile.readline()
		self.NAtoms  = int(clean_and_split(line)[0])

		line         = DataFile.readline()
		self.NBonds  = int(clean_and_split(line)[0])

		line         = DataFile.readline()
		self.NAngles = int(clean_and_split(line)[0])

		line = DataFile.readline() # empty line

		line             = DataFile.readline()
		self.NAtomTypes  = int(clean_and_split(line)[0])

		line             = DataFile.readline()
		self.NBondTypes  = int(clean_and_split(line)[0])

		line             = DataFile.readline()
		self.NAngleTypes = int(clean_and_split(line)[0])

		line = DataFile.readline() # empty line

		line             = DataFile.readline()
		lineArgs         = clean_and_split(line)
		self.Box_lo[0]   = float(lineArgs[0])
		self.Box_hi[0]   = float(lineArgs[1])

		line             = DataFile.readline()
		lineArgs         = clean_and_split(line)
		self.Box_lo[1]   = float(lineArgs[0])
		self.Box_hi[1]   = float(lineArgs[1])

		line             = DataFile.readline()
		lineArgs         = clean_and_split(line)
		self.Box_lo[2]   = float(lineArgs[0])
		self.Box_hi[2]   = float(lineArgs[1])	


	def __read_Atoms(self,DataFile):
		line    = DataFile.readline() # empty line
		self.At = [None] * self.NAtoms

		for iAt in range(0,self.NAtoms):
			line     = DataFile.readline()
			lineArgs = clean_and_split(line)

			idx      = int(lineArgs[0])-1
			iSpc     = int(lineArgs[1])-1
			Type     = int(lineArgs[2])-1
			charge   = float(lineArgs[3])
			xyz      = [float(r) for r in lineArgs[4:]] 

			AtomIn   = atom(idx, iSpc, Type, charge, xyz)

			self.At.insert(idx, AtomIn)


	def __read_Bonds(self,DataFile):
		line     = DataFile.readline() # empty line
		self.Bnd = [None] * self.NBonds

		for iBnd in range(0,self.NBonds):
			line     = DataFile.readline()
			lineArgs = clean_and_split(line)

			idx      = int(lineArgs[0])-1
			Type     = int(lineArgs[1])-1
			Atoms    = [int(i) for i in lineArgs[2:]] 

			BondIn   = bond(idx, Type, Atoms)

			self.Bnd.insert(idx, BondIn)


	def __read_Angles(self,DataFile):
		line     = DataFile.readline() # empty line
		self.Ang = [None] * self.NAngles

		for iAng in range(0,self.NAngles):
			line     = DataFile.readline()
			lineArgs = clean_and_split(line)

			idx      = int(lineArgs[0])-1
			Type     = int(lineArgs[1])-1
			Atoms    = [int(i) for i in lineArgs[2:]] 

			AngleIn  = angle(idx, Type, Atoms)

			self.Ang.insert(AngleIn.idx, AngleIn)


def data_to_molecule(data: data_lammps):
	"""
	Creates a molecule provided a data_lammps object.
	Note that the the data_lammps must only contain
	one species. In the case of an atom the function 
	will return an error.
	"""
	if (not data.NBonds > 1):
		sys.exit('data_to_molecule : ERROR - This function only works for data with bonds')

	Atoms = [None] * data.NAtoms
	iSpc  = data.At[0].iSpc

	for iAt in range(0, data.NAtoms):
		At = data.At[iAt]
		if (At.iSpc != iSpc):
			sys.exit('data_to_molecule : ERROR - This function only works for data with atoms belonging to a single molecule')

		Atoms.insert(iAt, At)

	Bonds  = np.copy(data.Bnd)
	Angles = np.copy(data.Ang)
	
	return molecule(iSpc, Atoms, Bonds, Angles) 







