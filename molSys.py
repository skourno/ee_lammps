import numpy    as np
import sys

from   lattice  import lattice
from   atom     import atom
from   molecule import bond, angle, molecule
from   domain   import domain
from   copy     import deepcopy

class Mol_System:
	NAtoms      = 0
	NBonds      = 0
	NAngles     = 0
	NAtomTypes  = 0
	NBondTypes  = 0
	NAngleTypes = 0
	NSpcTypes   = 0

	NSpecies    = []

	At          = []
	Bnd         = []
	Ang         = []

	Box         = domain(np.zeros(3,np.double))

	def __init__(self, Box: domain):
		self.Box = np.copy(Box)
		pass

	def insert_spc(self,Spc):
		"""
		This function should be used when inserting molecules or 
		atoms in a Mol_System. An Spc can be either an atom or a 
		molecule.
		"""
		SpcIn = np.copy(Spc)
		if   type(SpcIn) == molecule:
			self.__insert_mol(Spc)
		elif type(SpcIn) == atom:
			self.__insert_atom(Spc)

	def __insert_atom(self, atomIn: atom):
		atomIn.xyz   = self.Box.fold(atomIn.xyz)

		if (atomIn.Type > self.NSpcTypes-1):
			self.__expand_NSpecies(atomIn.Type)

		atomIn.idx   = self.NAtoms
		self.At.append(atomIn)
		self.NAtoms += 1


	def __insert_mol(self,mol: molecule):
		for AtomIn in mol.At:
			self.__insert_atom(AtomIn)

		for BndIn in mol.Bnd:
			self.Bnd.append(BndIn)
			self.NBonds += 1

		for AngIn in mol.Ang:
			self.Ang.append(AngIn)
			self.NAngles += 1

	def __expand_NSpecies(self,MaxTypeIdx):
		HowManyMore     = MaxTypeIdx - self.NSpcTypes + 1
		padding         = [0] * HowManyMore
		self.NSpecies   = self.NSpecies + padding
		self.NSpcTypes += HowManyMore




 

def Create_config_on_lattice(latt: lattice, Directions, shuffle_mols=False, NShuffles=0):
	"""
	Create a configuration within a Mol_System object on a lattice latt. The 
	lattice size can't be of size smaller than the requested number of species 
	to be placed on the lattice.
	The Directions argument is a list of tuples. Each tuple must contain in 
	element 0 the molecule/atom that will be placed in the lattice and in 
	element 1 the number of species that the user wishes to place.
	Two OPTIONAL argument can be provided. If 'shuffle_mols' is True then when 
	all the species have been inserted, we perform random position swaps between 
	them. The default number of NShuffles is MolSys.NAtoms, but it can be changes
	through the last optional argument.
	"""

	NSpecies    = []
	NAtomsTot   = 0
	NBondsTot   = 0
	NAnglesTot  = 0

	# to begin we need to count the total number of atoms requested
	for Dir in Directions:
		Spc            = deepcopy(Dir[0])
		NSpeciesOfThis = int(Dir[1])

		if   type(Spc) == atom:
			NAtSpc  = 1
			NBonds  = 0
			NAngles = 0
		elif type(Spc) == molecule:
			NAtSpc = Spc.NAtoms
			NBonds = Spc.NBonds
			NAngles= Spc.NAngles
		else:
			sys.exit("Create_config_on_lattice : ERROR - Invalid species of type %s" %str(type(Spc)))

		NAtomsTot        += NAtSpc  * NSpeciesOfThis
		NBondsTot        += NBonds  * NSpeciesOfThis
		NAnglesTot       += NAngles * NSpeciesOfThis
		NSpecies.append(NSpeciesOfThis) 

	if np.sum(NSpecies) > latt.NCells:
		sys.exit("Create_config_on_lattice : ERROR - the specified lattice is not big enough")

	SimBox    = domain(latt.edge)
	MolSys    = Mol_System(SimBox)

	for Dir in Directions:
		Spc            = deepcopy(Dir[0])
		NSpeciesOfThis = int(Dir[1])

		for iSpc in range(0,NSpeciesOfThis):
			xyz = next(latt)
			Spc.change_coords(xyz,SimBox)
			MolSys.insert_spc(Spc)

	MolSys.NSpecies = np.copy(NSpecies)

	if (NAtomsTot != MolSys.NAtoms):
		sys.exit("Create_config_on_lattice : ERROR - miscounted the number of atoms")

	return MolSys














