import numpy    as np
import sys
import random

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

	NSpcOf      = []
	NSpc        = 0

	At          = []
	Bnd         = []
	Ang         = []

	Spc         = []

	Box         = domain(np.zeros(3,np.double))

	def __init__(self, Box: domain):
		self.Box = deepcopy(Box)
		pass

	def insert_spc(self,Spc):
		"""
		This function should be used when inserting molecules or 
		atoms in a Mol_System. An Spc can be either an atom or a 
		molecule.
		"""
		SpcIn   = deepcopy(Spc)
		iSpc    = self.NSpc

		if   type(SpcIn) == molecule:
			SpcIn.iSpc   = iSpc
			iSpc         = self.__insert_mol(SpcIn)
		elif type(SpcIn) == atom:
			SpcIn.iSpc   = iSpc
			idx          = self.__insert_atom(SpcIn)
		else:
			sys.exit("Mol_System.insert_spc : ERROR - Unrecognised species")
		
		self.Spc.append(SpcIn)
		self.NSpc += 1

	def rand_swaps(self,NSwaps,TargetSlice=slice(0,-1,1)):
		for iSwap in range(NSwaps):
			Spc1   = random.choice(self.Spc)
			Spc2   = random.choice(self.Spc[TargetSlice])

			print(iSwap, Spc1.iSpc,Spc2.iSpc)

			xyz1   =  Spc1.center()
			xyz2   =  Spc2.center()

			Spc1.change_coords(xyz2, self.Box)
			Spc2.change_coords(xyz1, self.Box)

	def __insert_atom(self, atomIn: atom):

		if (atomIn.iSpc > self.NSpcTypes):
			self.__expand_SpcArrays(atomIn.iSpc)

		if (atomIn.Type+1 > self.NAtomTypes):
			self.NAtomTypes += atomIn.Type - self.NAtomTypes + 1 

		atomIn.idx   = self.NAtoms
		self.At.append(atomIn)
		self.NAtoms += 1

		return atomIn.idx

	def __insert_bond(self, BndIn: bond):
		if (BndIn.Type+1 > self.NBondTypes):
			self.NBondTypes += BndIn.Type - self.NBondTypes + 1 

		BndIn.idx    = self.NBonds
		self.Bnd.append(BndIn)
		self.NBonds += 1
		return BndIn.idx

	def __insert_angle(self, AngIn: angle):
		if (AngIn.Type+1 > self.NAngleTypes):
			self.NAngleTypes += AngIn.Type - self.NAngleTypes + 1

		AngIn.idx    = self.NAngles
		self.Ang.append(AngIn)
		self.NAngles+= 1
		return AngIn.idx

	def __insert_mol(self, mol: molecule):
		for AtomIn in mol.At:
			AtomIn.iSpc = mol.iSpc
			old_idx     = AtomIn.idx
			new_idx     = self.__insert_atom(AtomIn)

		idx_shift = new_idx - old_idx

		for BndIn in mol.Bnd:
			BndIn.At    += idx_shift
			idx          = self.__insert_bond(BndIn)

		for AngIn in mol.Ang:
			AngIn.At    += idx_shift
			idx          = self.__insert_angle(AngIn)

		return mol.iSpc

	def __expand_SpcArrays(self,MaxTypeIdx):
		HowManyMore       = MaxTypeIdx - self.NSpcTypes + 1
		padding           = [0] * HowManyMore
		self.NSpcOf       = self.NSpcOf + padding
		self.NSpcTypes   += HowManyMore

def Create_config_on_lattice(latt: lattice, Directions, shuffle_mols=False, NSwaps=0):
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

	if np.sum(NSpecies) > latt.NNodes:
		sys.exit("Create_config_on_lattice : ERROR - the specified lattice is not big enough")

	SimBox    = domain(latt.edge)
	MolSys    = Mol_System(SimBox)
	iter_latt = iter(latt)      

	for Dir in Directions:
		Spc            = deepcopy(Dir[0])
		NSpeciesOfThis = int(Dir[1])

		for iSpc in range(0,NSpeciesOfThis):
			xyz = next(iter_latt) #+ 1.e-4 # add a small perturbation to avoid atoms on the surface of the lattice
			Spc.change_coords(xyz,SimBox)
			MolSys.insert_spc(Spc)

	MolSys.NSpeciesOf = np.copy(NSpecies)

	if (NAtomsTot != MolSys.NAtoms):
		print(NAtomsTot, MolSys.NAtoms)
		sys.exit("Create_config_on_lattice : ERROR - miscounted the number of atoms")

	if (shuffle_mols):
		MolSys.rand_swaps(NSwaps)


	return MolSys














