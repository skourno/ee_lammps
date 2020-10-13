import numpy as np

from operator import add, mul
from atom     import atom
from domain   import domain
from copy     import deepcopy

class bond:
	idx   = -1
	Type  = -1
	At    = []

	def __init__(self, idx, Type, Atoms):
		self.idx    = idx
		self.Type   = Type
		self.At     = np.copy(Atoms)

	def __call__(self):
		return [ At[0].idx , At[1].idx ]


class angle:
	idx   = -1
	Type  = -1
	At    = []

	def __init__(self, idx, Type, Atoms):
		self.idx    = idx
		self.Type   = Type
		self.At     = np.copy(Atoms)	

	def __call__(self):
		return [ At[0].idx , At[1].idx, At[2].idx ]


class molecule:
	idx  = -1
	At   = []
	Bnd  = []
	Ang  = []

	NAtoms  = 0
	NBonds  = 0
	NAngles = 0

	def __init__(self, idx, Atoms, Bonds, Angles):
		self.idx   = idx
		self.At    = deepcopy(Atoms)
		self.Bnd   = deepcopy(Bonds)
		self.Ang   = deepcopy(Angles)

		self.NAtoms   = len(self.At)
		self.NBonds   = len(self.Bnd)
		self.NAngles  = len(self.Ang)

	def center(self):
		r_sum  = np.zeros(3,np.double)

		for atom in self.At:
			r_sum = list( map(add, r_sum, atom.xyz) )

		return [ ri / self.NAtoms for ri in r_sum ]


	def change_coords(self, xyz_new, box: domain):
		refPoint   = self.center()
		nRefPoint  = [ - ri for ri in refPoint ]
		displace   = list( map(add, xyz_new, nRefPoint) )

		for atom in self.At:
			atom.xyz = list( map(add, atom.xyz, displace) )
			atom.xyz = box.fold(atom.xyz)


	def displace(self, displace, box: domain):
		for atom in self.At: 
			atom.xyz = list( map(add, atom.xyz, displace) )
			atom.xyz = box.fold(atom.xyz)








