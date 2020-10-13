import numpy as np
import sys

from operator import add
from domain   import domain

class atom:
	"""An atom class. Iniitialize with global index, index of Species, type, charge and xyz coords vec"""
	idx    = -1
	iSpc   = -1
	Type   = -1
	xyz    = []
	charge = 0.0

	def __init__(self, idx, iSpc, Type, charge, xyz):
		self.idx    = idx
		self.iSpc   = iSpc
		self.Type   = Type
		self.charge = charge
		self.xyz    = np.copy(xyz)


	def change_coords(self, xyz_new, box: domain):
		self.xyz   = np.copy(xyz_new)
		self.xyz   = box.fold(self.xyz)


	def displace(self, displace, box: domain):
		self.xyz = list( map(add, self.xyz, displace) )
		self.xyz = box.fold(self.xyz)
