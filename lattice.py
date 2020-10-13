import numpy as np
import sys

from operator import mul

class point:
	"""A simple point class, initialize with coordinates xyz"""
	xyz = np.zeros(3,np.double)

	def __init__(self,xyz):
		if (len(xyz) != 3):
			sys.exit('point.__init__ : ERROR - Points need to be 3-dimensional')

		self.xyz = np.copy(xyz)

	def __call__(self):
		return self.xyz

class unitCell:
	"""
	The unitCell is the building block of a lattice. To initialize
	the kind of the unit cell needs to be defined. The choices are
	cub (stands for cubic) / bcc / fcc 
	"""
	kind    = ''   # stores the unitCell kind cub/bcc/fcc
	p       = []   # stores the coordinates of the unitCell (edge is always 1.0)
	NNodes  = 0    # Nodes per unit cel 

	def __init__(self,kind):
		self.kind = kind

		if  (kind == 'cub'):
			self.NNodes = 1
			xyz         = np.zeros(3,np.double)
			self.p.append(point(xyz))

		elif(kind == 'bcc'):
			self.NNodes = 2
			xyz1        = [0.0 , 0.0, 0.0]
			xyz2        = [0.5 , 0.5, 0.5]
			self.p      = [point(xyz1), point(xyz2)]

		elif(kind == 'fcc'):
			self.NNodes = 4
			xyz1        = [0.0 , 0.0, 0.0]
			xyz2        = [0.5 , 0.5, 0.0]
			xyz3        = [0.5 , 0.0, 0.5]
			xyz4        = [0.0 , 0.5, 0.5]
			self.p      = [point(xyz1), point(xyz2), point(xyz3), point(xyz4)]
		else:
			sys.exit('unitCell.__init__ : ERROR - Unsupported unit Cell. Choose between cub/bcc/fcc')

	def __call__(self, iNode):
		"""
		Calling a unitCell requires a Node index and returns 
		the corresponding coords
		"""
		return self.p[iNode].xyz

	def __iter__(self):
		for iNode in range(0, self.NNodes):
			yield self.p[iNode].xyz


class lattice:
	"""
	A lattice is initialized using the desired kind (cub/bcc/fcc),
	a vector with the edge length of the unitary cell and the
	minimum number of nodes in the lattice.
	Note that the lattice will always be cubic in node numbers, but
	it can be trapezoidal by using different unitCell edges in each 
	direction.
	"""
	kind   = ""
	NCells    = 0
	NCells_1D = 0
	NNodes    = 0
	uCell     = 0
	edge      = 0.0
	scaleUC   = 0.0

	def __init__(self, kind, scaleUC, NNodesMin):
		self.uCell     =  unitCell(kind)
		self.NCells_1D =  int(np.ceil( (NNodesMin / self.uCell.NNodes )**(1.0/3.0) ))
		self.NCells    =  self.NCells_1D**3
		self.NNodes    =  self.uCell.NNodes * self.NCells
		self.scaleUC   =  np.copy(scaleUC)
		self.edge      =  [ rUC * (self.NCells_1D-1) for rUC in scaleUC ]  

	def __call__(self, iCx, iCy, iCz, iNode):
		xyz   = np.zeros(3,np.double)
		iCell = [iCx,iCy,iCz]

		xyz   = np.multiply(iCell,self.scaleUC) + self.uCell(iNode)

		return xyz


	def __iter__(self):
		for iCx in range(0, self.NCells_1D):
			for iCy in range(0, self.NCells_1D):
				for iCz in range(0, self.NCells_1D):
					for iNode in range(0, self.uCell.NNodes):
						yield self(iCx, iCy, iCz, iNode)


	def __next__(self):
		for iCx in range(0, self.NCells_1D):
			for iCy in range(0, self.NCells_1D):
				for iCz in range(0, self.NCells_1D):
					for iNode in range(0, self.uCell.NNodes):
						yield self(iCx, iCy, iCz, iNode)

		
