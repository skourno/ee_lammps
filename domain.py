import numpy as np

class domain:
	"""Only rectangular simulation boxes allowed"""
	edge = np.zeros(3,np.double)

	def __init__(self, edge):
		self.edge = np.copy(edge)

	def fold(self, xyz):
		R = np.copy(xyz)
		R = np.where(R > self.edge, R - self.edge *   np.floor(R/self.edge)     , R ) 
		R = np.where(R < 0.0      , R + self.edge * (-np.floor(R/self.edge) + 1), R )
		return R
