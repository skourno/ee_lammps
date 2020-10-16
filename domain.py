import numpy as np
import sys

class domain:
	"""Only rectangular simulation boxes allowed"""
	edge = np.zeros(3,np.double)

	def __init__(self, edge):
		self.edge = np.copy(edge)

	def fold(self, xyz):
		Img = np.zeros(0,np.integer) 
		R   = np.copy(xyz)

		R_gt_edge = np.array([r > e   for r,e in zip(R,self.edge)] )
		R_lt_zero = np.array([r < 0.0 for r   in R] )

		Img = np.floor(R/self.edge)
		Img = np.array([ int(i) for i in Img ])

		R   = np.where(R_gt_edge, R - self.edge *   np.floor(R/self.edge), R )
		R   = np.where(R_lt_zero, R + self.edge *  -np.floor(R/self.edge), R )

		return R, Img

	def unfold(self, xyz, Img):
		return np.array([r + e*i for r,e,i in zip(xyz, self.edge, Img)])

