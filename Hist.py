import sys  
import numpy as np
class Histogram:
	# Initialize the class
	def __init__(self, min_lim,max_lim,width_bin):
		self.min       = min_lim
		self.max       = max_lim

		# the width specified by the user cannot be used always. We chose
		# the closent value that is feasible based on the specified boundaries
		NBins          = round((max_lim - min_lim) / width_bin) + 1
		feasible_width = (max_lim - min_lim)       / (NBins-1)

		self.width_bin = feasible_width
		self.NBins     = NBins
		self.binValue  = np.zeros(NBins,np.double)

	# this class is callable with argument the index of a bin "iBin". 
	# Result is the value stored in the bin
	def __call__(self,iBin):
		if (iBin < 0 or iBin > self.NBins-1):
			print("\n")
			print("Histogram.__call__ : ERROR - Bin index is out of bounds")
			print("ERROR info         : This equality should hold %d <= %d <= %d \n" %(0, iBin, self.NBins-1))
			sys.exit()

		return self.binValue[iBin]


	# find the correct bin that "value" belongs to
	def idx_of(self,value):
		return round((value - self.min)/self.width_bin)

	# write the histogram in a file. Accompany the write with a tag reference at the beginning
	def write_histo(self,tag,file):
		file.write( tag )
		file.write( "\n")
		for iBin in range(self.NBins):
			file.write('%3d    %6.5f \n' %(iBin, self(iBin)))

		file.write("\n")
		file.flush()
