import sys
# -------------------------------------------------------------
# A class to read the input file of the program and store the
# relevant data.
# -------------------------------------------------------------
class input_data:
	def __init__(self,InputFile):

		NLinesRead = 0

		with open(InputFile) as input:
			for line in input:
				if (line == "\n"):
					NLinesRead += 1
					continue # empty line

				# ignore anything that follows after a '#' as a comment 
				line        = line.partition('#')[0]
				line        = line.rstrip()
				lineArgs    = line.split()
				NLinesRead += 1

				if  (lineArgs[0] == 'ee_histo' ):
					error = self.ee_histo(lineArgs)
				elif(lineArgs[0] == 'use_wl'   ):
					error = self.use_wl(lineArgs)
				elif(lineArgs[0] == 'use_tmmc' ):
					pass
				elif(lineArgs[0] == 'iseed'    ):
					self.iseed     = lineArgs[1]
				elif(lineArgs[0] == 'read_data'):
					self.DataFile  = lineArgs[1]
				elif(lineArgs[0] == 'sim_time' ):
					pass
				elif(lineArgs[0] == 'set_ionPair_tp'):
					pass
				elif(lineArgs[0] == 'set_temp'):
					pass
				elif(lineArgs[0] == 'init_dos'):
					pass
				elif(lineArgs[0] == 'write'   ):
					pass
				else:
					print('input_data.__init__ : ERROR - Unrecognised input file entry\n')
					print('                      %s is not recognised\n' %lineArgs[0]    )
					sys.exit()     

				if (error):
					sys.exit('input_data.__init__ : ERROR - Erroneous entry at line %d\n' %NLinesRead)

  
  # ----------------------------------
	def ee_histo(self,lineArgs):
		if (len(lineArgs) != 4):
			return True
		self.boundary_min = lineArgs[1]
		self.boundary_max = lineArgs[2]
		self.desired_inc  = lineArgs[3]

		return False 

	# ----------------------------------
	def use_wl(self,lineArgs):
		if (len(lineArgs) < 2):
			return True

		if  (lineArgs[1] == 'yes'):
			if (len(lineArgs) != 6):
				return True

			self.use_wl      = True
			self.lnf         = lineArgs[2]
			self.lnf_scaler  = lineArgs[3]  
			self.ratio_crit  = lineArgs[4]  
			self.lnf_crit    = lineArgs[5] 
		elif(lineArgs[1] == 'no' ):
			self.use_wl = False



                


				
