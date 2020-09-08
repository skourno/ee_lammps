import numpy as np
import sys
# -------------------------------------------------------------
# A class to support TMMC simulation runs. A transition matrix
# is stored and can be used to estimate the free energy landscale
# implied by the subensembles. 
# -------------------------------------------------------------
class TMMC_histogram:


      # Initialize just with the number of sub-ensembles in the simulation
      def __init__(self,n_sub):
          self.n_sub     = n_sub
          self.visits    = np.zeros(n_sub,np.int64)
          self.wts       = np.zeros(n_sub,np.double)
          self.CM        = np.zeros((n_sub  ,3),np.double) # collection matrix
          self.TM        = np.ones( (n_sub-1,2),np.double) # transition matrix


      # the transition probability between two neighboring sub-ensembles 
      # from i_sub to a neighbor is stored in the collection matrix. Instead
      # of giving as input the second sub-ensemble, a direction iDir is
      # given instead. Applicable values are iDir=-1,0,+1. -1 means a decrease
      # of the sub-ensemble coordinate, +1 an increase and 0 a move that
      # keeps us in sub-ensemble i_sub
      def update_collection_matrix(self,i_sub,iDir,trans_prob):
          if (iDir < -1 or iDir > 1):
            sys.exit('update_collection_matrix : ERROR - Invalid direction')
            
          i_sub_local = 1
          j_sub_local = i_sub_local + iDir
          
          self.CM[i_sub,j_sub_local] +=  trans_prob 
          self.CM[i_sub,i_sub_local] += -trans_prob + 1.0 




      # use this to compute a weight/free energy and transition matrix
      # estimate based on the current collection matrix. 
      def update_TMMC_weights(self):

        self.wts    = np.zeros(n_sub,np.double) # re-initialize
        self.wts[0] = 0.0                       # set i_sub=0 as the reference state


        # loop over all entries of the transition matrix
        for i_sub in range(self.n_sub-1):
          self.TM[i_sub,0]  =  self.CM[i_sub  ,2] / np.sum(self.CM[i_sub  ,:]) # trans i_sub -> i_sub+1
          self.TM[i_sub,1]  =  self.CM[i_sub+1,0] / np.sum(self.CM[i_sub+1,:]) # trans i_sub <- i_sub+1

          # compute the weight of i_sub+1
          self.wts[i_sub+1] =  self.wts[i_sub] + numpy.log( self.TM[i_sub,0] / self.TM[i_sub,1] )




      # ----------------------------------
      def incr_visits(self,i_sub):
          self.visits[i_sub] += 1



      # output a weight estimate in file_out
      def write_to_a_file(self,file_out):
          file_out.write("\n")
          for i_sub in range(self.n_sub):
             file_out.write('%6.5f %6.5f %6.5f \n' %(i_sub*0.1, self.visits[i_sub]/mean_visits, \
                       self.wts[i_sub]))
          file_out.write("\n")
          file_out.flush()
