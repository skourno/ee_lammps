import numpy as np
import sys

from Hist   import Histogram  # parent class
from typing import Dict       # used to enforce the parent class typing


# -------------------------------------------------------------
class TMMC_histogram(Histogram):
  """
  A class to support TMMC simulation runs. A transition matrix
  is stored and can be used to estimate the free energy landscale
  implied by the subensembles. 
  """
  wts       =  0.0
  NSubs     =  0
  NVisits   =  0 
  CM        =  0.0      # collection matrix
  TM        =  0.0      # transition matrix
  activated =  False    # if True the simulation uses the wts of TMMC to roam the sub-ensembles
  allSubsLogged = False # if True it is safe to reduce the CM to weights

  # ----------------------------------
  def __init__(self,EEHisto: Histogram):
    """Initialize the class with a histogram of size NSubs"""
    NSubs               = EEHisto.NBins
    self.wts            = np.zeros(NSubs,np.double)
    self.wts            = np.copy(EEHisto.binValue)
    self.NSubs          = NSubs
    self.NVisits        = np.zeros(NSubs,np.int64)
    self.CM             = np.zeros((NSubs  ,3),np.double) 
    self.TM             = np.ones( (NSubs  ,2),np.double) 
    self.activated      = False                           
    self.allSubsLogged  = False                           

    # inherited histogram variables
    self.min            = EEHisto.min
    self.max            = EEHisto.max
    self.width_bin      = EEHisto.width_bin
    self.binValue       = self.wts

  # ----------------------------------
  def __call__(self,i_sub):
    """
    Self-calling the class will return the 
    weight of the requested sub-ensemble i_sub
    """
    return self.wts[i_sub]

  # ----------------------------------
  def check_if_all_subensembles_have_been_logged(self):
    """
    Returns True if trans probabalities for all the 
    sub-ensembles have been logged
    """
    allSubsLogged = True # assume this is true and check if false
    for i_sub in range(self.NSubs):
      if (np.sum(self.CM[i_sub,:]) == 0.0):
        # in this case no trans_prob has been logged for i_sub yet
        allSubsLogged = False
        break
    return allSubsLogged


  # ---------------------------------- 
  def update_collection_matrix(self,i_sub,iDir,trans_prob):
    """
    The transition probability between two neighboring sub-ensembles 
    from i_sub to a neighbor is stored in the collection matrix. Instead
    of giving as input the second sub-ensemble, a direction iDir is
    given instead. Applicable values are iDir=-1,0,+1. -1 means a decrease
    of the sub-ensemble coordinate, +1 an increase and 0 a move that
    keeps us in sub-ensemble i_sub
    """
    if (iDir < -1 or iDir > 1):
      sys.exit('update_collection_matrix : ERROR - Invalid direction')
      
    i_sub_local = 1
    j_sub_local = i_sub_local + iDir
    
    self.CM[i_sub,j_sub_local] +=  trans_prob 
    self.CM[i_sub,i_sub_local] += -trans_prob + 1.0 

  # ----------------------------------
  def update_TMMC_weights(self):
    """
    Use this to compute a weight/free energy and transition matrix
    estimate based on the current collection matrix. 
    """
    # first avoid doing anything if 
    if (not self.allSubsLogged):
      self.allSubsLogged = self.check_if_all_subensembles_have_been_logged()
      if (not self.allSubsLogged):
        return # leaves the current estimate to be all zeros

    self.wts[0] = 0.0  # set i_sub=0 as the reference state

    # loop over all entries of the transition matrix
    for i_sub in range(self.NSubs-1):
      self.TM[i_sub,1]  =  self.CM[i_sub+1,0] / np.sum(self.CM[i_sub+1,:]) # trans i_sub <- i_sub+1
      self.TM[i_sub,0]  =  self.CM[i_sub  ,2] / np.sum(self.CM[i_sub  ,:]) # trans i_sub -> i_sub+1


      if (self.TM[i_sub,0] == 0):
        self.allSubsLogged = False # adresses a small bug appearing when not all directions havent been logged
        return

      # compute the weight of i_sub+1
      self.wts[i_sub+1] =  self.wts[i_sub] + np.log( self.TM[i_sub,1] / self.TM[i_sub,0] )
      
      #if (irank == 0):
      #  print(i_sub+1, self.TM[i_sub,1], self.TM[i_sub,0], self.wts[i_sub+1], flush=True)



  # ----------------------------------
  def incr_visits(self,i_sub):
    self.NVisits[i_sub] += 1
  
  # ----------------------------------
  def reset_visits(self):
    self.NVisits[:]      = 0

  # ----------------------------------
  def write(self,tag,file):
    """Write the TMMC histogram in a file"""
    dev_from_mean      = np.zeros(self.NSubs) # initialize the deviation from mean visits column

    if (self.activated):
      mean_visits      = np.mean(self.NVisits)
      dev_from_mean[:] = self.NVisits[:] / mean_visits # compute deviation from mean visits
      TMMCInfoString   = "The weights of TMMC are used for roaming the sub-ensembles"
    else:
      TMMCInfoString   = "The weights of TMMC are inactive"

    file.write("# %s\n" %tag           )
    file.write('# %s\n' %TMMCInfoString)
    file.write("#\n")

    if (not self.allSubsLogged):
      file.write("# Insufficient sampling - Weights cannot be reduced from the CM\n")
      file.write("# > Printing the CM instead ...\n")
      file.write("#  iSub   dir_left       stay       dir_right\n")
      file.write("#\n")
      for i_sub in range(self.NSubs):
        file.write('%6d %12.3e %12.3e %12.3e \n' \
                 %(i_sub, self.CM[i_sub,0], self.CM[i_sub,1], self.CM[i_sub,2] ))   
    else:
      file.write("#   subEnsCoord   devFromMeanVisits      DOS       CM/dir_left      CM/stay      CM/dir_right\n")
      file.write("#\n")
      for i_sub in range(self.NSubs):
        subEnsCoord = self.min + i_sub*self.width_bin
        file.write('%12.3f %16.3f %16.5f %14.3e %14.3e %14.3e\n' \
                 %(subEnsCoord, dev_from_mean[i_sub], self(i_sub), self.CM[i_sub,0], self.CM[i_sub,1], self.CM[i_sub,2]))
    
    file.write("\n")
    file.flush()
