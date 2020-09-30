import numpy as np

from Hist   import Histogram  # parent class
from typing import Dict       # used to enforce the parent class typing

class WL_histogram(Histogram):
  wts           = 0.0                        # the array of WL weights
  lnf           = 0.0                        # WL correction factor
  lnf_scaler    = 1.0                        # the WL correction factor is decreases according to lnf *= lnf_scaler               
  ratio_crit    = 0.0                        # flatness ratio that results in a decrease of lnf through the scaler                      
  NSubs         = 0                          # number of sub-ensembles in the simulation
  lnf_crit      = 0.0                        # when lnf <= crit_lnf the WL_histogram is considered equilibrated
  isItDone      = False                      # becomes True when lnf <= crit_lnf
  NVisits       = 0                          # number os visits logged for each sub-ensemble

  # 
  def __init__(self,EEHisto: Histogram,lnf,lnf_scaler,ratio_crit, lnf_crit):
    """
    Initializes a Wang-Landau histogram

    Arguments : 
    EEHisto     -  A histogram that defines the Expanded sub-ensemble
    lnf         -  Initial value of the WL correction factor
    lnf_scaler  -  The WL correction factor is decreases according to lnf *= lnf_scaler
    ratio_crit  -  Flatness ratio that results in a decrease of lnf through the scaler 
    lnf_crit    -  When lnf <= crit_lnf the WL_histogram is considered equilibrated

    The variable EEHisto should be of the Histogram type
    """

    self.NSubs         = EEHisto.NBins              
    NSubs              = self.NSubs
    self.wts           = np.zeros(NSubs,np.double)
    self.wts           = np.copy(EEHisto.binValue)         
    self.lnf           = lnf                        
    self.lnf_scaler    = lnf_scaler                               
    self.ratio_crit    = ratio_crit                                     
    self.lnf_crit      = lnf_crit                   
    self.isItDone      = False                      
    self.NVisits       = np.zeros(NSubs,np.int64)   

    # inherited histogram variables
    self.min           = EEHisto.min
    self.max           = EEHisto.max
    self.width_bin     = EEHisto.width_bin
    self.binValue      = self.wts

  #-----------------------------------------------------
  def __call__(self,i_sub):
    """
    Self-calling the class will return the weight 
    of the requested sub-ensemble i_sub
    """
    return self.wts[i_sub]

  #-----------------------------------------------------
  def reset_visits(self):
    self.NVisits[:] = 0

  #-----------------------------------------------------
  def incr_visits(self,i_sub):
    self.NVisits[i_sub] += 1

  #-----------------------------------------------------
  def penalize(self,i_sub):
    self.incr_visits(i_sub)
    self.wts[i_sub] -= self.lnf
    self.check_for_convergence()

  #-----------------------------------------------------
  def update_penalty(self):
    self.lnf *= self.lnf_scaler
    self.reset_visits()

  #-----------------------------------------------------
  def check_for_convergence(self):
    """
    Checks the flatness of the WL histogram and 
    updates lnf and isItDone accordingly
    """
    ratio_current = np.min(self.NVisits/np.mean(self.NVisits))
    if (ratio_current > self.ratio_crit):
      self.update_penalty()
    if (self.lnf <= self.lnf_crit):
      self.isItDone = True

  #-----------------------------------------------------
  def write(self,tag,file):
    """Write the WL histogram in a file"""
    dev_from_mean    = np.zeros(self.NSubs) # initialize the deviation from mean visits column

    file.write("# %s\n" %tag)
    file.write('# lnf  =  %.3g \n' %(self.lnf))
    file.write("#\n")
    file.write("#   subEnsCoord   devFromMeanVisits      DOS\n")
    file.write("#\n")

    dos         = np.zeros(self.NSubs,np.double)
    dosMax      = np.amax(self.wts[:])
    dos         = self.wts[:] - dosMax


    mean_visits = np.mean(self.NVisits)
    if (mean_visits == 0.0):
      pass
    else:
      dev_from_mean[:] = self.NVisits[:] / mean_visits

    for i_sub in range(self.NSubs):
      subEnsCoord = self.min + i_sub*self.width_bin
      file.write('%12.3f %16.3f %16.5f \n' \
               %(subEnsCoord, dev_from_mean[i_sub], dos[i_sub]))

    file.write("\n")
    file.flush()