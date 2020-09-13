import numpy as np

from Hist   import Histogram  # parent class
from typing import Dict       # used to enforce the parent class typing

class WL_histogram(Histogram):
  # Initializes a Wang-Landau histogram
  # the variable wts_hist (stands for weights) should be of the Histogram type
  def __init__(self,EEHisto: Histogram,lnf,lnf_scaler,ratio_crit, lnf_crit):

    self.wts           = EEHisto.binValue           # the array of WL weights
    self.lnf           = lnf                        # WL correction factor
    self.lnf_scaler    = lnf_scaler                 # the WL correction factor is decreases according to lnf *= lnf_scaler               
    self.ratio_crit    = ratio_crit                 # flatness ratio that results in a decrease of lnf through the scaler                      
    self.NSubs         = EEHisto.NBins              # number of sub-ensembles in the simulation
    self.lnf_crit      = lnf_crit                   # when lnf <= crit_lnf the WL_histogram is considered equilibrated
    self.isItDone      = False                      # becomes True when lnf <= crit_lnf
    NSubs              = self.NSubs
    self.NVisits       = np.zeros(NSubs,np.int64)   # number os visits logged for each sub-ensemble

    # inherited histogram variables
    self.min           = EEHisto.min
    self.max           = EEHisto.max
    self.width_bin     = EEHisto.width_bin

  # self-calling the class will return the weight of the requested sub-ensemble i_sub
  def __call__(self,i_sub):
    return self.wts[i_sub]

  def reset_visits(self):
    self.NVisits[:] = 0

  def incr_visits(self,i_sub):
    self.NVisits[i_sub] += 1

  def penalize(self,i_sub):
    self.incr_visits(i_sub)
    self.wts[i_sub] -= self.lnf
    self.check_for_convergence()

  def update_penalty(self):
    self.lnf *= self.lnf_scaler
    self.reset_visits()

  # checks the flatness of the WL histogram and update lnf and isItDone accordingly
  def check_for_convergence(self):
    ratio_current = np.min(self.NVisits/np.mean(self.NVisits))
    if (ratio_current > self.ratio_crit):
      self.update_penalty()
    if (self.lnf <= self.lnf_crit):
      self.isItDone = True

  # write the WL histogram in a file
  def write(self,tag,file):
    dev_from_mean    = np.zeros(self.NSubs) # initialize the deviation from mean visits column

    file.write("# %s\n" %tag)
    file.write('# lnf  =  %.3g \n' %(self.lnf))
    file.write("#\n")
    file.write("#   subEnsCoord   devFromMeanVisits      DOS\n")
    file.write("#\n")
    mean_visits = np.mean(self.NVisits)

    if (mean_visits == 0.0):
      pass
    else:
      dev_from_mean[:] = self.NVisits[:] / mean_visits

    for i_sub in range(self.NSubs):
      subEnsCoord = self.min + i_sub*self.width_bin
      file.write('%12.3f %16.3f %16.5f \n' \
               %(subEnsCoord, dev_from_mean[i_sub], self(i_sub)))

    file.write("\n")
    file.flush()