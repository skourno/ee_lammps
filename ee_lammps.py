#-----------------------------------------------------------------------
#                              ee_lammps.py
#-----------------------------------------------------------------------
#
# A program to perform expanded ensemble simulations of electrolyte
# solutions. Lammps is used to simulate and explore the individual
# sub-ensembles, while a biased MC boltzman criterion is used to
# attempt subensemble changes.
#
#-----------------------------------------------------------------------
#
#  Written by:
#
#  version 0.0 - Srikanth Ravipati  (rsri131@gmail.com)
#  version 1.0 - Spiros Kournopulos (skourno@gmail.com)
#
#  github : https://github.com/skourno/ee_lammps.git
#
#-----------------------------------------------------------------------

import numpy as np
import sys

from mpi4py         import MPI
from Hist           import Histogram
from WL             import WL_histogram
from tmmc           import TMMC_histogram
from lammps         import lammps
from sim_lmp        import simulation, run_lammps_sim
from input          import read_input
from test_particle  import set_tpIonPair

comm = MPI.COMM_WORLD

# The name and path to the input file is given as an argument
InputFilePath=sys.argv[1]

# Read the input file
inData       = read_input(InputFilePath)

# Initialize the necessary global variables and an instance of lammps(lmp)
SimData      = simulation(inData)




WLHist.wts    = wts_init # give an initial estimate of the weights

# run an initial equilibration run
run_lammps_sim(lmp,SimData.time_equil)

# setup an ion pair as the test particles
set_tpIonPair(inData.nameCation, inData.nameAnion, SimData, lmp)

delta_q_test  = EEHist.width_bin
q_current     = 1.0
acceptTrans   = False

# initialize the following auxiliary variables
idx_test_dir = 0 # index that decides the direction
idx_shift_Na = 0 # index to help with choosing a new Na test particle
idx_shift_Cl = 0 # index to help with choosing a new Cl test particle


for i_loop in range(100000):
  if (np.mod(i_loop,200) == 0 and comm.Get_rank() == 0):
    tag = "Simulation time: %8d fs" %(time_equil + i_loop*time_sub_sim)
    WLHist.write(tag, WL_FILE_out)

    TMHist.update_TMMC_weights()     # compute the current estimate of TMMC for the weights
    TMHist.write(tag, TM_FILE_out) 
  
  comm.Barrier()
  lmp.command("run ${time_sub_sim}")

  e_old        = lmp.extract_compute("thermo_pe",0,0)
  old_idx_hist = EEHist.idx_of(q_current)

  # If all ion pairs have full fractional charges, 
  # randomly select ions to be the test pair
  if (old_idx_hist == EEHist.NBins-1 and comm.Get_rank() == 0):
    idx_shift_Na = np.random.randint(0,NIonsPairs)
    idx_testCat  = i_atom_Na_i + idx_shift_Na 
    idx_shift_Cl = np.random.randint(0,NIonsPairs)
    idx_testAn  = i_atom_Cl_i + idx_shift_Cl

  idx_testCat = comm.bcast(idx_testCat,root=0)
  idx_testAn = comm.bcast(idx_testAn,root=0)
  flag = lmp.set_variable("idx_testCat",idx_testCat)
  flag = lmp.set_variable("idx_testAn",idx_testAn)
  
  # Processor that has rank zero decides which direction to move
  # and broadcasts to everybody else
  if (comm.Get_rank() == 0):
    rand_num = np.random.rand()
    if (rand_num < 0.5):
      idx_test_dir = +1
    else:
      idx_test_dir = -1

  idx_test_dir = comm.bcast(idx_test_dir,root=0)

  # If the attempted transition is going out of limits,
  # reject
  if ((old_idx_hist == 0 and idx_test_dir == -1) or \
      (old_idx_hist == EEHist.NBins-1 and  idx_test_dir == 1)):
      WLHist.penalize(old_idx_hist)  

      # update the TMMC collection matrix
      trans_prob = 0.0 
      TMHist.update_collection_matrix(old_idx_hist, idx_test_dir, trans_prob)
  else:
      # Make the temporary changes to the charge values
      # Note that these changes have to be reverted if the
      # attempted transition gets rejected
      delta_q      = idx_test_dir * delta_q_test
      q_testCat    = q_current + delta_q
      q_testAn     = -q_current - delta_q
      new_idx_hist = EEHist.idx_of(q_testCat)

      comm.Barrier()
      flag = lmp.set_variable("q_testCat", q_testCat)
      flag = lmp.set_variable("q_testAn", q_testAn)
      lmp.command("set atom ${idx_testCat} charge ${q_testCat}")
      lmp.command("set atom ${idx_testAn} charge ${q_testAn}")

      # Get the energy by doing a false run
      lmp.command("run 0")
      e_new = lmp.extract_compute("thermo_pe",0,0)

      # compute the trans_probability
      delta_e    = (e_new - e_old)
      expDE      = np.exp(-Beta * delta_e)
      trans_prob = np.minimum(expDE,1.0)

      # update the collection matrix of TMMC
      TMHist.update_collection_matrix(old_idx_hist, idx_test_dir, trans_prob)

      # Processor that has rank zero decides to either
      # accept or reject the transition and broadcasts
      # the decision to others
      if (comm.Get_rank() == 0):
        delta_w     = WLHist(new_idx_hist) - WLHist(old_idx_hist)
        arg         = np.exp(-Beta * delta_e + delta_w)
        acceptTrans = (arg > np.random.rand())

      acceptTrans = comm.bcast(acceptTrans,root=0)

      if (acceptTrans):
        WLHist.penalize(new_idx_hist)
        q_current = q_testCat
      else:
        WLHist.penalize(old_idx_hist)
        q_testCat = q_current
        q_testAn = -q_current
        comm.Barrier()
        flag = lmp.set_variable("q_testCat", q_testCat)
        flag = lmp.set_variable("q_testAn", q_testAn)
        lmp.command("set atom ${idx_testCat} charge ${q_testCat}")
        lmp.command("set atom ${idx_testAn} charge ${q_testAn}")
        # Need this false run to update the forces correctly
        lmp.command("run 0")

WL_FILE_out.close()

