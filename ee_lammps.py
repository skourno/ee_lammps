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
from sim_lmp        import run_lammps_sim, setup_tip4p_with_Ions, set_lammps_dump
from input          import input_data
from test_particle  import test_IonPair
from simData        import simData


InputFilePath = sys.argv[1]               # The name and path to the input file is given as an argument
inData        = input_data(InputFilePath) # Read the input file
sim           = simData(inData)           # Pass the input to a system info host

comm          = MPI.COMM_WORLD            # initialize mpi
lmp           = sim.init_lammps_Sim()     # Initialize the simulation

# Import the configuration as specified in the simData
sim.import_config(lmp)            

# setup a tip4p with ions simulation        
setup_tip4p_with_Ions(lmp, sim.Temp)

comm.Barrier()
run_lammps_sim(lmp,sim.NSteps_equil)    # run an initial equilibration run

if (inData.write_dump_parsed):
  # set the dump file print info
  set_lammps_dump(lmp, inData.wstep_dump, inData.outFile_dump)

# setup an ion pair as the test particles
testPart      = test_IonPair(inData.cationName, inData.iTypeTestCat,\
                             inData.anionName,  inData.iTypeTestAn, \
                             sim, lmp) 
#accTrans      = False
#
# initialize the following auxiliary variables
#idx_test_dir = 0 # index that decides the direction
#idx_shift_Na = 0 # index to help with choosing a new Na test particle
#idx_shift_Cl = 0 # index to help with choosing a new Cl test particle
#
#
NStepsSE  = sim.NSteps_subEns
NStepsP   = sim.NSteps_prod
NLoops    = int(NStepsP / NStepsSE)
NStepsRan = 0 + sim.NSteps_equil

for i_loop in range(NLoops):
  # update the TMMC Histogram
  if (sim.use_tmmc_bool and np.mod(NStepsRan,sim.NStepsUpdateTM)):
    sim.TMHist.update_TMMC_weights()

  timeStamp = "Simulation step: %8d" %(NStepsRan)
  if (sim.write_wl   and np.mod(NStepsRan,sim.NWStepWL) == 0 and comm.Get_rank() == 0):
    sim.WLHist.write(timeStamp, sim.WL_FILE_out)
  
  if (sim.write_tmmc and np.mod(NStepsRan,sim.NWStepTM) == 0 and comm.Get_rank() == 0):
    tag = "Simulation step: %8d" %(NStepsRan)
    sim.TMHist.write(timeStamp, sim.TM_FILE_out) 
  
  comm.Barrier()
  run_lammps_sim(lmp,NStepsSE)

  pe_old   = lmp.extract_compute("thermo_pe",0,0)
  iSub_old = 0

  sys.exit("executed")
#  old_idx_hist = EEHist.idx_of(q_current)
#
#  # If all ion pairs have full fractional charges, 
#  # randomly select ions to be the test pair
#  if (old_idx_hist == EEHist.NBins-1 and comm.Get_rank() == 0):
#    idx_shift_Na = np.random.randint(0,NIonsPairs)
#    idx_testCat  = i_atom_Na_i + idx_shift_Na 
#    idx_shift_Cl = np.random.randint(0,NIonsPairs)
#    idx_testAn  = i_atom_Cl_i + idx_shift_Cl
#
#  idx_testCat = comm.bcast(idx_testCat,root=0)
#  idx_testAn = comm.bcast(idx_testAn,root=0)
#  flag = lmp.set_variable("idx_testCat",idx_testCat)
#  flag = lmp.set_variable("idx_testAn",idx_testAn)
#  
#  # Processor that has rank zero decides which direction to move
#  # and broadcasts to everybody else
#  if (comm.Get_rank() == 0):
#    rand_num = np.random.rand()
#    if (rand_num < 0.5):
#      idx_test_dir = +1
#    else:
#      idx_test_dir = -1
#
#  idx_test_dir = comm.bcast(idx_test_dir,root=0)
#
#  # If the attempted transition is going out of limits,
#  # reject
#  if ((old_idx_hist == 0 and idx_test_dir == -1) or \
#      (old_idx_hist == EEHist.NBins-1 and  idx_test_dir == 1)):
#      WLHist.penalize(old_idx_hist)  
#
#      # update the TMMC collection matrix
#      trans_prob = 0.0 
#      TMHist.update_collection_matrix(old_idx_hist, idx_test_dir, trans_prob)
#  else:
#      # Make the temporary changes to the charge values
#      # Note that these changes have to be reverted if the
#      # attempted transition gets rejected
#      delta_q      = idx_test_dir * delta_q_test
#      q_testCat    = q_current + delta_q
#      q_testAn     = -q_current - delta_q
#      new_idx_hist = EEHist.idx_of(q_testCat)
#
#      comm.Barrier()
#      flag = lmp.set_variable("q_testCat", q_testCat)
#      flag = lmp.set_variable("q_testAn", q_testAn)
#      lmp.command("set atom ${idx_testCat} charge ${q_testCat}")
#      lmp.command("set atom ${idx_testAn} charge ${q_testAn}")
#
#      # Get the energy by doing a false run
#      lmp.command("run 0")
#      e_new = lmp.extract_compute("thermo_pe",0,0)
#
#      # compute the trans_probability
#      delta_e    = (e_new - e_old)
#      expDE      = np.exp(-Beta * delta_e)
#      trans_prob = np.minimum(expDE,1.0)
#
#      # update the collection matrix of TMMC
#      TMHist.update_collection_matrix(old_idx_hist, idx_test_dir, trans_prob)
#
#      # Processor that has rank zero decides to either
#      # accept or reject the transition and broadcasts
#      # the decision to others
#      if (comm.Get_rank() == 0):
#        delta_w     = WLHist(new_idx_hist) - WLHist(old_idx_hist)
#        arg         = np.exp(-Beta * delta_e + delta_w)
#        acceptTrans = (arg > np.random.rand())
#
#      acceptTrans = comm.bcast(acceptTrans,root=0)
#
#      if (acceptTrans):
#        WLHist.penalize(new_idx_hist)
#        q_current = q_testCat
#      else:
#        WLHist.penalize(old_idx_hist)
#        q_testCat = q_current
#        q_testAn = -q_current
#        comm.Barrier()
#        flag = lmp.set_variable("q_testCat", q_testCat)
#        flag = lmp.set_variable("q_testAn", q_testAn)
#        lmp.command("set atom ${idx_testCat} charge ${q_testCat}")
#        lmp.command("set atom ${idx_testAn} charge ${q_testAn}")
#        # Need this false run to update the forces correctly
#        lmp.command("run 0")

#WL_FILE_out.close()

