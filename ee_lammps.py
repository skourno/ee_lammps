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

# initialize the following auxiliary variables
acceptTrans  = False
idx_test_dir = 0      # index that decides the direction

NStepsSE  = sim.NSteps_subEns
NStepsP   = sim.NSteps_prod
NLoops    = int(NStepsP / NStepsSE)
NStepsRan = 0

wl_used   = sim.use_wl_bool
tmmc_used = sim.use_tmmc_bool
NSubs     = sim.EEHist.NBins

if (comm.Get_rank() == 0):
  print("----------------------------------------------------------")
  print(" > Starting the Expanded Ensemble exploration ...")
  print("\n")
  print("    MD steps - Sub index - EE coord - PE (kcal/mol) - idx testCat - idx")
  print("\n")


for i_loop in range(NLoops):

  # update the TMMC Histogram
  if (tmmc_used and np.mod(NStepsRan,sim.NStepsUpdateTM)):
    sim.TMHist.update_TMMC_weights()

  timeStamp = "Simulation step: %8d" %(NStepsRan)
  if (sim.write_wl   and np.mod(NStepsRan,sim.NWStepWL) == 0 and comm.Get_rank() == 0):
    sim.WLHist.write(timeStamp, sim.WL_FILE_out)
  
  if (sim.write_tmmc and np.mod(NStepsRan,sim.NWStepTM) == 0 and comm.Get_rank() == 0):
    tag = "Simulation step: %8d" %(NStepsRan)
    sim.TMHist.write(timeStamp, sim.TM_FILE_out) 
  
  comm.Barrier()
  run_lammps_sim(lmp,NStepsSE)
  NStepsRan += NStepsSE

  pe_old   = lmp.extract_compute("thermo_pe",0,0)
  iSub_old = sim.EEHist.idx_of(testPart.ee_coord())


  if (testPart.Type =='Ion Pair'): 
    charge = testPart.ee_coord()
    if (charge == testPart.fullCharge): 
      # fully charged and we can suffle the test particles
      testPart.shuffle_testPart(lmp,comm)

  if (comm.Get_rank() == 0):
    print("%10d %10d %10.2f %15.1f %12d %12d" \
      %(NStepsRan+sim.NSteps_equil, iSub_old, testPart.charge, pe_old, testPart.idxCat, testPart.idxAn))


  # Processor that has rank zero decides which direction to move
  # and broadcasts to everybody else
  if (comm.Get_rank() == 0):
    rand_num = np.random.rand()
    if (rand_num < 0.5):
      idx_test_dir = +1
    else:
      idx_test_dir = -1
  idx_test_dir = comm.bcast(idx_test_dir,root=0)

  # If the attempted transition is going out of limits, reject it
  if ((iSub_old == 0       and  idx_test_dir == -1) or \
      (iSub_old == NSubs-1 and  idx_test_dir ==  1)):
    
    if (wl_used):
      sim.WLHist.penalize(iSub_old)  
    
    if (tmmc_used):
      trans_prob = 0.0 
      sim.TMHist.update_collection_matrix(iSub_old, idx_test_dir, trans_prob)
  else:
    # Make the temporary changes to the test particles
    # Note that these changes have to be reverted if the
    # attempted transition gets rejected
    testPart.subEns_change(lmp,idx_test_dir)
    iSub_new  = sim.EEHist.idx_of(testPart.ee_coord())

    # Get the energy by doing a false run
    comm.Barrier()
    lmp.command("run 0")
    pe_new = lmp.extract_compute("thermo_pe",0,0)

    # compute the trans_probability
    delta_pe    = (pe_new - pe_old)
    expDE       = np.exp(-sim.Beta * delta_pe)
    trans_prob  = np.minimum(expDE,1.0)
    
    if (tmmc_used):
      # update the collection matrix of TMMC
      sim.TMHist.update_collection_matrix(iSub_old, idx_test_dir, trans_prob)
    
    # Processor that has rank zero decides to either
    # accept or reject the transition and broadcasts
    # the decision to others

    if (comm.Get_rank() == 0):
      delta_w     = sim.WLHist(iSub_new) - sim.WLHist(iSub_old)
      arg         = np.exp(-sim.Beta * delta_pe + delta_w)
      acceptTrans = (arg > np.random.rand())
    
    acceptTrans = comm.bcast(acceptTrans,root=0)

    if (acceptTrans):
      if (wl_used):
        sim.WLHist.penalize(iSub_new)
    else:
      if (wl_used):
        sim.WLHist.penalize(iSub_old)

      # change back to the previous sub ens
      testPart.subEns_change(lmp,-idx_test_dir) 

      # Need this false run to update the forces correctly
      comm.Barrier()
      lmp.command("run 0")


if (sim.write_wl):
  sim.WL_FILE_out.close()

if (sim.write_tmmc):
  sim.TM_FILE_out.close()

