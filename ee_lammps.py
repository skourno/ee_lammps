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

from mpi4py   import MPI
from Hist     import Histogram
from WL       import WL_histogram
from tmmc     import TMMC_histogram
from lammps   import lammps
from init_sim import Initialize_tip4p_with_Ions_simulation

comm = MPI.COMM_WORLD

# General EE histogram info
charge_min  = 0.0
charge_max  = 1.0
desired_inc = 0.1
EEHist      = Histogram(charge_min, charge_max, desired_inc)

# WL info
lnf_initial = 1.e-3    # initial lnf
lnf_scaler  = 0.5    
ratio_crit  = 0.8    
lnf_crit    = 1.e-5  
WLHist      = WL_histogram(EEHist,lnf_initial,lnf_scaler,ratio_crit,lnf_crit)

# TMMC Initialization
TMHist      = TMMC_histogram(EEHist)

WL_FILE_out = open("dos_WL.dat","w") # File to write DOS from WL
TM_FILE_out = open("dos_TM.dat","w") # File to write DOS from TMMC

# Needs to be passed as arguments through an input file (TODO)
iseed         = 623469
n_atoms_water = 6402
i_atom_Na_i   = 6403
i_atom_Na_f   = 6460
i_atom_Cl_i   = 6461
i_atom_Cl_f   = 6518
natoms        = 6518
NIonsPairs    = 58
time_equil    = 5000 # length of equilibration run in fs (no sub-ens change attempts during this time)
time_sub_sim  = 20   # each sube-ens is simulated for this many fs between sub-ens changee attempts
Temp          = 298.0
DataFileName  = "tip4p05.data"


wts_init      = np.zeros(WLHist.NSubs,np.double)
wts_init[ 0]  =  0.000
wts_init[ 1]  = -1.969
wts_init[ 2]  = -10.781
wts_init[ 3]  = -23.906
wts_init[ 4]  = -43.688
wts_init[ 5]  = -70.219
wts_init[ 6]  = -105.344
wts_init[ 7]  = -146.063
wts_init[ 8]  = -194.531
wts_init[ 9]  = -248.875
wts_init[10]  = -311.344

WLHist.wts    = wts_init # give an initial estimate of the weights


lmp = Initialize_tip4p_with_Ions_simulation(DataFileName, Temp, iseed)

# run an initial equilibration run
lmp.command("variable time_equil string -1")
flag = lmp.set_variable("time_equil",time_equil)
lmp.command("run ${time_equil}")

# make auxiliary ion groups
lmp.command("group na type 3")
lmp.command("group cl type 4")

# Initializing auxiliary variables
lmp.command("variable idx_test_Na  string -1")
lmp.command("variable idx_test_Cl  string -1")
lmp.command("variable q_test_Na    string  0")
lmp.command("variable q_test_Cl    string  0")
lmp.command("variable time_sub_sim string -1")

# Fixing test ions/atom numbers
idx_test_Na   = 6460
idx_test_Cl   = 6518
flag          = lmp.set_variable("idx_test_Na" ,idx_test_Na)
flag          = lmp.set_variable("idx_test_Cl" ,idx_test_Cl)
flag          = lmp.set_variable("time_sub_sim",time_sub_sim)


Beta          = 1.0 / (Temp * 8.314 / 4184.0)
delta_q_test  = EEHist.width_bin
q_current     = 1.0
acceptTrans   = False

# initialize the following auxiliary variables
idx_test_dir = 0 # index that decides the direction
idx_shift_Na = 0 # index to help with choosing a new Na test particle
idx_shift_Cl = 0 # index to help with choosing a new Cl test particle


for i_loop in range(100000):
  if (np.mod(i_loop,10) == 0 and comm.Get_rank() == 0):
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
    idx_test_Na  = i_atom_Na_i + idx_shift_Na 
    idx_shift_Cl = np.random.randint(0,NIonsPairs)
    idx_test_Cl  = i_atom_Cl_i + idx_shift_Cl

  idx_test_Na = comm.bcast(idx_test_Na,root=0)
  idx_test_Cl = comm.bcast(idx_test_Cl,root=0)
  flag = lmp.set_variable("idx_test_Na",idx_test_Na)
  flag = lmp.set_variable("idx_test_Cl",idx_test_Cl)
  
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
      q_test_Na    = q_current + delta_q
      q_test_Cl    = -q_current - delta_q
      new_idx_hist = EEHist.idx_of(q_test_Na)

      comm.Barrier()
      flag = lmp.set_variable("q_test_Na", q_test_Na)
      flag = lmp.set_variable("q_test_Cl", q_test_Cl)
      lmp.command("set atom ${idx_test_Na} charge ${q_test_Na}")
      lmp.command("set atom ${idx_test_Cl} charge ${q_test_Cl}")

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
        q_current = q_test_Na
      else:
        WLHist.penalize(old_idx_hist)
        q_test_Na = q_current
        q_test_Cl = -q_current
        comm.Barrier()
        flag = lmp.set_variable("q_test_Na", q_test_Na)
        flag = lmp.set_variable("q_test_Cl", q_test_Cl)
        lmp.command("set atom ${idx_test_Na} charge ${q_test_Na}")
        lmp.command("set atom ${idx_test_Cl} charge ${q_test_Cl}")
        # Need this false run to update the forces correctly
        lmp.command("run 0")

WL_FILE_out.close()

