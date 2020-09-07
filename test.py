# A script to perform WL expanded ensemble (biased) MD simulation
# to obtain free-energy estimates of charge perturbation for a pair
# of ions using LAMMPS.

# Written by Srikanth Ravipati (rsri131@gmail.com)

# Not used in a production run; only for debugging
def find_fractional_charge(q_atoms,type_atoms,idc,natoms):
    q_Na =  1.0
    id_atom_Na = -1
    q_Cl = -1.0
    id_atom_Cl = -1
    for i_atom in range(natoms):
       if (type_atoms[i_atom] == 3 or type_atoms[i_atom] == 4):
          if (abs(q[i_atom]) < 1.0):
             if (type_atoms[i_atom] == 3):
                q_Na = q[i_atom]
                id_atom_Na = idc[i_atom]
             else:
                q_Cl = q[i_atom]
                id_atom_Cl = idc[i_atom]
    return id_atom_Na, q_Na, id_atom_Cl, q_Cl
### To debug at any time
#####x = lmp.extract_atom("x",3)
#q = lmp.extract_atom("q",2)
#t = lmp.extract_atom("type",0)
#idx = lmp.extract_atom("id",0)
#data = find_fractional_charge(q,t,idx,natoms)



import numpy as np
from mpi4py import MPI
from Hist import Histogram
from WL import WL_histogram

comm = MPI.COMM_WORLD

# General EE info
# n_subensembles, lnf_init, acc_ration, lnf_crit
Hist_EE = Histogram(0.0,1.0,0.1,11)

# WL info
# n_subensembles, lnf_initial, acc_ratio_crit, lnf_crit
Hist_WL = WL_histogram(Hist_EE.size,1.0,0.75,1e-4)

# File to write DOS from WL
ee_hist_WL = open("ee_hist_WL.dat","w")


# Needs to be passed as arguments
iseed = 623469
n_atoms_water = 6402
i_atom_Na_i = 6403
i_atom_Na_f = 6460
i_atom_Cl_i = 6461
i_atom_Cl_f = 6518
natoms = 6518
idx_test_Na = 6460
idx_test_Cl = 6518

delta_q_test = Hist_EE.width_bin

Temp = 298.0
Beta = 1.0 / (Temp * 8.314 / 4184.0)

from lammps import lammps
lmp = lammps(cmdargs=["-screen","none"])
lmp.command("variable seed equal 582783")
lmp.command("units real")
lmp.command("atom_style full")
lmp.command("read_data tip4p05.data")
lmp.command("include tip4p05.ff")
lmp.command("variable Text equal 298.0")
lmp.command("variable Pext equal 1.0")
lmp.command("velocity all create ${Text} 1234")
lmp.command("neighbor 2.0 bin")
lmp.command("neigh_modify every 1 delay 0 check yes")

lmp.command("thermo_style custom step temp press epair evdwl ecoul elong")
lmp.command("thermo_modify flush yes")
lmp.command("thermo 100") 

lmp.command("timestep 2.0")

lmp.command("fix constrain all shake 1.0e-4 100 0 b 1 a 1")
lmp.command("fix integrate all nvt temp ${Text} ${Text} 100.0")
        #iso ${Pext} ${Pext} 1000.0")
lmp.command("fix removeMomentum all momentum 1 linear 1 1 1")

lmp.command("group na type 3")
lmp.command("group cl type 4")


# Initializing variables
lmp.command("variable idx_test_Na string -1")
lmp.command("variable idx_test_Cl string -1")
lmp.command("variable q_test_Na string Inf")
lmp.command("variable q_test_Cl string Inf")

# Fixing test ions/atom numbers
flag = lmp.set_variable("idx_test_Na",idx_test_Na)
flag = lmp.set_variable("idx_test_Cl",idx_test_Cl)


lmp.command("run 5000")
q_current = 1.0
idx_test_dir = 19999999999999
accept = False

idx_shift_Na = -10000000000
idx_shift_Cl = -10000000000

for i_loop in range(100000):
   if (np.mod(i_loop,500) == 0 and comm.Get_rank() == 0):
      Hist_WL.write_to_a_file(ee_hist_WL)
   comm.Barrier()
   lmp.command("run 20")
   e_old = lmp.extract_compute("thermo_pe",0,0)
   old_idx_hist = Hist_EE.idx_hist(q_current)

   # If all ion pairs have full fractional charges, 
   # randomly select ions to be the test pair
   if (old_idx_hist == Hist_EE.size-1 and comm.Get_rank() == 0):
      idx_shift_Na = np.random.randint(0,58)
      idx_test_Na  = i_atom_Na_i + idx_shift_Na 
      idx_shift_Cl = np.random.randint(0,58)
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
         idx_test_dir = 1
      else:
         idx_test_dir = -1
   idx_test_dir = comm.bcast(idx_test_dir,root=0)

   # If the attempted transition is going out of limits,
   # reject
   if ((old_idx_hist == 0 and idx_test_dir == -1) or \
       (old_idx_hist == Hist_EE.size-1 and  idx_test_dir == 1)):
       Hist_WL.penalize(old_idx_hist)  
   else:
       # Make the temporary changes to the charge values
       # Note that these changes have to be reverted if the
       # attempted transition gets rejected
       delta_q   = idx_test_dir * delta_q_test
       q_test_Na = q_current + delta_q
       q_test_Cl = -q_current - delta_q
       new_idx_hist = Hist_EE.idx_hist(q_test_Na)

       comm.Barrier()
       flag = lmp.set_variable("q_test_Na", q_test_Na)
       flag = lmp.set_variable("q_test_Cl", q_test_Cl)
       lmp.command("set atom ${idx_test_Na} charge ${q_test_Na}")
       lmp.command("set atom ${idx_test_Cl} charge ${q_test_Cl}")

       # Get the energy by doing a false run
       lmp.command("run 0")
       e_new = lmp.extract_compute("thermo_pe",0,0)

       # Processor that has rank zero decides to either
       # accept or reject the transition and broadcasts
       # the decision to others
       if (comm.Get_rank() == 0):
          delta_e = (e_new - e_old)
          arg = np.exp(-Beta * delta_e + Hist_WL.wts[new_idx_hist] - \
                Hist_WL.wts[old_idx_hist])
          accept = (arg > np.random.rand())
       accept = comm.bcast(accept,root=0)

       if (accept):
          Hist_WL.penalize(new_idx_hist)
          q_current = q_test_Na
       else:
          Hist_WL.penalize(old_idx_hist)
          q_test_Na = q_current
          q_test_Cl = -q_current
          comm.Barrier()
          flag = lmp.set_variable("q_test_Na", q_test_Na)
          flag = lmp.set_variable("q_test_Cl", q_test_Cl)
          lmp.command("set atom ${idx_test_Na} charge ${q_test_Na}")
          lmp.command("set atom ${idx_test_Cl} charge ${q_test_Cl}")
          # Need this false run to update the forces correctly
          lmp.command("run 0")

ee_hist_WL.close()
