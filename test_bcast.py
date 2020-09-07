from mpi4py import MPI

comm = MPI.COMM_WORLD
test_int = -1
test_bool = False
if (comm.Get_rank() == 0):
   test_int = 0
   test_bool = True

test_int = comm.bcast(test_int,root=0)
test_bool = comm.bcast(test_bool,root=0)
print(comm.Get_rank(), test_int, test_bool)
