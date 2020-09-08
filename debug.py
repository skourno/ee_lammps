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
