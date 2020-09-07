grep Pos $1 > tmp.dat
n_moles_water=$2
n_atoms_water=$(($2*4))
n_atoms_water_3atom=$(($2*3))
n_bonds_water=$(($2*2))
n_angles_water=$(($2))


n_atoms_na=$3
n_atoms_cl=$4

echo " "
echo " "
echo $(($n_atoms_water_3atom+$n_atoms_na+$n_atoms_cl)) atoms
echo $(($n_bonds_water)) bonds
echo $((n_angles_water)) angles

echo " "

echo 4 atom types
echo 1 bond types
echo 1 angle types

echo " "

echo 0.0  $5 xlo  xhi
echo 0.0  $5 ylo  yhi
echo 0.0  $5 zlo  zhi

echo " "

n_atoms_na_0=$(($n_atoms_water+1))
n_atoms_na_f=$(($n_atoms_water+$n_atoms_na))

n_atoms_cl_0=$(($n_atoms_na_f+1))
n_atoms_cl_f=$(($n_atoms_cl_0+$n_atoms_cl))
head -n $n_atoms_water tmp.dat > water.dat
sed -i '' 'n; n; n; d' water.dat
awk '{print NR "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' water.dat > water_lmp.dat
awk -F" " '{gsub("1","OW  -1.1128",$3)}1' water_lmp.dat > water.dat
awk -F" " '{gsub("2","HW1 0.5564",$3)}1' water.dat > water_lmp.dat
awk -F" " '{gsub("3","HW2 0.5564",$3)}1' water_lmp.dat > water.dat

awk -F" " '{gsub("OW","2",$3)}1' water.dat > water_lmp.dat
awk -F" " '{gsub("HW1","1",$3)}1' water_lmp.dat > water.dat
awk -F" " '{gsub("HW2","1",$3)}1' water.dat > water_lmp.dat


sed -n "$n_atoms_na_0,$n_atoms_na_f p" tmp.dat > na.dat
awk -v nshift=$n_atoms_water_3atom -v nmol_shift=$n_moles_water '{print NR+nshift "\t" NR+nmol_shift "\t" 3 "\t" "1.000" "\t" $4 "\t" $5 "\t" $6}' na.dat > na_lmp.dat

sed -n "$n_atoms_cl_0,$n_atoms_cl_f p" tmp.dat > cl.dat
awk -v nshift=$(($n_atoms_water_3atom+$n_atoms_na)) -v nmol_shift=$(($n_moles_water+$n_atoms_na)) '{print NR+nshift "\t" NR+nmol_shift "\t" 4 "\t" "-1.000" "\t" $4 "\t" $5 "\t" $6}' cl.dat > cl_lmp.dat

cat water_lmp.dat > atoms.dat
cat na_lmp.dat >> atoms.dat
cat cl_lmp.dat >> atoms.dat

echo " "
echo "Atoms"
echo " "
cat atoms.dat

echo " "
echo "Bonds"
echo " "

j_bond=0
n_atoms=0
for ((i_mol=1;i_mol<=$n_moles_water;i_mol++));do
    for ((i_bond=1;i_bond<=2;i_bond++));do
        j_bond=$(($j_bond+1))
	echo $j_bond 1 $(($n_atoms+1)) $(($n_atoms+1+$i_bond)) 
    done 
    n_atoms=$(($n_atoms+3))
done


echo " "
echo "Angles"
echo " "

j_bend=0
n_atoms=0
for ((i_mol=1;i_mol<=$n_moles_water;i_mol++));do
    j_bend=$(($j_bend+1))
    echo $j_bend 1 $(($n_atoms+2))  $(($n_atoms+1)) $(($n_atoms+3)) 
    n_atoms=$(($n_atoms+3))
done

#rm na.dat cl.dat cl_lmp.dat na_lmp.dat water.dat water_lmp.dat atoms.dat 
