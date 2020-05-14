#!/bin/bash
# PDB=1NEM
# #PDB=1F27
# ligand=A-27-BDG-A-28-NEB-A-29-BDR-A-30-IDG
# #ligand=A-37-BTN
# #config_path=./${PDB}_${ligand}/${PDB}_${ligand}_vina.conf
# #receptor_path=./${PDB}_${ligand}/${PDB}_receptor.pdbqt
# #ligand_path=./${PDB}_${ligand}/${PDB}_ligand_${ligand}.pdbqt
# #save_docking_path=./${PDB}_${ligand}/${PDB}_ligand_${ligand}_tmd.pdbqt
# #save_log_path=./${PDB}_${ligand}/${PDB}_ligand_${ligand}_tmd.log

# #cmd="../bin/tmd --config ${config_path} --receptor ${receptor_path} --ligand ${ligand_path} --out ${save_docking_path} --log ${save_log_path}"
# #cmd="../bin/tmd --config ${config_path} --receptor ${receptor_path} --ligand ${ligand_path} --out ${save_docking_path} --log ${save_log_path} --randomize_only"
# cmd="../bin/tmd ./sybyl_pdbqt/${PDB}_receptor.mol2 ./sybyl_pdbqt/1NEM_BDR_1_lig_wH.tmd"
# echo $cmd
# $cmd

while read line
do


num_process=`pgrep tmd | wc -l`
while [ $num_process -ge 4 ]
do
echo "current num_process: $num_process"
sleep 30
num_process=`pgrep tmd | wc -l`
done


config_path=./ligandRNA/${line}/${line}_tmd.conf
receptor_path=./ligandRNA/${line}/${line}_rec_wH.mol2
ligand_path=./ligandRNA/${line}/${line}_lig_wH.tmd
save_docking_path=./ligandRNA/${line}/${line}_lig_wH_out.pdbqt
save_log_path=./ligandRNA/${line}/${line}_lig_wH_out.log

cmd="../bin/tmd --config ${config_path} --receptor ${receptor_path} --ligand ${ligand_path} --out ${save_docking_path} --log ${save_log_path}"
echo $cmd
$cmd &
#sleep 1000000000000000000000000
echo "##############################################################################"
done < ./ligandRNA/in.dat
