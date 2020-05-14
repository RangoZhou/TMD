import numpy as np
#import math
import sys
import os
import os.path
from shutil import copy
#import subprocess

def get_sybyl_dict(path):
    sybyl_dict = {}
    flag = False
    with open(path) as f:
        for line in f:
            if line.strip() == "":
                continue
            if line[0:9] == "@<TRIPOS>":
                flag = False
            if line[0:13] == "@<TRIPOS>ATOM":
                flag = True
                continue
            if flag == True:
                mol2_items = line.strip().split()
                atom_name = mol2_items[1]
                # res_id = mol2_items[6]
                res_name = mol2_items[7]
                sybyl_type = mol2_items[5]
                identifier = res_name + '-' + atom_name
                if identifier in sybyl_dict:
                    print(identifier)
                    exit()
                sybyl_dict[identifier] = sybyl_type
    return sybyl_dict

in_dat = sys.argv[1]
assert os.path.isfile(in_dat), 'Error: no '+in_dat+' file found!'

f = open(in_dat)
cases = [line.strip() for line in f if line[0] != "#" and line.strip() != ""]
f.close()
print(cases)

for case in cases:
    mol2_path = "./{}/{}_lig_wH.mol2".format(case,case)
    pdbqt_path = "./{}/{}_lig_wH.pdbqt".format(case,case)
    save_pdbqt_path = "./{}/{}_lig_wH.tmd".format(case,case)
    assert os.path.isfile(mol2_path), 'Error: no '+mol2_path+' file found!'
    assert os.path.isfile(pdbqt_path), 'Error: no '+pdbqt_path+' file found!'

    #sybyl_dict atom_name-res_name : sybyl_type
    sybyl_dict = get_sybyl_dict(mol2_path)
    # print(sybyl_dict)
    # exit()
    context = []
    with open(pdbqt_path) as f:
        for line in f:
            line = line.strip()
            if line[0:4]=="ATOM" or line[0:6]=="HETATM":
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                element = line[77:79].strip()
                if element!="H" and element!="HD":
                    identifier = res_name + '-' + atom_name
                    sybyl_type = sybyl_dict[identifier]
                    sline = '{:<82}{:<5}\n'.format(line,sybyl_type)
                    context.append(sline)
                else:
                    sline = '{:<82}{:<5}\n'.format(line,"H")
                    context.append(sline)
            else:
                context.append(line+'\n')

    with open(save_pdbqt_path,'w') as f:
        for line in context:
            f.write(line)

# ligand_pdb_path = predock_folder+pdb+"_ligand_"+l+".pdb"
# #receptor_pdb_path = predock_folder+pdb+"_receptor.pdb"
# if docking_folder.find("rigid_ligand")!=-1:
#     in_ligand_path = save_folder+pdb+"_rigid_ligand_"+l+".pdbqt"
#     save_ligand_path = save_folder+pdb+"_rigid_ligand_"+l+"_vina.pdbqt"
#     docking_log_path = save_folder+pdb+"_rigid_ligand_"+l+"_vina.log"
# else:
#     in_ligand_path = save_folder+pdb+"_ligand_"+l+".pdbqt"
#     save_ligand_path = save_folder+pdb+"_ligand_"+l+"_vina.pdbqt"
#     docking_log_path = save_folder+pdb+"_ligand_"+l+"_vina.log"
# docking_conf_path = save_folder+pdb+"_"+l+"_vina.conf"
# in_receptor_path = save_folder+pdb+"_receptor.pdbqt"