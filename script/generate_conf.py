import numpy as np
#import math
import sys
import os
import os.path
from shutil import copyfile
#import subprocess


def read_mol2(file):
    flag = False
    mol_dict = {'atoms':[], 'mol2_types':[], 'residues':[], 'coords':[], 'bond_contents':[]}
    with open(file) as f:
        for line in f:
            l = line.strip()
            if l != '' and l[0] != '#':
                if l[0:13] == '@<TRIPOS>BOND':
                    flag = True
                    continue
                if l[0:9] == '@<TRIPOS>':
                    flag = False
                if flag:
                    mol_dict['bond_contents'].append("<TRIPOS>BOND "+l)
    with open(file) as f:
        for line in f:
            l = line.strip()
            if l != '' and l[0] != '#':
                if l[0:13] == '@<TRIPOS>ATOM':
                    flag = True
                    continue
                if l[0:9] == '@<TRIPOS>':
                    flag = False

                if flag:
                    items = l.split()
                    atom_name = items[1]
                    x, y, z = float(items[2]), float(items[3]), float(items[4])
                    mol2_type = items[5]
                    residue = items[7]

                    mol_dict['atoms'].append(atom_name)
                    mol_dict['mol2_types'].append(mol2_type)
                    mol_dict['residues'].append(residue)
                    mol_dict['coords'].append(np.array((x, y, z)))

    return mol_dict

def read_pdbqt(file):
    mol_dict = {'atoms':[], 'pdbqt_types':[], 'residues':[], 'coords':[], 'charges':[]}
    with open(file) as f:
        for line in f:
            l = line.strip()
            if l[0:4] == "ATOM" or l[0:6] == "HETATM":
                x = float(l[30:38])
                y = float(l[38:46])
                z = float(l[46:54])
                #record str[0:6]
                #serial = int([6:11])
                atom_name = l[12:16].strip()
                #altloc  str(16
                resname = l[17:20].strip()
                if resname == 'HOH':
                    continue
                #chainid str[21]
                #resseq int[22:26]
                #insertion
                charge = float(l[70:76].strip())
                pdbqt_type = l[77:79].strip()
                mol_dict['atoms'].append(atom_name)
                mol_dict['pdbqt_types'].append(pdbqt_type)
                mol_dict['residues'].append(resname)
                mol_dict['coords'].append(np.array((x, y, z)))
                mol_dict['charges'].append(charge)
    return mol_dict


def get_site_box(in_receptor, in_ligand):
        ligand_mol2 = read_mol2(in_ligand)
        ligand_coordinates = ligand_mol2['coords']
        ligand_coordinates = np.array(ligand_coordinates)
        ligand_center = np.sum(ligand_coordinates, axis=0)/ligand_coordinates.shape[0]

        receptor_mol2 = read_mol2(in_receptor)
        receptor_coordinates = receptor_mol2['coords']
        receptor_coordinates = np.array(receptor_coordinates)
        receptor_center = np.sum(receptor_coordinates, axis=0)/receptor_coordinates.shape[0]
        receptor_upper_xyz = np.amax(receptor_coordinates,axis=0)
        receptor_lower_xyz = np.amin(receptor_coordinates,axis=0)

        ligand_upper_xyz = np.amax(ligand_coordinates,axis=0)
        ligand_lower_xyz = np.amin(ligand_coordinates,axis=0)

        total_buffer = ligand_upper_xyz - ligand_lower_xyz
        buffer_max_size = np.max(total_buffer)
        maximum_buffer = np.array([buffer_max_size,buffer_max_size,buffer_max_size])

        boxes = (receptor_upper_xyz - receptor_lower_xyz + maximum_buffer, ligand_upper_xyz - ligand_lower_xyz + maximum_buffer)

        centers = (receptor_center, ligand_center)

        return centers, boxes



in_receptor = sys.argv[1]
in_ligand = sys.argv[2]
save_conf_path = sys.argv[3]
box_mode = sys.argv[4]

#binding_center = None
#box_size = None
centers,boxes = get_site_box(in_receptor, in_ligand)
if box_mode == 'global':
    binding_center, box_size = centers[0], boxes[0]
elif box_mode == 'local':
    binding_center, box_size = centers[1], boxes[1]
else:
    print('box mode arguments wrong')
    with open(save_conf_path,"w") as conf_file:
        pass


with open(save_conf_path,"w") as conf_file:
    conf_file.write("center_x = {:.2f}\n".format(binding_center[0]))
    conf_file.write("center_y = {:.2f}\n".format(binding_center[1]))
    conf_file.write("center_z = {:.2f}\n".format(binding_center[2]))
    conf_file.write("size_x = {:.2f}\n".format(box_size[0]))
    conf_file.write("size_y = {:.2f}\n".format(box_size[1]))
    conf_file.write("size_z = {:.2f}\n".format(box_size[2]))
    #conf_file.write("exhaustiveness = {}\n".format(20))
    conf_file.write("cpu = {}\n".format(1))
# print("---------------------------------------------------")
