import numpy as np
#import math
import sys
import os
import os.path
from shutil import copy
#import subprocess

print("usage: python generate_vina_conf.py in_dat")
in_dat = sys.argv[1]
assert os.path.isfile(in_dat), 'Error: no '+in_dat+' file found!'

with open(in_dat) as f:
    for line in f:
        if line[0] == "#" or line.strip() == "":
            continue
        name = line.strip()
        docking_folder = "./"+name+"/"
        save_folder = docking_folder
        assert os.path.isdir(docking_folder), 'Error: no '+docking_folder+' directory found!'
        assert os.path.isdir(save_folder), 'Error: no '+save_folder+' directory found!'
        
        in_ligand_path = docking_folder+name+"_lig_wH.tmd"
        docking_conf_path = docking_folder+name+"_tmd.conf"
        in_receptor_path = docking_folder+name+"_rec_wH.mol2"
        
        ligand_coordinates = []
        with open(in_ligand_path) as in_ligand:
            for lig_line in in_ligand:
                if lig_line[0:6]=="HETATM" or lig_line[0:6]=="ATOM  ":
                    ligand_coordinates.append((float(lig_line[30:38]),float(lig_line[38:46]),float(lig_line[46:54])))
        
        
        ligand_coordinates = np.array(ligand_coordinates)
        #print(ligand_coordinates,ligand_coordinates.shape)
        center = np.sum(ligand_coordinates, axis=0)/ligand_coordinates.shape[0]
        print(name)
        print('center',center)
        #print('min xyz',np.amin(ligand_coordinates,axis=0))
        #print('max xyz',np.amax(ligand_coordinates,axis=0))
        #upper_buffer = np.amax(ligand_coordinates,axis=0) - center
        #lower_buffer = center - np.amin(ligand_coordinates,axis=0)
        #print('lower_buffer',lower_buffer)
        #print('upper_buffer',upper_buffer)
        #print("total_buffer",upper_buffer+lower_buffer)
        #maximum_buffer = np.maximum(lower_buffer,upper_buffer)
        #print("maximum of lower and upper buffer",maximum_buffer)
        #size_of_box = maximum_buffer*2+10
        #print("size of the box",size_of_box)
        #max_dimension_of_box = np.amax(size_of_box,axis=0)
        #if max_dimension_of_box < 10:
        #    max_dimension_of_box = 10
        #print("max dimension of the box",max_dimension_of_box)
        
        max_dimension_of_box = 8.0
        
        #size_list.append(max_dimension_of_box)
        with open(docking_conf_path,"w") as conf_file:
            conf_file.write("center_x = {:.2f}\n".format(center[0]))
            conf_file.write("center_y = {:.2f}\n".format(center[1]))
            conf_file.write("center_z = {:.2f}\n".format(center[2]))
            conf_file.write("size_x = {:.2f}\n".format(max_dimension_of_box))
            conf_file.write("size_y = {:.2f}\n".format(max_dimension_of_box))
            conf_file.write("size_z = {:.2f}\n".format(max_dimension_of_box))
            conf_file.write("cpu = {}\n".format(16))
        print("---------------------------------------------------")
