import os
import numpy as np
# from shutil import copy
import subprocess

# def remove_substructure_from_mol2(path, save_path):
#     mol2_contents = []
#     with open(path) as (mol2):
#         flag = True
#         for mol2_line in mol2:
#             if mol2_line[0:9] == '@<TRIPOS>':
#                 if mol2_line[0:21] == '@<TRIPOS>SUBSTRUCTURE':
#                     flag = False
#                 else:
#                     flag = True
#             if flag == True:
#                 mol2_contents.append(mol2_line)

#     with open(save_path, 'w') as save_write:
#         for mol2_write_line in mol2_contents:
#             save_write.write(mol2_write_line)


def generate_pdbqt(path, save_path, log_dir, molecule_type, check_hydrogen, reserve_charge):
    cmd = ''
    if molecule_type == 'protein':
        cmd = '~/.local/mgltools_x86_64Linux2_1.5.6/bin/pythonsh ~/.local/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r ' + path + ' -o ' + save_path
    if molecule_type == 'rna':
        cmd = '~/.local/mgltools_x86_64Linux2_1.5.6/bin/pythonsh ~/.local/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r ' + path + ' -o ' + save_path + ' -U nphs_lps'
    if molecule_type == 'ligand':
        cmd = '~/.local/mgltools_x86_64Linux2_1.5.6/bin/pythonsh ~/.local/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l ' + path + ' -o ' + save_path
    if check_hydrogen == True:
        cmd += ' -A checkhydrogens'
    if reserve_charge == True:
        cmd += ' -C'
    # rigid
    # if molecule_type == 'ligand':
        # cmd += ' -Z'
    error_outfile = open(os.path.join(log_dir, 'pdbqt_error.log'), 'a')
    normal_outfile = open(os.path.join(log_dir, 'pdbqt_normal.log'), 'a')
    error_outfile.write(path + ':' + molecule_type + '\n')
    normal_outfile.write(path + ':' + molecule_type + '\n')
    subprocess.run([cmd], shell=True, stdout=normal_outfile, stderr=error_outfile)
    normal_outfile.write('=====================================================================\n')
    error_outfile.write('=====================================================================\n')
    normal_outfile.close()
    error_outfile.close()


if __name__ == '__main__':
    save_ligand_dir = './ligand-pdbqt/'
    if not os.path.exists(save_ligand_dir):
        os.makedirs(save_ligand_dir)
    # convert ligand
    for i in range(1,101):
    #for i in range(1,2):
        ligand_file = './ligand-addH/{}_H.mol2'.format(i)
        save_path = '{}lig-{}.pdbqt'.format(save_ligand_dir, i)
        # remove_substructure_from_mol2(path=mol2_file, save_path=mol2_file)
        generate_pdbqt(path=ligand_file, save_path=save_path, log_dir=save_ligand_dir, molecule_type='ligand', check_hydrogen=False, reserve_charge=True)