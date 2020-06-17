import os
import numpy as np
# from shutil import copy
import subprocess
import sys

def remove_substructure_from_mol2(path, save_path):
    mol2_contents = []
    with open(path) as (mol2):
        flag = True
        for mol2_line in mol2:
            if mol2_line[0:9] == '@<TRIPOS>':
                if mol2_line[0:21] == '@<TRIPOS>SUBSTRUCTURE':
                    flag = False
                else:
                    flag = True
            if flag == True:
                mol2_contents.append(mol2_line)

    with open(save_path, 'w') as save_write:
        for mol2_write_line in mol2_contents:
            save_write.write(mol2_write_line)


def generate_pdbqt(path, save_path, molecule_type, check_hydrogen, reserve_charge):
    cmd = ''
    if molecule_type == 'protein':
        cmd = '~/.local/mgltools_x86_64Linux2_1.5.6/bin/pythonsh ~/.local/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r ' + path + ' -o ' + save_path
    if molecule_type == 'rna':
        cmd = '~/.local/mgltools_x86_64Linux2_1.5.6/bin/pythonsh ~/.local/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r ' + path + ' -o ' + save_path + ' -U nphs_lps'
    if molecule_type == 'ligand':
        cmd = '~/.local/mgltools_x86_64Linux2_1.5.6/bin/pythonsh ~/.local/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l ' + path + ' -o ' + save_path + " -U lps "
            # "-p H -p HD -p HS -p C -p A -p N -p NA -p NS -p OA -p OS -p F -p Mg -p MG -p P -p SA -p S -p Cl -p CL -p Ca -p CA -p Mn -p MN -p Fe -p FE -p Zn -p ZN -p Br -p BR -p I"
        # + " -U ''" does not work with reserve charge, if use this charge will be re calculated
    if check_hydrogen == True:
        cmd += ' -A checkhydrogens'
    if reserve_charge == True:
        cmd += " -C "
    # rigid
    # if molecule_type == 'ligand':
        # cmd += ' -Z'
    subprocess.run([cmd], shell=True)


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



def generate_tmd(pdbqt_path, mol2_path, save_path):
    with open(pdbqt_path) as f:
        pdbqt_contents = f.readlines()
    mol2_dict = read_mol2(mol2_path)
    with open(save_path, 'w') as f:
        atom_index = 0
        for line in pdbqt_contents:
            l = line.strip()
            if l[0:4] == "ATOM" or l[0:6] == "HETATM":
                mol2_type = mol2_dict["mol2_types"][atom_index]
                residue = mol2_dict["residues"][atom_index]
                resname = l[17:20].strip()
                # assert(residue==resname)
                assert(len(l)<=82)
                assert(len(mol2_type)<=5)
                f.write('{:<82}{:<5}\n'.format(l,mol2_type))
                atom_index += 1
            else:
                f.write('{}\n'.format(l))
    with open(save_path, 'a') as f:
        for bc in mol2_dict['bond_contents']:
            f.write('{}\n'.format(bc))


if __name__ == '__main__':
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    add_charge = sys.argv[3]
    charge_type = sys.argv[4]
    if len(sys.argv) != 5:
        print('arguments num not correct!')
        exit()
    mol2_output = input_file + '.tmp1.mol2'
    with open('/tmp/prepare_tmd.tmp', 'w') as f:
        if add_charge == "True":
            f.write('{} {} True {}\n'.format(os.getcwd()+'/'+input_file,os.getcwd()+'/'+mol2_output,charge_type))
        elif add_charge == "False":
            f.write('{} {} False {}\n'.format(os.getcwd()+'/'+input_file,os.getcwd()+'/'+mol2_output,charge_type))
        else:
            print('charge argument not correct!')
            exit()
    os.system('chimera --nogui ~/TMD/script/prepare_mol2.py')

    pdbqt_input = mol2_output
    pdbqt_output = pdbqt_input + '.tmp2.pdbqt'
    remove_substructure_from_mol2(path=pdbqt_input, save_path=pdbqt_input)
    generate_pdbqt(path=pdbqt_input, save_path=pdbqt_output, molecule_type='ligand', check_hydrogen=False, reserve_charge=True)

    babel_input = pdbqt_output
    babel_output = babel_input + '.tmp3.mol2'
    babel_command="obabel -ipdbqt {} -omol2 -O {}".format(babel_input,babel_output)
    os.system(babel_command)

    # final_mol2_input = babel_output
    # final_mol2_output = final_mol2_input + '.tmp4.mol2'
    # with open('/tmp/prepare_tmd.tmp', 'w') as f:
    #     if add_charge == "True":
    #         f.write('{} {} True\n'.format(os.getcwd()+'/'+final_mol2_input,os.getcwd()+'/'+final_mol2_output))
    #     elif add_charge == "False":
    #         f.write('{} {} False\n'.format(os.getcwd()+'/'+final_mol2_input,os.getcwd()+'/'+final_mol2_output))
    #     else:
    #         print('charge argument not correct!')
    #         exit()
    # os.system('chimera --nogui ~/TMD/script/prepare_mol2.py')

    generate_tmd(pdbqt_path=pdbqt_output, mol2_path=babel_output, save_path=output_file)

    #os.system('rm ./prepare_tmd.tmp')
    #os.system('rm '+mol2_output)
    #os.system('rm '+pdbqt_output)