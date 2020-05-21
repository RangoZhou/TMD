import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
# from rdkit import DataStructs
# from rdkit.Chem.rdCoordGen import AddCoords
import copy

import sys
import os

input_mol_path = sys.argv[1]
conf_num = int(sys.argv[2])
out_folder = sys.argv[3]
out_put_name = sys.argv[4]

if input_mol_path.find('.pdb') != -1:
    input_mol = Chem.rdmolfiles.MolFromPDBFile(input_mol_path,removeHs=False)
if input_mol_path.find('.mol2') != -1:
    input_mol = Chem.rdmolfiles.MolFromMol2File(input_mol_path,removeHs=False)

# if input_mol.find('1Y26') != -1:
#     input_mol.GetAtomWithIdx(0).SetNumExplicitHs(1)
# input_mol = Chem.RemoveHs(input_mol)
# input_mol = Chem.AddHs(input_mol)

initial_mol = copy.deepcopy(input_mol)
template_mol = copy.deepcopy(input_mol)
template_mol.RemoveAllConformers()

confids = AllChem.EmbedMultipleConfs(input_mol, numConfs=conf_num-1, numThreads=8)
conf_ids = [cid for cid in confids]
# initial_conf_id = input_mol.AddConformer(initial_mol.GetConformer(0))
# conf_ids = [initial_conf_id] + conf_ids
generate_mols = [ copy.deepcopy(template_mol) for i in range(len(conf_ids)+1)]
print('num of confs:', len(generate_mols))
generate_mols[0].AddConformer(initial_mol.GetConformer(0))
for i in range(1,len(conf_ids)+1):
    generate_mols[i].AddConformer(input_mol.GetConformer(i-1))
print('len generate_mols: ', len(generate_mols))
# rmslist = []
# AllChem.AlignMolConformers(input_mol, RMSlist=rmslist)
# print('rmslist:', rmslist)

# res = AllChem.MMFFOptimizeMoleculeConfs(input_mol)
# print("optimized energy")
# optimized_energy_list = [r[1] for r in res]
with open('{}/{}_record.txt'.format(out_folder, out_put_name), 'w') as f:
    for i, mol in enumerate(generate_mols):
        mol_MMFF_properties = rdkit.Chem.rdForceFieldHelpers.MMFFGetMoleculeProperties(mol)
        mol_MMFF_forcefield = rdkit.Chem.rdForceFieldHelpers.MMFFGetMoleculeForceField(mol,mol_MMFF_properties)
        energy = rdkit.ForceField.rdForceField.ForceField.CalcEnergy(mol_MMFF_forcefield)
        save_path = '{}/{}_{}.pdb'.format(out_folder, out_put_name, i)
        f.write('{}_{} {:.2f}\n'.format(out_put_name, i, energy))
        pdbwriter = Chem.rdmolfiles.PDBWriter(save_path)
        pdbwriter.write(mol)