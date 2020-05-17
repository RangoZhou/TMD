#!/bin/bash
#python script_file input_ligand num_conformer(include input) output_folder prefix_for_output_file
python rdkit_generate_conformations.py ../test/script_test/1F27_ligand.mol2 10 ../test/script_test/ 1f27
