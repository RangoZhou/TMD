#!/bin/bash
#python script_file input_ligand num_conformer(include input) output_folder prefix_for_output_file
#!/bin/bash

for name in 1AKX 1BYJ 1F27 1FUF 1J7T1 1J7T2 1KOC 1KOD 1LVJ 1PBR 1UUD 1XPF 1Y26 2BE01 2BE02 2BEE1 2BEE2 2FCZ1 2FCZ2 2FD01 2FD02 3SUH 3SUX 4FEJ 4FEL 4FEN 4FEO 4FEP 4JF2 4KQY 4LVW1 4LVW2 4LVX1 4LVX2 4LVY1 4LVY2 4LVZ 4LW0
do
    #mkdir -p ./${name}/conformer
    echo $name
    python rdkit_generate_conformations.py ./${name}/${name}_ligand.pdb 200 ./${name}/conformer/ ${name}_conf
done

for name in 1FMN
do
    #mkdir -p ./${name}/conformer
    echo $name
    python rdkit_generate_conformations.py ./${name}/${name}_ligand.mol2 200 ./${name}/conformer/ ${name}_conf
done