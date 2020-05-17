#import subprocess
#import sys
import os
#import Midas
from chimera import runCommand
from WriteMol2 import writeMol2


ligand_ids = [str(i) for i in range(1,101)]

for lig_id in ligand_ids:
    inFile = "./generated-confs-pdb/{}.pdb".format(lig_id)

    runCommand("open " + inFile)

    runCommand("delete element.H")
    runCommand("addh")
    runCommand("addcharge all method gas")
    s="write format pdb 0 ./ligand-addH/{}_H.pdb".format(lig_id)
    print(s)
    runCommand(s)
    s="write format mol2 0 ./ligand-addH/{}_H.mol2".format(lig_id)
    print(s)
    runCommand(s)

    runCommand("close 0")