#import subprocess
#import sys
import os
#import Midas
from chimera import runCommand
from WriteMol2 import writeMol2

input_file = ''
output_file = ''
with open('/tmp/prepare_tmd.tmp') as f:
    items = f.readlines()[0].split()
    input_file = items[0]
    output_file = items[1]
    add_charge = items[2]
# print(input_file, output_file)
# os.system('pwd')
runCommand("open "+input_file)
if add_charge == "True":
    runCommand("delete element.H")
    runCommand("addh")
    runCommand("addcharge all method gas")
s="write format mol2 0 {}".format(output_file)
runCommand(s)
runCommand("close 0")