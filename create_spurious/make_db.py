import os
from constants import *

command = "makeblastdb -in {} -out db/{}/database -parse_seqids -dbtype prot".format(proteome_path, exp_name)
print("Executing:", command)
os.system(command)
