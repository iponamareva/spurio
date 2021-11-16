import os
import sys
import subprocess
import numpy as np
from constants import exp_name

import num_sequences

NUM = getattr(num_sequences, exp_name + "_NUM")
print("NUM sequences retrieved:", NUM)

stepsize = 1000
for v in [1]:#,3,5,10]:
    for i in np.arange(0, NUM, stepsize):
        log_file_name = "logs/" + exp_name + "/blastp_search_" + str(i) + "_" + str(i + stepsize) + ".log" 
        os.system('bsub -M 6000 -n 2 -o {} \'python blastp_search.py -s {} -e {}\''.format(log_file_name, int(i), int(i+stepsize)))
