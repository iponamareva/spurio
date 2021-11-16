import os
import sys
import subprocess
import numpy as np
from constants import exp_name

import num_sequences

NUM = getattr(num_sequences, exp_name + "_NUM")
print("NUM sequences retrieved:", NUM)

if not os.path.exists("logs/" + exp_name):
    os.makedirs("logs/" + exp_name)

stepsize = 1000
for v in [1]:#,3,5,10]:
    for i in np.arange(0, NUM, stepsize):
        log_file_name = "logs/" + exp_name + "/find_fake_" + str(i) + "_" + str(i + stepsize) + ".log"
        os.system('bsub -M 2000 -n 2 -o {} \'python find_fake_prots.py  -s {} -e {}\''.format(log_file_name, int(i), int(i+stepsize)))
