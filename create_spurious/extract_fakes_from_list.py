import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Blast import NCBIXML
import os
from math import *
from tools import *
from constants import *

command = "cat fake_files/{}* > output/{}_fake_proteins_table.tab".format(exp_name, exp_name)
print("Running concatenation:", command)
os.system(command)

fake_prot_table = "output/" + exp_name + "_fake_proteins_table.tab"
path_all_seqs = "output/" + exp_name + "_all_generated_seqs.fasta"

in_fasta_dic = SeqIO.to_dict(SeqIO.parse(path_all_seqs, 'fasta'))
in_size = len(in_fasta_dic)

descriptions, seqs, l = [], [], []
for i, key in enumerate(in_fasta_dic):
    descriptions.append(in_fasta_dic[key].description)
    seqs.append(in_fasta_dic[key].seq)
    l.append(len(in_fasta_dic[key].seq))

#### NOW WE NEED TO TAKE ONLY THOSE WHICH ARE FAKE ####
#### this table already exists, we made it before #####
#### we ran tblastn against all  sequences
#and threw away those who had matches with proteins from saccharomyces_cerevisiae
## now the table contains only those who does not have those hits

f = open(fake_prot_table)
fake_set = set()
for line in f:
    fake_set.add(line.split()[0])
print("Number of fake protein numbers: " + str(len(fake_set)))

f = open("output/" + exp_name + "_fake_proteins.fasta", 'w')
fake_seqs, fake_seqs_accs = [], []

for i in range(in_size):
    number = descriptions[i].split('_')[-1]
    if number in fake_set:
        seq = seqs[i]
        print(">fake_prot_" + str(number), file=f)
        for j in np.arange(0, len(seq), 80):
            print(seq[j:min(j + 80, len(seq))], file=f)

        fake_seqs.append(seq)
        fake_seqs_accs.append(number)

build_len_hist(fake_seqs, name="output/" + exp_name + "_fake_seqs")
save_generated(fake_seqs, accs=fake_seqs_accs, sort=False, path="output/" + exp_name + "_all_fake_seqs.fasta", id_='allfake_seq_')

#########################################
#### LET'S GENERATE A BATCH MANUALLY ####
#########################################

# print("total number of sequences: " + str(len(fake_seqs)))
# batch_f = open("batch_005.fasta", 'w')

# l = []
# for seq in fake_seqs:
#     l.append(len(seq))
# indexes = np.argsort(l)

# fake_seqs_accs = np.array(fake_seqs_accs)
# fake_seqs_accs = fake_seqs_accs[indexes]
# fake_seqs = np.array(fake_seqs)
# fake_seqs =fake_seqs[indexes]

# # Define coordinates of bins (pivots)
# N_BINS = 8
# pivots = []
# pivots.append(0)
# temp_thresholds = [60 + i*20 for i in range(N_BINS)]
# print("thresholds")
# print(temp_thresholds)
# #### take 25 out of each bin ######

# j = 0
# for i in range(len(fake_seqs)):
#     if len(fake_seqs[i]) > temp_thresholds[j]:
#         pivots.append(i)
#         j += 1
#         if (j == N_BINS):
#             break
# print("pivots found:")
# print(pivots)

# SAMPLE_SIZE = 25
# for i in range(N_BINS):
#     indexes = np.arange(pivots[i], pivots[i + 1], 1)
#     np.random.shuffle(indexes)

#     for i in range(SAMPLE_SIZE):

#         print('>fake_seq_' + str(fake_seqs_accs[indexes[i]]), file=batch_f)

#         seq = fake_seqs[indexes[i]]
#         for j in np.arange(0, len(seq), 80):
#             print(seq[j:min(j + 80, len(seq))], file=batch_f)


# last_index = pivots[len(pivots) - 1]

# for i in range(last_index, len(fake_seqs)):
#     print('>fake_seq_' + str(fake_seqs_accs[i]), file=batch_f)
#     seq = fake_seqs[i]
#     for j in np.arange(0, len(seq), 80):
#         print(seq[j:min(j + 80, len(seq))], file=batch_f)










