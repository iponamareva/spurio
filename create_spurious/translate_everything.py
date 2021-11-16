import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Blast import NCBIXML
import os
from math import *
from tools import *
from constants import *

if not os.path.exists("output"):
    os.makedirs("output")

##########################################
### Extracting FASTA #####################
##########################################

chromosomes_fasta = SeqIO.to_dict(SeqIO.parse(genome_path, 'fasta'))

##########################################
### Translate everything #################
##########################################

translated_chroms = []
corresponding_nuc_seqs = []
gene_names = []
generated_seqs = []

for chrom in chromosomes_fasta:
    seq = chromosomes_fasta[chrom]

    translated_chroms.append(seq.seq.translate())
    translated_chroms.append(seq.seq[1:].translate())
    translated_chroms.append(seq.seq[2:].translate())
    translated_chroms.append(seq.seq.reverse_complement().translate())
    translated_chroms.append(seq.seq[:-1].reverse_complement().translate())
    translated_chroms.append(seq.seq[:-2].reverse_complement().translate())


for entry in translated_chroms:
    stops = [pos for pos, char in enumerate(entry) if char == STOP]

    for i in range(len(stops)):
        if i == 0:
            continue
        if (stops[i] - stops[i - 1] > AA_TH):
            generated_seqs.append(entry[stops[i - 1] + 1 : stops[i]])


print("Number of generated seqs : " + str(len(generated_seqs)))
NUM_seq_file = open('output/num_sequences.py', 'a')
print(exp_name + "_NUM = " + str(len(generated_seqs)), file=NUM_seq_file)

build_len_hist(generated_seqs, name="output/" + exp_name + "_all_generated", log=False, open=False, bins=60)
save_generated(generated_seqs, path="output/" + exp_name + "_all_generated_seqs.fasta", id_="all_seqs")


