import os
from itertools import islice
from Bio import SeqIO
from Bio.Blast.Applications import NcbitblastnCommandline as ncl
import num_sequences
from tools import check_dirs_for_exp
from constants import * 

NUM = getattr(num_sequences, exp_name + "_NUM")
print("NUM sequences retrieved:", NUM)

check_dirs_for_exp(exp_name)

output_dir = "blastp_search/" + exp_name + "/queries"
query_fasta = "output/" + exp_name + "_all_generated_seqs.fasta"

overwrite = True
allqueries = SeqIO.parse(query_fasta, 'fasta')

for n_entry in range(NUM):
    SeqIO.write(next(allqueries)  , output_dir+'/querytemp_'+str(n_entry)+'.fasta', 'fasta')
