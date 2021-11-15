import os
from itertools import islice
from Bio import SeqIO
from Bio.Blast.Applications import NcbitblastnCommandline as ncl
from myconstants import *

outputfolder = 'blastp_search/queries'
queryfasta = 'pombe_all_generated_seqs.fasta'

#path_fasta = 'entry.fasta-3.fasta'
overwrite = True

allqueries = SeqIO.parse(queryfasta,'fasta')

for n_entry in range(NUM):
    SeqIO.write(next(allqueries)  , outputfolder+'/querytemp_'+str(n_entry)+'.fasta', 'fasta')
