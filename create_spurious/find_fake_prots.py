import os
import sys
import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastformatterCommandline

NUM = 113776

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--start", type=int, default=0)
    parser.add_argument("-e", "--end", type=int, default=10)
    args = parser.parse_args(argv)
    return args.start, args.end

if __name__ == "__main__":
    start, stop = main(sys.argv[1:])
    
    path_list = 'fake_files/fake_proteins_table_' +str(start) + '.tab'
    f = open(path_list, 'a')
    
    for n_entry in range(start, stop):

        blast_xml_filename ='blastp_search/blast/blastout_' + str(n_entry) + '.xml'

        with open(blast_xml_filename) as blast_file:
            blast_records = NCBIXML.parse(blast_file)
            for blast_record in blast_records:
                if (len(blast_record.alignments) == 0):
                   print(str(n_entry) + '\t' + 'no_alignments', file=f)
                for alignment in blast_record.alignments:
                    if (len(alignment.hsps) == 0):
                            print(str(n_entry) + '\t' + 'no_hsps', file=f)
                        
