import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Blast import NCBIXML
import os
from math import *


def build_len_hist(array, name='hist.png', log=False, open=False, bins=60, th=-1):
    
    ''' Builds Histogram of lengths'''

    if name != 'hist.png':
        name = name + "_hist.png"
    l = []
    for elem in array:
        ### should put it as parameter
        if (th != -1 and len(elem) < th):
            l.append(len(elem))
        elif th == -1:
            l.append(len(elem))

    if log:
        l = np.log10(l)

    plt.figure(facecolor='white')
    plt.hist(l, bins=bins)
    plt.savefig(name, dpi=200)
    if open==True:
        os.system('open ' + name)

def save_generated(array, accs=None, sort=False, path='saved_seqs.fasta', id_='>generated_seq_'):
    
    ''' Saves generated dequences from a list '''

    f = open(path, 'w')
    if id_ != '>generated_seq_':
        id_ = '>' + id_ + '_'
    if accs==None:
        if sort==True:
            array.sort(key=len)
        for i, seq in enumerate(array):
            print(id_ + str(i), file=f)
            for j in np.arange(0, len(seq), 80):
                print(seq[j:min(j + 80, len(seq))], file=f)
    else:
        if (len(array) != len(accs)):
            print('Different lenghts of array and acc array, exiting with exit code 1 from func save_generated')
            return 1

        for i in range(len(array)):
            print(id_ + str(accs[i]), file=f)
            seq = array[i]
            for j in np.arange(0, len(seq), 80):
                print(seq[j:min(j + 80, len(seq))], file=f)

    print("Saved all generated sequences, path=" + path)
    if accs!=None:
        print('Sequences were saved with provided acc numbers, not sorted')
    else:
        print('Sequences were enumerated automatically')

def batch_generator(array, mode='full', th1=-1, th2=-1, path='batch_default_path.txt', id_='>generated_seq_'):
    
    # Not used anymore

    if mode == 'full':
        save_generated(array, path=path, id_=id_)
        return 0
    if id_ != '>generated_seq_':
        id_ = '>' + id_ + '_'
    f = open(path, 'w')
    if mode == 'thresholds':
        i = 0
        array.sort(key=len)
        for seq in array:
            if th2 == -1:
                if len(seq) < th1:
                    continue
                i += 1
                print(id_ + str(i), file=f)
                for j in np.arange(0, len(seq), 80):
                    print(seq[j:min(j + 80, len(seq))], file=f)
    print("Saved requested batch in mode '" + mode + "', path=" + path)

