import os
import argparse
import sys
from itertools import islice
from Bio import SeqIO
from Bio.Blast.Applications import NcbitblastnCommandline as ncl

outputfolder = 'blastp_search'
overwrite = True

#### First -- need to make a database, which is S.cerevisiae ###

def main(argv):
    '''Take command line arguments for fasta sequence and evalue threshold.'''
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--start", type=int, default=0)
    parser.add_argument("-e", "--end", type=int, default=10)
    args = parser.parse_args(argv)
    return args.start, args.end

if __name__ == "__main__":
    start, stop = main(sys.argv[1:])

    for n_entry in range(start, stop):

        fullblast_outname_asn =  outputfolder+'/blast'+'/blastout_'+str(n_entry)+'.asn'
        fullblast_outname_xml =  outputfolder+'/blast'+'/blastout_'+str(n_entry)+'.xml'

        outname = outputfolder+'/blastout_'+str(n_entry)+'.txt'

        for outfile in [fullblast_outname_asn, fullblast_outname_xml, outname]:
            if os.path.exists(outfile) and overwrite:
                print('Removing old outputfile: ', outfile)
                os.system('rm {}'.format(outfile))
            elif os.path.exists(outfile) and not overwrite:
                print('Found existing output. Wont overwrite unless explicitly told to. Exiting.')
                sys.exit()


        dbposition_blank = 'db/database'

        tblastn_cline = ncl(cmd = 'blastp',query = outputfolder+'/queries/querytemp_'+str(n_entry)+'.fasta', db = dbposition_blank, out = fullblast_outname_asn , outfmt= 11, max_target_seqs = 10000, evalue = 1, num_threads = 8)
        stdout, stderr = tblastn_cline()
        os.system('blast_formatter -archive {} -out {} -outfmt 5'.format(fullblast_outname_asn, fullblast_outname_xml))
        os.system('blast_formatter -archive {} -out {} -outfmt \'7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qseq sseq\''.format(fullblast_outname_asn, outname))


