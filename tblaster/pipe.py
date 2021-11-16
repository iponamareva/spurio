# Here, we store the function that takes care of everything related to the
# tblastn search (which is mayor step 1 of 3 in our pipeline).

def search_prot_against_db(n_entry,queryfasta, dbposition, querytype, overwrite, database_path, outputfolder):
    '''Take a query fasta, database and index of fasta entry of interest.
    This script will create a temporary fasta containing only one entry,
    and hand it over to tblastn, which will create a file tblastn_out_n.txt
    
    Currently, Spurio uses two databases: dbfilter_part1 and part2. Part one
    contains the genomes of mycoplasma and spiroplasma organisms, which need to
    be treated separately due to their different genetic code.

    If only one database is to be used, this file will need some changes.'''
    outputfolder = outputfolder + '/blast'

    import os
    from itertools import islice
    from Bio import SeqIO
    from Bio.Blast.Applications import NcbitblastnCommandline as ncl
    from . import blastXMLmerge

    # Load fasta file containing all queries that we ever want to look at.
    allqueries = SeqIO.parse(queryfasta,'fasta')

    # Write the n_th entry to helper file querytemp.fasta. querytemp.fasta contains only our 1 sequence of interest, and makes
    # it easier to work with it, since we don't have to deal with our huge queries.fasta anymore.
    SeqIO.write(next(islice(allqueries, n_entry, None))  , outputfolder+'/querytemp_'+querytype+'_'+str(n_entry)+'.fasta', 'fasta')

    #retrieve_path = outputfolder+'/querytemp_'+querytype+'_'+str(n_entry)+'.fasta'
    #cur_q = list(SeqIO.parse(retrieve_path,'fasta'))
    #print("THIS IS ID " + cur_q[0].id)
    #if len(cur_q[0].id.split('|')) > 1:
    #    prot_name = cur_q[0].id.split('|')[1]
    #    if reviewer(prot_name, PE_th, evi_dic) == 0:
    #        return 0

    # Define blast output files for the two searches against both of our databases.
    
    print("debugging in pipe py")
    print(str(n_entry))

    fullblast_outname_asn =  outputfolder+'/blastout_'+querytype+'_'+str(n_entry)+'.asn'
    fullblast_outname_xml =  outputfolder+'/blastout_'+querytype+'_'+str(n_entry)+'.xml'

    # This file is only produced for manual inspection
    outname = outputfolder+'/blastout_'+querytype+'_'+str(n_entry)+'.txt'

    # Clean Directory
    for outfile in [fullblast_outname_asn,fullblast_outname_xml, outname]:
        if os.path.exists(outfile) and overwrite:
            print('Removing old outputfile: ', outfile)
            os.system('rm {}'.format(outfile))
        elif os.path.exists(outfile) and not overwrite:
            print('Found existing output. Wont overwrite unless explicitly told to. Exiting.')
            sys.exit()

    # Database stuff has gone out of control. This is part of the path to both parts of our database (normal and mycoplasma).
    # That means, our 'database' thing that we created in executer is completely useless now.

    # dbposition_blank = 'db/database'
    dbposition_blank = database_path
    # not sure how to manage this part (ira)
   
    # We run mycoplasma bacteria separately, because they have a different genetic code. To correct the e-values,
    # we have to tell tblastn how large our database is in total. At the moment, I checked that manually and hardcoded it.
    # (evil!). For a single database, this argument is obsolete.

    nucs_db = int(5448272552)
    # зачем-то это все-таки нужно???


    ####################
    ### RUN TBLASTN! ###
    ####################

    print('    a) Running tblastn. This step can take several minutes.')

    # Search query against first database 
    tblastn_cline = ncl(cmd = 'tblastn',query = outputfolder+'/querytemp_'+querytype+'_'+str(n_entry)+'.fasta', db = dbposition_blank + '/database' , out = fullblast_outname_asn , outfmt= 11, max_target_seqs = 4000, evalue = 1)
    #tblastn_cline = ncl(cmd = 'tblastn',query = outputfolder+'/querytemp_'+querytype+'_'+str(n_entry)+'.fasta', db = dbposition_blank + '/database')


    print("debugging in pipe py 1")
    print(tblastn_cline)

    stdout, stderr = tblastn_cline()
    print("finished")

    # Reformat the results
    os.system('blast_formatter -archive {} -out {} -outfmt 5'.format(fullblast_outname_asn, fullblast_outname_xml))

    os.system('blast_formatter -archive {} -out {} -outfmt \'7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qseq sseq\''.format(fullblast_outname_asn, outname))

    # os.system('cat {} {} > {}'.format(outname+'main',outname+'myco',outname))
    # os.system('rm {} {}'.format(outname+'main',outname+'myco'))

    ###########################################################
    ### We do very strange stuff to merge our 2 xml files. ####
    ###########################################################

    ### If only one database was used, a single xml file has been produced 
    ### already and we don't need this step. 


def create_bedrequest(n_fasta, filename,n0=1,n1=8,n2=9):
    '''

    95 percent sure this function is not in use.


    Take a tblastn output (outputfmt 7), identify ORFS of interest and
    create a file 'bedrequests.bed', which can be used by bedtools to
    return the actual nucleotide sequences of our matches.
    '''

    import re, os, numpy, pdb

    with open(filename) as f:
        for line in f:
            # Remove \n, split lines by colum
            line = line[:-1].split()

            # If our line is not a comment line...
            if not line[0] == '#':

                # We determine what strand (+/-) our ORF lies on
                strand = '+' if int(line[n1])<int(line[n2]) else '-'
                start, stop = min(int(line[n1]),int(line[n2])),max(int(line[n1]),int(line[n2]))

                # And write ID, Start, Stop and Strand (+/-) to bedrequests.bed, which will later be used by 'bedrequest findFasta'
                os.system('printf \'{}\t{}\t{}\t{}\t{}\t{}\n\' >> {}/bedrequests_{}.bed'.format(line[n0],start-1,stop, 'dummy', 'dummy', strand, outputfolder, n_fasta))


def exec_bedrequest(n_fasta, bedfile, sourceposition):
    '''
    Also probably not in use. The bedtools story has been shifted to project.py

    Take bedrequests.bed file previously created by 'create_bedrequest'
    function. Run Bedtools to return a list of nucleotide sequences.
    '''

    import os

    os.system('bedtools getfasta -s -fi {} -bed {} -fo {}/nuc_seqs_{}.txt'.format(sourceposition, bedfile, outputfolder, str(n_fasta)))


###########################################################################################################
####################### RUN OUR PIPELINE ##################################################################
###########################################################################################################

#import os, sys
#assert (sys.version_info[0] >= 3), "Must be using Python 3"

# Define parameters of search
#nfasta, db_short, source_short, query_short = 0, 'bactall2', 'bactall', 'uniquery.fa'
#n_fasta, db_short, source_short, query_short = 0, 'bact_all', 'allbact', 'query.fasta'


# Set up paths, create folders
#query = '/nfs/research/bateman/hoeps/pipe/tblastn_bact/query/{}'.format(query_short)
#database = '/nfs/research/bateman/hoeps/pipe/tblastn_bact/db/{}/{}'.format(db_short, db_short)
#source = '/nfs/research/bateman/hoeps/pipe/tblastn_bact/db/{}/source/{}.fa'.format(db_short, source_short)
#for n_fasta in range(1):#
#    run_pipeline(n_fasta, db_short, query_short, query, database, source)
