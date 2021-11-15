'''This is the main script to create fancy pictures and prepare protein data
for further analysis.'''

SCALERNAME = 'models/MinMaxScaler_1111'

from sklearn.externals import joblib

import argparse
import os
import sys
import numpy as np
from datetime import datetime

import matplotlib as mpl
mpl.use('Agg')

import classifier.gpm_model.helperfunctions as helperfunctions
import classifier.gpm_model.gpm_predict as gpm_predict
from project.project import create_and_run_bedrequest, find_stops
from tblaster.pipe import search_prot_against_db
from project.tools import *
from constants import *

CURDATE = str(datetime.now())

# Assert Python version > 3.0
assert (sys.version_info[0] >= 3), "Must be using Python 3"

def main(argv):
    '''Take command line arguments for fasta sequence and evalue threshold.'''
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--start", type=int, default=0)
    parser.add_argument("-e", "--end", type=int, default=1)
    parser.add_argument("-v", "--evalue", type=int, default = 1)
    parser.add_argument("-r", "--source", default="db/source.fa")
    parser.add_argument("-q", "--query", type=str, default="myfiles/queries/query_0.txt")
    parser.add_argument("-qt", "--querytype", type=str, default='query')
    parser.add_argument("-db", "--database", type=str, default='db/database')
    parser.add_argument("-qn", "--querynumber", type=int, default=-1)
    parser.add_argument("-o", "--outputfolder", type=str, default='output_default')
    parser.add_argument("-mn", "--modelname", type=str, default='clf1.npy')
    args = parser.parse_args(argv)
    return args.start, args.end, args.evalue, args.query, args.source, args.querytype, args.database, args.querynumber, args.outputfolder, args.modelname


if __name__ == "__main__":
    start, stop, eval_cutoff, query, source, querytype, database, qn, OUTPUT_DIR, model_name = main(sys.argv[1:])


    CLF_DIR = 'classifier'
    BLAST_DIR = os.path.join(OUTPUT_DIR, "blast")
    BED_DIR = os.path.join(OUTPUT_DIR, "bedrequests")
    NUC_DIR = os.path.join(OUTPUT_DIR, "nuc_seqs")
    FIG_DIR = os.path.join(OUTPUT_DIR, "figs")
    NP_DIR = os.path.join(OUTPUT_DIR, "np")
    SUMM_DIR = os.path.join(OUTPUT_DIR, 'summaries_1')
    ERROR_DIR = os.path.join(OUTPUT_DIR, 'error_logs')
    STOPS_DIR = os.path.join(OUTPUT_DIR, 'stop_logs')
    MODEL_DIR = "models"

    print("*"*80)
    print("start time: " + str(datetime.now()))
    print("*"*80)

    
    DEFAULT_MODEL_PATH = os.path.join(MODEL_DIR, model_name)

    run_tblast, run_bed, run_histmaker, run_prediction = True, True, True, True
    overwrite_old_results = True

    error_logfile = os.path.join(ERROR_DIR, 'errors_summary_{}.log'.format(querytype))
    table_out_filename = os.path.join(SUMM_DIR, "table.txt")

    for directory in [OUTPUT_DIR, BLAST_DIR, BED_DIR, NUC_DIR, FIG_DIR, NP_DIR,
                      SUMM_DIR, ERROR_DIR, STOPS_DIR]:
        try:
            os.mkdir(directory)
        except FileExistsError:
            continue

    # Run pipeline: Search query through blast, use bedtools to extend results,
    # format and plot results ###
    

    for i in range(start, stop):
        iter_start_time = datetime.now()
        print("*"*80)
        print('Processing protein number {} from file {}'.format(i, query), flush=True)
        print("at  " + str(iter_start_time))
        print("*"*80)


        try:
            append_aa_end = 20  # We want pictures that show 20aa further than the stop of our protein

            # Define outputname and destination for tblastn result and nucleotide seqs.
            blast_out_filename = os.path.abspath(os.path.join(
                BLAST_DIR, "blastout_{}_{}.xml".format(querytype, i)))
            nuc_seq_out_filename = os.path.abspath(os.path.join(
                NUC_DIR, "nuc_seq_{}_{}.txt".format(querytype, i)))
            bed_out_filename = os.path.abspath(os.path.join(
                BED_DIR, "bedrequests_{}_{}.bed".format(querytype, i)))
            stops_out_filename = os.path.abspath(os.path.join(
                STOPS_DIR, "stops_log_{}_{}.txt".format(querytype, i)))

            # First big step:
            # Run our query against tblastn. (Time-intensive, around 1-3 min usually.)
            # We show our sequence to tblastn, which will return an xml and a text file
            # that both describe all homologues it finds in our genome database. These files
            # contain information about evalue, sequences, their position in the genome
            # and a couple more. To be found at project/blastouts/

            print("   searching protein against db")

            ###

            if run_tblast:
                search_prot_against_db(i, query, None, querytype, overwrite_old_results, database, OUTPUT_DIR)
            else:
                print('blast files already provided')
                
            # Intermediate step: define a few more paths and names

            qfasta_filename = os.path.abspath(os.path.join(BLAST_DIR, "querytemp_{}_{}.fasta".format(querytype, i)))

            with open(qfasta_filename) as qfasta_file:
                queryfastatemp = [line.rstrip() for line in qfasta_file]
            qlen = len(''.join(queryfastatemp[1:]))
            qname = queryfastatemp[0]

            # Second big step:
            # Use bedtools to extend every result to capture the full length + 20aa appendix
            # Short explanation: 'search_prot_against_db' has created an .xml file in project/blastouts/
            # which has information about where exactly in every genome lie the dna sequences of
            # interest. Our bedtool run takes these information and extracts these DNA sequences,
            # along with some neighbouring nucleotides (for more info, see the function).

            # Use bedtools to extend every result to capture the full length + 20aa appendix
            print("   creating and running bed request")

            if run_bed:
                create_and_run_bedrequest(i, blast_out_filename, source,
                                          bed_out_filename, qlen + append_aa_end,
                                          nuc_seq_out_filename)
            else:
                print('bed requestst already run')


            # Third big step
            # Format and plot results
            # We have everything we need now. Next, we have to glue our sequence fragments together,
            # find stop codons and draw a nice picture. See the function for more info.

            print("   finding the stops and saving data + building hists")
            if run_histmaker:
                print("    running histmaker")
                accs, X = find_stops(nuc_seq_out_filename, qname, i, qfasta_filename,
                           qlen, querytype, eval_cutoff, FIG_DIR, NP_DIR, stops_out_filename, table_out_filename)
            
            print("   running prediction")
            if run_prediction:
                print('    c) Running Gaussian classifier to predict spuriosity')
                clf_all = np.load(DEFAULT_MODEL_PATH)
                clf = clf_all[0]

                results_txtfile = os.path.join(SUMM_DIR,
                                               '{}_summary.txt'.format(querytype))

                #prob = gpm_predict.gpm_predict(X, accs, clf, min_, max_,
                #                               outlist=results_txtfile,
                #                               write_to_list=True)
                scaler = joblib.load(SCALERNAME)
                
                print('predicting spuriosity')
                prob = gpm_predict.gpm_predict(X, accs, clf, scaler,
                                               outlist=results_txtfile,
                                               write_to_list=True, prep=True)


                print('Protein {} (Query {}_{}) done!'
                      ' Spurious probability: {}'.format(accs[0],
                                                         querytype,
                                                         i,
                                                         np.round(prob[0], 6)))
                print('Results were saved to {}, {} and {}. \n'.format(FIG_DIR,
                                                                    NP_DIR,
                                                                    SUMM_DIR))


            iter_end_time = datetime.now()
            print("iteration took " + str(iter_end_time - iter_start_time))

        except RuntimeError as RTE:
            print('Error encountered with sequence {}. Skipping.'.format(i),
                  file=sys.stderr)
            with open(error_logfile, 'a') as logf:
                logf.write('Expected error in {} '
                           'number {}: {} \n'.format(querytype, i, RTE))
        except StopIteration:
            print(99)
        except Exception as e:
            print('Unexpected error. Stopping.', file=sys.stderr)
            with open(error_logfile, 'a') as logf:
                logf.write('Unexpected error in' 
                           '{} number {}: {} \n'.format(querytype, i, e))

