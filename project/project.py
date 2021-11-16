'''This file contains most of the important functions and things.
It is imported by executer.py, and the functions are executed
throughout the run of executer.py.'''

import os
import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastformatterCommandline

TBLASTN_COVERED = 4
TBLASTN_GAP = 0
TBLASTN_START_CODON = 99
ANYWHERE_STOP_CODON = 300
FRONT_STOP = 1
BODY_STOP = 2
BACK_STOP = 3    
ALIGNMENT_BODY = 1
PSEUDOCOUNT = 1

stops_out_file = "stops_log.txt"

def create_and_run_bedrequest(n_fasta, blast_xml_filename, fasta_filename,
                              bed_out_filename, qlen, outname):
    '''
    Take a tblastn output (outputfmt 7), identify ORFS of interest and
    create a file 'bedrequests.bed', which can be used by bedtools to
    return the actual nucleotide sequences of our matches.

    Parameters
    ----------
    n_fasta : int
    blast_xml_filename : str
    fasta_filename : str
    bed_out_filename : str
    qlen : int
    outname : str
    '''

    # tblastn returns an alignment of the translated dna sequences it found. We, however,
    # want that alignment PLUS all the dna that is around it. (that is, everything that is
    # black in the pictures.) So, we take every result, check where it came from (e.g. position
    # 123456 - 123600 in genome XYZ.). Then, we extract a dna that comes right before and after
    # that, e.g. (123400 - 123465) and (123591 - 123640). Depending on how long the homologue
    # protein is compared to the query protein, these "appendices" are shorter or longer.

    # In mayor step 3, the appendices and the homologue proteins are then glued together to
    # receive our pictures.

    # We want appendices to overlap with our homologue protein, so we can check
    # everything is right.
    secoverlap = 3

    healthy_seqs = {}

    # Load xml files created by tblastn. The xml files contain the information
    # about where, who and what the homologue proteins are.
    NcbiblastformatterCommandline(archive=blast_xml_filename,
                                  outfmt='7 qseqid sseqid evalue qseq sseq',
                                  out="{}.txt".format(blast_xml_filename.split('.')[0]))
    # For every matched sequence...
    with open(blast_xml_filename) as blast_file:
        blast_records = NCBIXML.parse(blast_file)
        with open(bed_out_filename, "w") as bed_file:
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        # Update list of known sequences
                        healthy_seqs[alignment.title] = hsp.sbjct

                        # Extract the critical information needed for bedrequest

                        # a) where in the genome is the homologue?
                        nuc_start = hsp.sbjct_start
                        nuc_end = hsp.sbjct_end
                        
                        # Plus or minus strand?
                        strand = '+' if (hsp.frame[1] > 0) else '-'

                        # b) which parts of the query does this homologue cover?
                        # (e.g. amino acids 9 to 104)
                        q_start = hsp.query_start
                        q_end = hsp.query_end

                        startorig = min(nuc_start, nuc_end)
                        stoporig = max(nuc_start, nuc_end)

                        # How far do we have to reach back and forth, to cover the
                        # whole query protein (plus 20aa at the end?)
                        walkback = int(q_start) - 1
                        walkfront = qlen - int(q_end)
                        
                        # Also, we want our 3aa overlap.
                        extend = 1
                        if extend == 1:
                            if strand == '+':
                                start = startorig - (3 * walkback)
                                stop = stoporig + (3 * walkfront)
                            elif strand == '-':
                                start = startorig - (3 * walkfront)
                                stop = stoporig + (3 * walkback)

                        # and write ID, Start, Stop and Strand (+/-) to
                        # bedrequests.bed, which will later be used by
                        # 'bedtools getfasta ...'

                        adjusted_hit_id = '|'.join(alignment.hit_id.split('|')[1:])

                        bed_out_template = '{}\t{{}}\t{{}}\t{{}}\t{{}}\t{{}}'.format(adjusted_hit_id)
                        bed_out_lines = []
                        vals = [(start - 1, startorig - 1 + (3 * secoverlap)),
                                (stoporig - (3 * secoverlap), stop)]
                        tail = ["{}|{}:{}|{}|{}|{}".format(adjusted_hit_id,
                                                        q_start,
                                                        q_end,
                                                        hsp.query,
                                                        hsp.sbjct,
                                                        hsp.expect),
                                'dummy',
                                strand]
                        for val1, val2 in vals:
                            line = bed_out_template.format(val1,
                                                           val2, *tail)
                            bed_out_lines.append(line)
                        if strand == '-':
                            bed_out_lines.reverse()
                        
                        if (vals[0][0] > 0) and (vals[0][1] > 0):
                            for line in bed_out_lines:
                                print(line, file=bed_file)

    print('    b) Bedfiles written. Running bedtools to retrieve nuc seqs now.')
    cmd = ('bedtools getfasta -s -fi {} '
           '-bed {} '
           '-fo {} -name').format(fasta_filename, bed_out_filename, outname).split()
    subprocess.run(cmd, check=True)


def uniprot_name_from_fasta_description(desc):
    components = desc.split()[0].split("|")
    if len(components) >= 2:
        name = components[1].split("/")[0]
    else:
        if desc[0] == ">":
            name = components[0][1:]
        else:
            name = desc
    return name

def load_sequence_files(fasta_filename):
    """Load FASTA file.

    If the FASTA file contains an odd number of sequences for a particular
    identifier, they are all discarded.

    Parameters
    ----------
    str
        The filename of a FASTA file.

    Returns
    -------
    [Bio.SeqRecord.SeqRecord]
    """
    files = list(SeqIO.parse(fasta_filename, 'fasta'))
    if len(files) < 2:
        raise RuntimeError('Only 0 or 1 tblastn hits found.')

    # Changing name. This became necessary after switching from bedtools version 2.20 to 2.26
    for file_ in files:
        file_.id = file_.id.split('::')[0]

    # IRA: this id looks like this: >ENA|CP000034|CP000034.1|5:268|IKYFS----

    # Remove single entries. They stem from tblastn results too close to genome borders.

    # IRA : need to check here whether the pair of seqences have same coordinates!

    last_filename = None
    first_n_of_filename = 0
    n_file = 0
    bad_ranges = []
    for i, file_ in enumerate(files):
        filename = file_.id
        if filename == last_filename:
            n_file += 1
        else:
            if n_file % 2 != 0:
                print('Warning 0: Found an issue in one of the two '
                      ' appendices to sequence {}. Removing {} '
                      'sequences.'.format(i, first_n_of_filename - i))
                bad_ranges.append((first_n_of_filename, i))
            first_n_of_filename = i
            n_file = 1
        last_filename = filename
    bad_ranges.reverse()
    for start, end in bad_ranges:
        del files[start:end]

    # A: append our bedtool results to tblast output. This is a sensible step.
    # Before merging, we make several quality checks to assure that the
    # sequences look exactly as they should.

    assert len(files) % 2 == 0, ('Odd number of sequences. At least one lonely '
                                 'sequence was not captured.')
    return files

def combine_sequences(sequences):
    """Combine pairs of sequences found by TBLASTN.
    Parameters
    ----------
    [Bio.SeqRecord.SeqRecord]

    Returns
    -------
    [str]
    Raises
    ------
    RuntimeError
        If the sequences in the file can't be combined.
    """

    newseq = []
    for i in range(0, len(sequences), 2):
        file_ = sequences[i]
        next_file = sequences[i + 1]
        if len(file_) % 3:
            raise RuntimeError('Length of front appendix '
                               'is not a multiple of 3')

        if len(next_file) % 3:
            raise RuntimeError('Length of end appendix is not a multiple of 3')

        if file_.id != next_file.id:
            raise RuntimeError('Front and End appendix do '
                               'not have the same ID.')

        front = str(file_.seq.translate())
        body = file_.id.split('|')[-2]
        query_body = file_.id.split('|')[-3]
        back = str(next_file.seq.translate())
        front_tail = front[-3:]
        body_head = body[:3]
        body_tail = body[-3:]
        back_head = back[:3]

        if not (_seq_contains_gaps_or_pseudocode(front_tail)
                or _seq_contains_gaps_or_pseudocode(body_head)
                or front_tail == body_head):
            raise RuntimeError('Seq {}: Front appendix and MSA body do not'
                               ' overlap correctly. (what should be there {} '
                               'vs what tblastn claims {}), aka {} <- what is'
                               ' really in the genome(bedtools)'.format(
                                   i, front_tail, body_head, file_.seq[-9:]))
        if not (_seq_contains_gaps_or_pseudocode(body_tail)
                or _seq_contains_gaps_or_pseudocode(back_head)
                or body_tail == back_head):
            raise RuntimeError('Seq {}: End appendix and MSA body do not '
                               'overlap correctly. (what should be there {} '
                               'vs what tblastn claims {}), aka {}  <- what '
                               'is really in the genome(bedtools)'.format(
                                   i, back_head, body_tail, next_file.seq[0:9]))


        # If everything seems fine, merge! This is the important step #
        # Just to be clear, we are taking the homologue sequence found by
        # tblastn as-is, including gaps.
        # We then glue our appendices to it. The reason we do it so strangely
        # is that we need the gaps.
        newseq.append((file_.id, front[:-3], body, query_body, back[3:]))
    return newseq

def save_stops_data_init(stops_out_filename):

    os.system('printf \'Name\tHit start\tHit end\.......\n\' >> {}'.format(stops_out_filename))

def save_stops_data(acc_list, p_start, p_end, stops_array, seqlen, stops_out_filename):

    slen = len(stops_array)
    coords = np.nonzero((stops_array>=1) & (stops_array <= 3))    
    os.system('printf \'{}\t{}\t{}\t{}\n\' >> {}'.format("|".join(acc_list), p_start, p_end, "\t".join(map(str, coords)), stops_out_filename))

def make_visualisation_matrix(combined_sequences, qlen, stops_out_filename):
    n_seq = len(combined_sequences)

    save_stops_data_init(stops_out_filename)
    print("    successful: stop output file created: " + stops_out_filename, flush=True)
    
    pic = np.zeros((n_seq, qlen+20))
    pic_wo_gaps = np.zeros((n_seq, qlen + 20))
    pic_stops = np.zeros((n_seq, qlen+20))

    counter = 0
    all_stop_codons = []
    sseqlens = np.zeros(n_seq)
    starts = np.zeros(n_seq)
    stops = np.zeros(n_seq)
    evals = np.zeros(n_seq)

    for file_id, front_tail, body, query_body, back_head in combined_sequences:

        # _ is like this : "ENA", "LM993812", "LM993812.1"
        *_, protein_coords, seq, query_seq, evalue_s = file_id.split("|")

        # Fix formatting of evalues from bedtools in some circumstances
        evalue_s = evalue_s.split('(')[0]
        start, stop = map(int, protein_coords.split(':'))
        evals[counter] = float(evalue_s)
        if evals[counter] == 0:
            evals[counter] = 1e-200
        starts[counter] = start
        stops[counter] = stop

        pic[counter, start - 1:stop] = -np.log10(evals[counter])
        pic_wo_gaps[counter, start - 1:stop] = -np.log10(evals[counter])
        pic_stops[counter, start - 1:stop] = TBLASTN_COVERED

        # Appendix before seq
        front_stop_codons = [i for i, letter in enumerate(front_tail)
                             if letter == '*']
        front_start_codons = [i for i, letter in enumerate(front_tail)
                              if letter == 'M']

        # Appendix after seq
        back_stop_codons = [i + stop - 1 for i, letter in enumerate(back_head)
                            if letter == '*']
        back_start_codons = [i + stop - 1 for i, letter in enumerate(back_head)
                             if letter == 'M']

        # Actual seq
        sseqlen = len(seq)
        sseqlen_no_gaps = len(seq.replace('-', ''))
        qseqlen = (stop - start) + 1

        #body_stop_codons = [round((j / sseqlen) * qseqlen) + start
        #                    for j, letter in enumerate(body) if (letter == '*')]

        body_stop_codons = [round((j / sseqlen) * qseqlen) + start
                            for j, letter in enumerate(body) if (letter == '*' and query_body[j] != '-')]
        body_start_codons = [round((j / sseqlen) * qseqlen) + start
                             for j, letter in enumerate(body) if letter == 'M']
        body_gaps = [round((j / sseqlen) * qseqlen) + start
                     for j, letter in enumerate(body) if letter == '-']

        # for better pictures
        # body  ---NNNN--KKK    => ---NN--KKK
        # query NNNN--MMMMMM       NNNNMMMMMM
        # Stop codons, if they appear in body which is aligned to a gap in query, will be glued to the left part
        # TODO: write better explanation:)

        body_stop_wo = []
        body_gaps_wo = []
        body_start_wo = []
        current_coordinate = start

        for j, letter in enumerate(body):
            if (letter == '*'):
                if (len(body_stop_wo) == 0 or body_stop_wo[len(body_stop_wo) - 1] != current_coordinate):
                    body_stop_wo.append(current_coordinate)

            if letter == '-':
                if (len(body_gaps_wo) == 0 or body_gaps_wo[len(body_gaps_wo) - 1] != current_coordinate):
                    body_gaps_wo.append(current_coordinate)
            
            if letter == 'M':
                if (len(body_start_wo) == 0 or body_start_wo[len(body_start_wo) - 1] != current_coordinate):
                    body_start_wo.append(current_coordinate) 

            #  we count as coordinates only those position in the alignment which correspond to AA, not gaps
            #  so we change this coordinate of * in the picture only when * was aligned to AA

            if (query_body[j] != '-'):
                current_coordinate += 1

        all_stop_codons = front_stop_codons + back_stop_codons + body_stop_codons
        all_stop_codons_wo = front_stop_codons + back_stop_codons + body_stop_wo
        all_start_codons = front_start_codons + back_start_codons + body_start_codons
        all_start_codons_wo = front_start_codons + back_start_codons + body_start_wo

        # For our histograms, we only accept "post-alignment stops" if the
        # alignment before covered at least 75% of the query.
        #if (len(files[i+1].seq.translate()[3:])-20) < (0.25 * (len(file_.id.split('|',5)[-2])-20)):
        #    occs_late.append(back_stop_codons)

        pic[counter, body_gaps] = TBLASTN_GAP
        pic[counter, body_start_codons] = TBLASTN_START_CODON
        pic[counter, front_start_codons] = TBLASTN_START_CODON
        pic[counter, back_start_codons] = TBLASTN_START_CODON
        pic[counter, all_stop_codons] = ANYWHERE_STOP_CODON

        pic_wo_gaps[counter, body_gaps_wo] = TBLASTN_GAP
        pic_wo_gaps[counter, all_start_codons_wo] = TBLASTN_START_CODON
        pic_wo_gaps[counter, all_stop_codons_wo] = ANYWHERE_STOP_CODON

        pic_stops[counter, body_stop_codons] = BODY_STOP
        pic_stops[counter, back_stop_codons] = BACK_STOP
        pic_stops[counter, front_stop_codons] = FRONT_STOP

        sseqlens[counter] = sseqlen_no_gaps
        counter += 1
        all_stop_codons = []

        ##### Here we need to add lines to the file, which describes the results for this query
        ##### Id of the sequence, start and stop of the hit, our coordinates (all of them -- need to think in which order)
        ##### and coordinates of the stops (relative and absolute)
        
        #  IRA
        save_stops_data(_, start, stop, pic_stops[counter - 1], -1, stops_out_filename)
        # print("stops log saved", flush=True)

    ###################################################################################################################################
    ### Now that we have our merged sequences and identified stop codons and evalues, we do a bit of shuffling and annoying stuff.  ###
    ### It might even be that we dont need some of these lines anymore. I can look through it when I have time.                     ###
    ###################################################################################################################################

    # Reorder sequences by evalue
    evals = np.array(evals)
    evals_sort = np.argsort(-1 * evals)

    sorted_pic_wo = pic_wo_gaps[evals_sort]
    sorted_pic = pic[evals_sort]#[np.argsort(np.max(pic2, axis = i))]
    sorted_pic_stops = pic_stops[evals_sort]
    
    starts = starts[evals_sort]
    stops = stops[evals_sort]
    sseqlens = sseqlens[evals_sort]

    return evals, evals_sort, sorted_pic_wo, sorted_pic, sorted_pic_stops, starts, stops, sseqlens


def count_stops_per_aa(mat_stops, starts, stops):
    '''
    Input: The colorful MSA in matrix form
    Output: Number of Stop codons per amino acid.
    '''
    # Remove 20 aa appended to end
    mat_stops = mat_stops[:,:-20]
 
    # Number of amino acids to ignore at each start and end
    aa_tolerance = 10   
    
    # We don't even need to continue, if there are no homologous seqs or the 
    # sequence is smaller than 2 times the tolerance. 
    if ((np.shape(mat_stops)[0] == 0) or 
                      (np.shape(mat_stops)[1] <= 2*aa_tolerance)):
        return 0 + PSEUDOCOUNT, 0


    
    # Body_mask is an array of 0 and 1, where 1 marks amino acids in the body.
    body_mask = np.zeros(np.shape(mat_stops))
    for hom_seq in range(len(mat_stops)):
        body_mask[hom_seq, int(starts[hom_seq] + aa_tolerance - 1)
                         : int(stops[hom_seq] - aa_tolerance)] = ALIGNMENT_BODY
    
    mat_masked = mat_stops * body_mask
    
    n_stop = np.sum(mat_masked == BODY_STOP)
    n_bodyaas = np.sum(mat_masked > 0) + PSEUDOCOUNT
    stops_per_aa = n_stop / n_bodyaas
    
    return n_bodyaas, stops_per_aa


def add_table_entry(n, name, table_out_filename):
    os.system('printf \'{}\t{}\n\' >> {}'.format(n, name, table_out_filename))

def find_stops(inputfile, qname, n_fasta, qfasta, qlen, querytype, eval_cutoff, fig_dir, np_dir, stops_out_filename, table_out_filename):
    # Load nucleotide appendices. The MSA body is in the ID.
    # IRA: how to save it better?
    # input file is nuc_seq... file
    # format: >ENA|CP000034|CP000034.1|5:268|IKYFS----

    files = load_sequence_files(inputfile)
    newseq = combine_sequences(files)

    # pic2 -- new version of picture, where subjecr sequences are not normalized
    # pic1 -- old Wolfram's picture with normalization

    evals, evals_sort, pic2, pic1, pic_stops, starts, stops, sseqlens = make_visualisation_matrix(newseq, qlen, stops_out_filename)
    print("    made vis matrix", flush=True)

    stop_masked_pic = np.ma.masked_where(pic2 != ANYWHERE_STOP_CODON, pic2)

    # Mask out elements which don't contain stop codons and then set them to be grey
    start_masked_pic = np.ma.masked_where(pic2 != TBLASTN_START_CODON, pic2)
    start_masked_pic[np.ma.nonzero(start_masked_pic)] = 40

    pic2_onlycol = np.ma.masked_where((pic2 == TBLASTN_GAP) |
                                      (pic2 == ANYWHERE_STOP_CODON) |
                                      (pic2 == TBLASTN_START_CODON), pic2)


    # Get the "Spurious_ORF_XY" name right...
    if querytype == "af":
        try:
            a = np.load('tblaster/queries/antifam_seqs/dict_np.npy')
            truename = qname.split(' ')[0].split('|')[-1]
            spuriosity = a[np.where(a[:, 0] == truename)[0][0],1]
        except:
            spuriosity = 'unknown_af'
    else:
        spuriosity = 'noafhit'

    y_cutoff = np.argmin(np.absolute(eval_cutoff - (-np.log10(evals[evals_sort]))))

    ### FIGURE 1: Colorful MSA
    plt.figure(facecolor='black')
    ax = plt.subplot()
    ax.set_facecolor('black')
    plt.imshow(stop_masked_pic, interpolation='None', cmap='bwr', vmin=-3, vmax=100)
    plt.imshow(start_masked_pic, interpolation='None', cmap='binary', vmin=-3, vmax=100)
    plt.imshow(pic2_onlycol, interpolation='None', vmin=-3, vmax=100)
    if not y_cutoff == 0:
        plt.axhline(y=y_cutoff, color='white', linestyle='--', linewidth=1)
    cbar = plt.colorbar(ticks=range(0, 100, 10), boundaries=np.arange(-3, 100.001, 0.001))
    cbar.set_clim(-3, 100)
    cbar.set_label('Similarity to query sequence [-log10(E-value)]')
    plt.xlabel('AA position')

    idxs_above_th = [-np.log10(evals[evals_sort]) > eval_cutoff]
    pic_stops_cutoff = pic_stops[idxs_above_th]

    sseqlens = sseqlens[idxs_above_th]
    starts = starts[idxs_above_th]
    stops = stops[idxs_above_th]
    sseqlens = np.reshape(sseqlens, (len(sseqlens), 1))
    starts = np.reshape(starts, (len(starts), 1))
    stops = np.reshape(stops, (len(stops), 1))
    n_homologous_seqs = np.shape(pic_stops_cutoff)[0]
    if np.shape(pic_stops_cutoff)[0] == 0:
        raise RuntimeError('No sequences left after '
                           'evalue cutoff. Earlier: {}'.format(
                               np.shape(pic_stops)[0]))

    # Extract the number of stop codons in the region of interest
    body_aa, stops_per_aa = count_stops_per_aa(pic_stops_cutoff, starts, stops)

    #Save the colorful MSA
    file_name_template = "{}_{}_{}_{}.{{}}".format(querytype,
                                                   n_fasta,
                                                   spuriosity,
                                                   uniprot_name_from_fasta_description(qname))
    plt.savefig(os.path.join(fig_dir, file_name_template.format("png")),
                dpi=800)
    plt.close()

    add_table_entry(n_fasta, uniprot_name_from_fasta_description(qname), table_out_filename)

    # Save a much reduced picture, containing only 0s and 1,2,3 at STOP
    # positions before (1), within (2) or after (3) the 'body' alignment.
    # This is used as input for simple_classifier.
    np.save(os.path.join(np_dir, file_name_template.format("npy")),
            np.hstack((stops_per_aa, n_homologous_seqs, qlen, body_aa)))

    return np.array([uniprot_name_from_fasta_description(qname).split('_')[-1]]), np.array([stops_per_aa, n_homologous_seqs, qlen, body_aa])

def _seq_contains_gaps_or_pseudocode(seq):
    return len(set(seq) & {'X', '-', '*'}) > 0
