This file contains a description of what we do to obtain spurious sequences for training.

Short version:
make_db.sh
translate_everything.py
make_queries.py
blastp_search.py (or parallel)
find_fake_prots.py (or parallel)
extract_fakes_from_list.py


1. Dowload .fasta file with genome of the organism of interest. For example, here the file is called is obtained from https://www.ebi.ac.uk/ena/data/view/GCA_000002945.2.

2. Download proteome, from https://www.uniprot.org/proteomes/

File constants.py describes the paths to this files : adjust if necessary. This file also contains exp_name, which will be the prefix of subsequent directories and files.

Run init.py: it will create the directories for logs and intermediate results.

3. Make a blastp database out of proteome sequences. We have a directory db here. Run make_db.py (and
adjust the paths in it if required)

4. Now we need to translate everything from the genome of our organism of
interest. Run translate_everything.py. By default, it creates 1 big file with
all the generated sequences. The number of generated sequences will be required later, so it is saved to the file num_sequences.py. Please, don't remove this file.

5.1. Now out of all generated sequences, we need to filter out those which have significant hits in our proteome database. First, we create lots of files for blastp seqrch by running make_queries.py. This may take a couple of minutes. If you work on cluster, it's better to do bsub "python3 make_queries.py"

5.2. Then we run blasp search: blastp_search.py This will take a long time.
If you are running it on cluster, you can run parallel_blastp.py

5.3. Then you have a lot of blastp output files in blastp_searsh/.

The idea is that we need to look into every output file and find significant hits there.
If there are some significant hits, we don't save the identifier of the sequence and move further. If the generated sequence does not look like a real protein in the database, we keep it.

Run find_fake_prots.py. Put the boundaries of your search there
(-s 0 -e 10000, for instanse) and run it. Also will take some time. Then:

Note: if you are running it on cluster, you can run
parallel_find_fake.py.

6. Then we take the table with our fake identifiers and generate a big output
file with al fake sequences.

extract_fakes_from_list.py

Now it will tell you where your outtput file with all the sequences is.

DONE

