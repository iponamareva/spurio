This file contains a description of what we do to obtain spurious sequences for training.

1. Dowload .fasta file with genome of the organism of interest. Our file is called
entry.fasta-3.txt and is obtained from
https://www.ebi.ac.uk/ena/data/view/GCA_000002945.2.

2. Download proteome, Our proteome is uniprot-filtered-reviewed_S_pombe.fasta.

3. Make a blastp database. We have a directory db here. Run ./make_db.sh (and
adjust the paths in it if required)

4. Now we need to translate everything from the genome of our organism of
interest. Run translate_everything.py. By default, it creates 1 big file with
all the generated sequences. Please, note that you need the number of
generated sequences for further analysis.

5.1. Now we need to filter out the sequences, which have significant hits in our
database. First, do the following:

mkdir fake_files
mkdir blastp_search
cd blastp_search
mkdir blast
mkdir queries

Then, we create lots of files for blastp seqrch by running
python make_queries.py
(may take a couple of minutes)

5.2. Then we run blasp search: blastp_search.py
This will take a long time.
If you are running it on cluster, you can run parallel_blastp.py

5.3. Then you have a lot of blastp output files in blastp_searsh/blast/.

Then we need to look into every output file and find significant hits there.
If there are some significant hits, we don't save the identifier.

Run find_fake_prots.py. Put the boundaries of your search there
(-s 0 -e 10000, for instanse) and run it. Also will take some time. Then:

cd fake_files
cat * > fake_proteins_table.tab
mv fake_proteins_table.tab ../

Note: if you are running it on cluster, you can simply run
parallel_find_fake.py.

If you had many files in your folder, it will concatenate them. If not, you
still need to do this :)

6. Then we take the table with our fake identifiers and generate a bit output
file with al fake sequences.

extract_fakes_from_list.py

Now fake_proteins.fasta is your outtput file with all the sequences.

DONE

