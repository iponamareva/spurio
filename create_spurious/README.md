# Fake seq generator

The code here helps to generate fake protein sequences for a chosen organism.

### How it works

The idea here is that we take the genome of the organism and translate it in 6 reading frames. The sequences that we obtain this way contain amino acids as well as stop codons. We split resulting amino acid sequences between the stop codons and take only sequences which size is larger than the threshold, for instance, 40 amino acids. This is a reasonable threshold for a length of a protein sequence, however, it can be adjusted.

After translating everything and obtaining candidates for fake protein sequences, we need to search these sequences against a blast database, which is constructed from the organism's proteome. We filter out the sequences which have hits to existing proteins.

### Code structure

`constants.py` contains constants needed for a run. You need to specify the experiment name `exp_name` and paths to the genome and the proteome files.

- `init.py` created directories for outputs and intermediate results
- `make_db.py` creates a blast database out of the organism's proteome.
- 'translate_everything.py` translated the genome in 6 reading frames.
- `make_queries.py` creates queries for the tblast search
- `blastp_search.py` performs blast search. There's a parallel version for running on cluster.
- `find_fake_prots.py` inspects the results of the blast search and filters out the sequences with significant hits to real proteins.
- `extract_fakes_from_list.py` generates the final output file with the results.

These files need to be executed in the same order.

### Extra

- `num_sequences.py` contains the resulting number of generated sequences, some scripts need it to run searches and iterate through the files. Don't delete it.
  

### How to get the data

1. Dowload .fasta file with genome of the organism of interest. For example, here the file is called is obtained from https://www.ebi.ac.uk/ena/data/view/GCA_000002945.2.

2. Download proteome, from https://www.uniprot.org/proteomes/

File constants.py describes the paths to this files : adjust if necessary. This file also contains exp_name, which will be the prefix of subsequent directories and files.

### Details
1. Run init.py: it will create the directories for logs and intermediate results.

2. Put the experiment name and paths to the data in `constants.py`

3. Make a blastp database out of proteome sequences. Run `make_db.py`.

4. Now we need to translate everything from the genome of our organism of
interest. Run translate_everything.py. By default, it creates 1 big file with
all the generated sequences. The number of generated sequences will be required later, so it is saved to the file `num_sequences.py`. Please, don't remove this file.

5. Now out of all generated sequences, we need to filter out those which have significant hits in our proteome database. First, we create lots of files for blastp seqrch by running `make_queries.py`. This may take a couple of minutes. If you work on cluster, it's better to do bsub `python3 make_queries.py`.

6. Then we run blasp search: `blastp_search.py`. This will take a long time.
If you are running it on cluster, you can run `parallel_blastp.py`. In this case, make sure that all the jobs are finished before moving further.

7. Next, we need to look into every output file and find significant hits with tthe proteome database. If there are some significant hits, we don't save the identifier of the sequence and move further. If the generated sequence does not look like a real protein in the database, we keep it. Run `find_fake_prots.py`. Put the boundaries of your search there
(-s 0 -e 10000, for instanse) and run it. This step also will take some time. On cluster, you can run `parallel_find_fake.py`.

8. Then we take the table with our fake identifiers and generate a big output
file with al fake sequences by running `extract_fakes_from_list.py`. It will tell you where your outtput file with all the sequences is.

9. DONE

