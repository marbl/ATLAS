# Outlier detection in BLAST hits
The tool to detect which of the top BLAST hits are relevant to the query sequence. 

To run the code, you will need Python3 (tested on 3.5), [SciPy](https://www.scipy.org/), and [NumPy](http://www.numpy.org/). If you have Python2 environment, checkout python2_old branch and use it. Note that this is older, unoptimized code that does the task but can be slow.

If you consider using this tool, please cite our publication which describes the methods used for separating true BLAST matches from matches to similar but unrelated sequences.

[Shah, N., Altschul, S. F., & Pop, M. (2017). Outlier Detection in BLAST Hits. In LIPIcs-Leibniz International Proceedings in Informatics (Vol. 88). Schloss Dagstuhl-Leibniz-Zentrum fuer Informatik](http://drops.dagstuhl.de/opus/volltexte/2017/7651).
The journal version (extended paper) is [Shah, Nidhi, Stephen F. Altschul, and Mihai Pop. "Outlier detection in BLAST hits." Algorithms for Molecular Biology 13.1 (2018): 7](https://almob.biomedcentral.com/articles/10.1186/s13015-018-0126-3).
For any queries, please either ask on github issue page or send an email to Nidhi Shah (nidhi@cs.umd.edu).

### Running BLAST

The first step is to run [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) to search query sequences against a database. Because our tool parses BLAST output, we expect the BLAST output file in a specific format as shown in the following command. 
```
blastn -query query.fasta  -db database  -outfmt " 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen " -task megablast -out blast.out [-num_threads 8]
```
We advise not to use max_targe_seqs to limit the number of hits returned by BLAST because of its unintuitive behavior. Please read our [letter](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166) to Bioinformatics journal. 

### Finding outliers in BLAST hits
For each query sequence, our tool extracts top BLAST hits that it thinks are related to the query sequence; we refer to these relavant hits as outliers/candidate DB sequences here and in the publication. We can also output a subset of BLAST output i.e. keeping only those hits that are related to the query sequence. We also introduced a fallback mechanism to rely on percent identity cutoff in cases when our method was not able to separate related hits. 
```
python3 score_blast.py -h
usage: score_blast.py [-h] -q QUERY_FILE -b BLAST_FILE [-a RAISETO]
                      [-out OUTPUT_FILE] [-blast BLAST_OP]
                      [-max MAX_BLAST_HITS] [-qc QC_THRESHOLD]
                      [-pid PID_THRESHOLD]

A tool to decide which of the top BLAST hits are related to the query sequence

optional arguments:
  -h, --help            show this help message and exit
  -q QUERY_FILE, --query_file QUERY_FILE
                        A fasta file of query sequences
  -b BLAST_FILE, --blast_file BLAST_FILE
                        BLAST output file with output format: -outfmt " 6
                        qseqid sseqid pident length mismatch gapopen qstart
                        qend sstart send evalue bitscore qseq sseq "
  -a RAISETO, --raiseto RAISETO
                        A hyperparameter to control weight given to conserved
                        coloumns in the multiple sequence alignment step
                        (default value = 2.7)
  -out OUTPUT_FILE, --output_file OUTPUT_FILE
                        output file name (default - outlier.txt)
  -blast BLAST_OP, --blast_op BLAST_OP
                        output file name for subsetting BLAST hits (default -
                        subset_blast.txt)
  -max MAX_BLAST_HITS, --max_blast_hits MAX_BLAST_HITS
                        Maximum number of BLAST hits per query to be considered (default = 300)
  -qc QC_THRESHOLD, --qc_threshold QC_THRESHOLD
                        Minimum query coverage for the hit to qualify (value
                        between 0 and 1, default = 0.9)
  -pid PID_THRESHOLD, --pid_threshold PID_THRESHOLD
                        Select hits >= pid (percent identity) when we cannot
                        find candidate hits through MSA technique
                        (default=100)
```


### Annotating query sequences
Once we find candidate relavant BLAST hits (outliers) for each query sequence, we just transfer the lowest common ancestor (LCA) of their taxonomy annotation to the query sequence. We have written a simple script that takes output file from the previous step (outlier.txt) and the tab separated taxonomy of sequences in the database as input. (Check out example/database_taxonomy.tsv for the file format.)
```
python3 assign_taxon.py -h
usage: assign_taxon.py [-h] -t DB_FILE [-o OUTLIER_FILE] [-c COL_IND]
                       [-out OUTPUT_FILE]

A script that finds consensus taxonomy for the top BLAST hits that are related
to query sequence

optional arguments:
  -h, --help            show this help message and exit
  -t DB_FILE, --db_file DB_FILE
                        file containing taxonomic information for database
                        sequences; refer to the db_taxonomy file in example
                        folder
  -o OUTLIER_FILE, --outlier_file OUTLIER_FILE
                        file containing outlier BLAST hits for each query
                        sequence (default - outlier.txt)
  -c COL_IND, --col_ind COL_IND
                        column number of candidate DB seqs (outliers) default
                        = 1
  -out OUTPUT_FILE, --output_file OUTPUT_FILE
                        output file name (default - consensus_taxonomy.txt)
```
The output file contains taxonomic annotation for all query sequences that had at least one BLAST hit detected as outlier by our tool. 
<!---
Please refer to our WABI [paper](http://drops.dagstuhl.de/opus/frontdoor.php?source_opus=7651) and AMB [paper](https://almob.biomedcentral.com/articles/10.1186/s13015-018-0126-3) for understanding the algorithm.
--->
For any questions, please email nidhi@cs.umd.edu.

