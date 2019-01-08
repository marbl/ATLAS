## Outlier detection in BLAST htis: Finding relevant BLAST hits and sample-specific database partitioning pipeline
This code is used to find related BLAST hits to the query sequence without relying on ad-hoc criteria on bitscore, evalue, or percent identity. We construct a multiple sequence alignment (using pairwise alignment from BLAST) of the query sequence and the hits, and use a statistical scoring to find the relevant top hits. 

### Publication:
- [Shah, N., Altschul, S. F., & Pop, M. (2017). Outlier Detection in BLAST Hits. In LIPIcs-Leibniz International Proceedings in Informatics (Vol. 88). Schloss Dagstuhl-Leibniz-Zentrum fuer Informatik](http://drops.dagstuhl.de/opus/volltexte/2017/7651).
- The journal version (extended paper) is [Shah, Nidhi, Stephen F. Altschul, and Mihai Pop. "Outlier detection in BLAST hits." Algorithms for Molecular Biology 13.1 (2018): 7](https://almob.biomedcentral.com/articles/10.1186/s13015-018-0126-3).

-------------------------------------
Installation
-------------------------------------
```
  git clone git@github.com:shahnidhi/outlier_in_BLAST_hits.git
  cd outlier_in_BLAST_hits
```
Requirements:
-------------------
Before installing the software you need to make sure the following programs are installed on your machine. You can use --user option with pip3 if you don't have root access on the machine. It is usually easier to create a virtual environment with Python 3 and install these modules.

1. [Python: Version 3.x](https://www.python.org/downloads/) (tested on 3.7) 
2. [Scipy](https://www.scipy.org/index.html)
- pip3 install scipy
3. [NetworkX](https://github.com/networkx/networkx)
- pip3 install networkx
4. [python-louvain](https://github.com/taynaud/python-louvain)
- pip3 install python-louvain

To check if everything is correctly set up, run run_test.sh from example folder.
```
cd example
sh run_test.sh
cd ..
```
### Running BLAST

The first step is to run [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) to search query sequences against a database. Because our tool parses BLAST output, we expect the BLAST output file in a specific format as shown in the following command. 
```
blastn -query query.fasta  -db database  -outfmt " 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen " -task megablast -out blast.out [-num_threads 8]
```
We advise not to constrain BLAST search to limit the number of hits by using max_target_seqs/num_alignments parameter if you have enough computational resources, or keep at least a reasonable number (minimum 100 hits) that will capture all the hits you care about. 
<!--
Please read our [letter](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166) to Bioinformatics journal and their response [here](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty1026/5259186).)
-->

### Brief Summary 
There are two main components of our BLAST outlier detection program:
1. Find relevant BLAST hits (outliers) for the query sequences
<p> For each query sequence, our tool uses a Bayesian approach to extract the subset of top BLAST hits that are most closely related to the query sequence; we refer to this subset of top related hits as outliers/candidate DB sequences. This subset of BLAST results with the outlier/candidate sequences is outputted as a text file (subset_blast.txt). We also implement a percent identity cutoff as a fallback mechanism for cases when our method cannot identify a subset of closely-related hits. If you are interested in just running this step of the pipeline, look [here](https://github.com/shahnidhi/outlier_in_BLAST_hits/blob/master/outlier_detection_readme.md) </p>
2. Database partitioning 
<p> Here, we use significant outliers detected in the first step to construct a "confusion" graph where nodes are database sequences and edges are how many reads are closely related to both of them, i.e. number of times they are present together in an outlier set. We partition this graph to detect densely connected clusters/parititions. This allows us to characterize, in a sample-specific manner, the extent of taxonomic ambiguity within classification. In theory, the members (DB sequences)of a tightly knit partition are indistinguishable based on the reads in the dataset. It hopes to provides as finer resolution as possible without losing  accuracy. Thus, making reads to partition assignment is especially powerful in cases when you cannot make species-level classification, these partitions can often tell which few species are likely the origin. </p>

Taxonomic annotation of sample
<p> For taxonomic annotation, our program outputs the latest common ancestor (LCA) of the outliers outliers for each query sequence. Query sequences can also be analyzed in a taxonomically agnostic manner by just looking at the partition number assigned to each query sequence. </p>
  
### Run outlier detection and database partitioning pipeline. 
```
python3 src/run_pipeline.py -h
usage: run_pipeline.py [-h] -q QUERY_FILE -b BLAST_FILE -t DB_FILE [-o DIR]
                       [-a RAISETO] [-max MAX_BLAST_HITS] [-qc QC_THRESHOLD]
                       [-pid PID_THRESHOLD]

BLAST relevant hits detection and partition assignment pipeline

optional arguments:
  -h, --help            show this help message and exit
  -q QUERY_FILE, --query_file QUERY_FILE
                        A fasta file of query sequences
  -b BLAST_FILE, --blast_file BLAST_FILE
                        BLAST output file with output format: -outfmt " 6
                        qseqid sseqid pident length mismatch gapopen qstart
                        qend sstart send evalue bitscore qseq sseq "
  -t DB_FILE, --db_file DB_FILE
                        file containing taxonomic information for database
                        sequences; refer to the db_taxonomy file in example
                        folder
  -o DIR, --dir DIR     Output directory to put results
  -a RAISETO, --raiseto RAISETO
                        A hyperparameter to control weight given to conserved
                        coloumns in the multiple sequence alignment step
                        (default value = 2.7)
  -max MAX_BLAST_HITS, --max_blast_hits MAX_BLAST_HITS
                        Maximum number of BLAST hits per query (default = 500)
  -qc QC_THRESHOLD, --qc_threshold QC_THRESHOLD
                        Minimum query coverage for the hit to qualify (value
                        between 0 and 1, default = 0.9)
  -pid PID_THRESHOLD, --pid_threshold PID_THRESHOLD
                        Select hits >= pid (percent identity) when we cannot
                        find candidate hits through MSA technique
                        (default=100)
```
### How to interpret the output?
It generates a bunch of files in the output folder. These files can be useful for debugging or understanding different steps of the pipeline better. 
The files you will be interested in 
1. results_outliers.txt - a tsv with reads and relevant DB sequences
2. results_partition_map_FINAL.txt - partition number to DB sequence mapping
3. results_read_to_partition_assignment.txt - partition assignment for reads
4. consensus_taxonomy_based_on_outliers.txt - LCA of outliers
5. consensus_taxonomy_based_on_partition.txt - LCA of members of the partition assigned
6. subset_blast.txt - keeps only the relevant hits in the BLAST output

Note: This tool is under active development. For any queries, comment, or bug reports, please post them on Github issue page or send an email to Nidhi Shah (nidhi@cs.umd.edu).
