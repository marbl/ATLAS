#This is a test script for phase 1: Outlier detection in BLAST htis
# Constructing BLAST database
makeblastdb -in database.fasta -out database.fasta -dbtype nucl

# Doing BLAST search
blastn -query query.fasta  -db database.fasta  -outfmt " 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen " -task megablast -out blast.out -num_threads 4

# Detecting outlier BLAST hits for query sequences
python ../src/score_blast.py -q query.fasta -b blast.out -out outlier.txt

#Assigning consensus taxonomic annotation to the query sequences
python ../src/assign_taxon.py -t database_taxonomy.tsv -o outlier.txt -out consensus_taxonomy_based_on_outliers.txt