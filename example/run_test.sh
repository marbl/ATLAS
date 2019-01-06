#This is a test script for full pipeline - outlier detection + database partitioning

# Constructing BLAST database
makeblastdb -in database.fasta -out database.fasta -dbtype nucl

# Doing BLAST search
blastn -query query.fasta  -db database.fasta  -outfmt " 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen " -task megablast -out blast.out -num_threads 4

#Check outputs in test_op directory, specifically check files consensus_taxonomy_based_on_outliers.txt, consensus_taxonomy_based_on_partition.txt, and results_read_to_partition_assignment.txt.
python ../src/run_pipeline.py -q query.fasta -b blast.out  -t database_taxonomy.tsv -o test_op