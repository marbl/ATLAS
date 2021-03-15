# Taxonomic classification of long sequences

If you want to annotate long sequences, such as contigs from metagenomic samples, we suggest the following pipeline. 
1. Predict genes or ORFs. We suggest using prodigal for this step. <br />
`prodigal -i [CONTIGS.fa] -p meta -a [CONTIGS_GENES.faa] -d [CONTIGS_GENES.fna] -o [CONTIGS_GENES.gbk]`
2. Predict which genes belong to 40 single copy prokaryotic marker genes. We use HMMs from fetchMG to identify which genes are marker genes. 
You can ignore this step if you don't want to use this information. <br />
`fetchMG.pl -m extraction [CONTIGS_GENES.faa] -o [CONTIGS_fetchmg] `
3. For all genes, search NCBI nt database, and find the best hits for each gene sequence using ATLAS. <br />
`blastn -query [CONTIGS_GENES.fna] -db [NCBI_nt/nt] -outfmt " 6 qseqid sseqid pident length mismatch gapopen qstart 
qend sstart send evalue bitscore qseq sseq qlen slen staxid" -out [CONTIGS_GENES.BLAST.out]` <br />
Please use staxid in the outfmt so it is easier to identify taxid for all hits. To find top BLAST hits for each gene sequence, run <br />
`python ATLAS/src/score_blast.py -q [CONTIGS_GENES.fna] -b [CONTIGS_GENES.BLAST.out] -qc 0.8 -pid 90 -out [results]`
4. We are going to use merged.dmp, names.dmp, and nodes.dmp file from NCBI taxonomy database. <br />
Download taxdump.tar.gz and unzip the folder - `wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz` <br />
 Next, get taxonomic annotation for all genes and contigs. <br />
`python ATLAS/src/contig_annotation.py -s [results_outliers.txt] -b [CONTIGS_GENES.BLAST.out] -m merged.dmp -n names.dmp -o nodes.dmp 
      -f CONTIGS_fetchmg/CONTIGS_GENES.all.marker_genes_scores.table -out annotation.json`
      
