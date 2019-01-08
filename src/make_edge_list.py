#!/usr/bin/env python3
import os
import sys
import argparse
from itertools import groupby
import subprocess
def fasta_iter(fasta_name):
    queries = {}
    queries_len = {}
    with open(fasta_name, 'r') as fasta_file:
        faiter = (x[1] for x in groupby(fasta_file, lambda line: line[0] == ">"))
        for header in faiter:
            header = next(header)[1:].strip()
            hsplit = header.split(" ")
            seq = "".join(sequence.strip() for sequence in next(faiter))
            queries[hsplit[0]] = seq
            queries_len[hsplit[0]] = len(seq)
    return queries, queries_len

def get_non_qual_hits(blast_file, queries_len, qc_threshold, pid_threshold):
    non_qual_hits = {}
    with open(blast_file) as f:
        for line in f:
            val = line.strip().split('\t')
            seqlen = abs(int(val[7]) - int(val[6])) + 1
            if seqlen < qc_threshold * queries_len[val[0]]:
                non_qual_hits[(val[0],val[1])] = 1
                continue
            if float(val[2]) < pid_threshold:
                non_qual_hits[(val[0],val[1])] = 1
    return non_qual_hits

def process_edges(outlier_file, edge_output_file, non_qual_hits):
    fw = open(edge_output_file + '.tmp', 'w')
    # fw.write('#start_node(db_seq1)\tend_node(db_seq2)\n')
    with open(outlier_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            val = line.strip().split('\t')
            outliers = val[1].split(';')
            outliers_keep = [i for i in outliers if (val[0],i) not in non_qual_hits]
            if len(outliers_keep) > 1:
                fw.write('\n'.join([x+'\t'+y for x in outliers_keep for y in outliers_keep if x!=y]) + '\n')
    fw.close()
    cmd = 'echo "#weight_of_edge start_node(db_seq1)\tend_node(db_seq2)" > ' + edge_output_file + ' && sort ' + edge_output_file + '.tmp' + ' | uniq -c >> ' + edge_output_file
    try:
            p = subprocess.check_output(cmd, shell=True)
            p = subprocess.check_output('rm ' + edge_output_file + '.tmp', shell=True)
    except subprocess.CalledProcessError as err:
        print >> sys.stderr, str(err.output)
        sys.exit(1)
    

def main():

    parser = argparse.ArgumentParser(description="Creates a graph edge list of DB sequences where edge weight is # queries they are present together in the candidate set")
    parser.add_argument("-b","--subset_blast_file", help="The subsetted BLAST output (subset_blast.txt)", default = "subset_blast.txt", required = False)
    parser.add_argument("-q","--query_file", help="A fasta file of query sequences",required=True)
    parser.add_argument("-o","--outliers_file", help="The candidate DB sequences for queries (results_outliers.txt) file", required = True)
    parser.add_argument("-out","--edge_output_file", help="Edge list with weights = #co-occurence", default="edges.list", required=False)
    parser.add_argument("-pid","--pid_threshold", help="Consider hits >= pid (percent identity) only", default=95, required=False)
    parser.add_argument("-qc","--qc_threshold", help="Minimum query coverage for the hit to qualify (value between 0 and 1, default = 0.9) ", default = 0.9, required=False)
    args = parser.parse_args()
    queries, queries_len = fasta_iter(args.query_file)
    non_qual_hits = get_non_qual_hits(args.subset_blast_file, queries_len, float(args.qc_threshold), float(args.pid_threshold))
    process_edges(args.outliers_file, args.edge_output_file, non_qual_hits)
if __name__ == "__main__":
    main()
