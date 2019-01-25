#!/usr/bin/env python3
"""
-------------------------------------------------------------------------------
This script generates relavant top BLAST hits, which we refer as outlier hits,
for each query sequence. 
Written by N. R. Shah [nidhi@cs.umd.edu] for WABI/AMB submission 2017.
-------------------------------------------------------------------------------
"""

import math
import argparse
from scipy import special
from itertools import groupby
import numpy as np

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

def calc_entropy(c_values, ctotal):
    if ctotal ==0:
        return 0
    prob_list = [float(a)/ctotal for a in c_values]
    entropy = 0
    for j in range(0, len(prob_list)):
        if prob_list[j] == 0.0:
            continue
        entropy += (-1)*prob_list[j]*math.log(prob_list[j],4)
    return entropy

def calc_col_score(partial_cumulative_sum, total_cumulative_sum, gamma_half_values, gamma_values):
    col_score = [gamma_half_values[x] for x in partial_cumulative_sum]
    col_score_sub = [gamma_half_values[0] for x in partial_cumulative_sum]
    return sum(col_score) - sum(col_score_sub) - gamma_values[total_cumulative_sum]

def execute(listofseqs, width, raiseto, gamma_half_values, gamma_values, listofpidnames):
    entropylist = []
    wholescore = 0
    a_list = np.cumsum([[1 if x == 'A' else 0 for x in row ] for row in listofseqs], axis = 0)
    c_list = np.cumsum([[1 if x == 'C' else 0 for x in row] for row in listofseqs], axis = 0)
    t_list = np.cumsum([[1 if x == 'T' else 0 for x in row] for row in listofseqs], axis = 0)
    g_list = np.cumsum([[1 if x == 'G' else 0 for x in row] for row in listofseqs], axis = 0)
    # To count all the characters in each column - just check the final row in the cumulative sum array
    c_total = a_list[-1, :] + c_list[-1, :] + t_list[-1, :] + g_list[-1, :] 
    entropy_list = [calc_entropy([a_list[-1, i], c_list[-1, i], t_list[-1, i], g_list[-1, i]],c_total[i]) for i in range(0, width)]
    colscore_new = [0 if ctotal_i == 0 else (entropy_list[i]**(raiseto))*calc_col_score([a_list[-1, i], c_list[-1, i], t_list[-1, i], g_list[-1, i]], ctotal_i, gamma_half_values, gamma_values) for i, ctotal_i in enumerate(c_total)]
    wholescore = sum(colscore_new)
    scorearray = []
    scorearray.append(wholescore)
    for k in range(1,len(listofseqs)):
        ya_sum = a_list[-1, :] - a_list[k-1, :]
        yc_sum = c_list[-1, :] - c_list[k-1, :]
        yt_sum = t_list[-1, :] - t_list[k-1, :]
        yg_sum = g_list[-1, :] - g_list[k-1, :]
        c_total_x = a_list[k-1, :] + c_list[k-1, :] + t_list[k-1, :] + g_list[k-1, :]
        c_total_y = c_total - c_total_x
        colscore_x = sum([0 if ctotal_i == 0 else (entropy_list[i]**(raiseto))*calc_col_score([a_list[k-1, i], c_list[k-1, i], t_list[k-1, i], g_list[k-1, i]], ctotal_i, gamma_half_values, gamma_values) for i, ctotal_i in enumerate(c_total_x)])
        colscore_y = sum([0 if ctotal_i == 0 else (entropy_list[i]**(raiseto))*calc_col_score([ya_sum[i] , yc_sum[i], yt_sum[i], yg_sum[i]], ctotal_i, gamma_half_values, gamma_values) for i, ctotal_i in enumerate(c_total_y)])
        score = colscore_x + colscore_y - wholescore
        scorearray.append(score)
        #This can even tell when query is sufficiently different than the rest i.e. when the peak in the scorearray is observed at index 1
        if len(scorearray) >= 2:
            if scorearray[-2] > scorearray[-1] and scorearray[-2] > scorearray[-3] and scorearray[-2] > 0:
            # This returns index in the scorearray with maximum score and the score array itself
                return len(scorearray)-2, scorearray
    return -1, scorearray

def write_summary(tuplescore, listofnames, output_file_outliers, output_file_summary, fblast, blast_lines, listofpidnames, ignore_flag):
    index, scorearray = tuplescore
    length = len(listofpidnames)
    if length == 0:
        listofpidnames.append('None')
    if index == 1:
        output_file_summary.write(''.join(['\t'.join([listofnames[0], 'NA', ';'.join(listofpidnames), listofnames[0], str(scorearray[index]), ignore_flag]),'\n']))
        if length != 0:
            output_file_outliers.write(''.join([listofnames[0], '\t', ';'.join(listofpidnames), '\n']))
            fblast.write(''.join(['\n'.join(blast_lines[0:length]),'\n']))
        return
    if index == -1:
    	# This is correctly selecting the right set of BLAST hits because the hits skipped by lower sequence coverage are not added to blast_lines and once we find a hit that is less than pid_threshold, we don't care about all hits occuring after that for that particular query sequence.
        output_file_summary.write(''.join(['\t'.join([listofnames[0], 'NA', ';'.join(listofpidnames), 'NA', 'NA', ignore_flag]),'\n']))
        if length != 0:
            output_file_outliers.write(''.join([listofnames[0], '\t', ';'.join(listofpidnames), '\n']))
            fblast.write(''.join(['\n'.join(blast_lines[0:length]),'\n']))
        return
    outliers = listofnames[1:index]
    output_file_summary.write(''.join(['\t'.join([listofnames[0], ';'.join(outliers), 'NA', 'NA', str(scorearray[index]), ignore_flag]),'\n']))
    output_file_outliers.write(''.join([listofnames[0], '\t', ';'.join(outliers), '\n']))
    fblast.write(''.join(['\n'.join(blast_lines[0:len(outliers)]), '\n']))
    return

def main():
    parser = argparse.ArgumentParser(description="A tool to decide which of the top BLAST hits are related to the query sequence")
    parser.add_argument("-q","--query_file", help="A fasta file of query sequences",required=True)
    parser.add_argument("-b","--blast_file", help="BLAST output file with output format: -outfmt \" 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq \"", required=True)
    parser.add_argument("-a","--raiseto", help="A hyperparameter to control weight given to conserved coloumns in the multiple sequence alignment step (default value = 2.7)", default=2.7, required=False)
    parser.add_argument("-out","--output_file", help="output file prefix (default - results)", default="results", required=False)
    parser.add_argument("-blast","--blast_op", help="output file name for subsetting BLAST hits (default - subset_blast.txt)", default="subset_blast.txt", required=False)
    parser.add_argument("-max","--max_blast_hits", help="Maximum number of BLAST hits per query (default = 300)", default=300, required=False)
    parser.add_argument("-qc","--qc_threshold", help="Minimum query coverage for the hit to qualify (value between 0 and 1, default = 0.9) ", default=0.9, required=False)
    parser.add_argument("-pid","--pid_threshold", help="Select hits >= pid (percent identity) when we cannot find candidate hits through MSA technique (default=100)", default=100, required=False)
    args = parser.parse_args()
    raiseto = float(args.raiseto)

    #Pre-computes gamma values 
    max_blast_hits = int(args.max_blast_hits)
    qc_threshold = float(args.qc_threshold)
    pid_threshold = float(args.pid_threshold)
    gamma_half_values = [special.gammaln(i+0.5) for i in range(0,max_blast_hits+5)]
    gamma_values = [special.gammaln(i+2) for i in range(0,max_blast_hits+5)]

    #Open file handles for summary and subset blast output
    output_file_outliers = open(str(args.output_file)+'_outliers.txt','w')
    output_file_outliers.write('#query_sequence\tcandidate_DB_seqs\n')
    output_file_summary = open(str(args.output_file)+'_outliers_summary.txt','w')
    output_file_summary.write('#query_sequence\tcandidate_DB_seqs(outliers)\tcandidate_DB_seqs_qualifying_percent_identity_cutoff\tquery_unrelated2DB\tscore_of_cut\tBLAST_hits_insufficient?\n')
    fblast = open(str(args.blast_op), 'w')

    #Read BLAST file
    queries, queries_len = fasta_iter(str(args.query_file))
    ignore_queries = {}
    ignore_flag = False
    # The lowest bitscore should at most be 0.8 of max bitscore 
    ignore_bitscore_threshold = 0.8
    listofseqs = []
    listofnames = []
    blast_lines = []
    current = 'None'
    max_bitscore = 0
    curr_bitscore = 0
    with open(str(args.blast_file)) as f:
        for line in f:
            val = line.strip().split('\t') 
            if current == 'None':
                current = val[0]
                query = queries[val[0]]
                query_name = val[0]
                listofseqs = []
                listofnames = []
                listofpidnames = []
                flagadd = True
                blast_lines = []
                listofseqs.append(list(query.upper()))
                listofnames.append(val[0])
                counter = 0
                max_bitscore = val[11]
            if val[0] != current:
                if float(curr_bitscore) < ignore_bitscore_threshold * float(max_bitscore):
                    ignore_flag = False
                    ignore_queries[current] = 1
                else:
                    ignore_flag = True
                tuplescore = execute(listofseqs, queries_len[query_name], raiseto, gamma_half_values, gamma_values, listofpidnames)
                write_summary(tuplescore, listofnames, output_file_outliers, output_file_summary, fblast, blast_lines, listofpidnames, str(ignore_flag))
                current = val[0]
                query = queries[val[0]]
                query_name = val[0]
                listofseqs = []
                listofnames = []
                listofpidnames = []
                flagadd = True
                blast_lines = []
                listofseqs.append(list(query.upper()))
                listofnames.append(val[0])
                counter = 0
            if val[0] == current and counter <= max_blast_hits:
                seqlen = abs(int(val[7])-int(val[6]))+1

                #Checking whether the BLAST hit at least covers 90% of the query length 
                if seqlen < qc_threshold * queries_len[query_name]:
                    continue
                if flagadd == True and float(val[2]) >= pid_threshold:
                    listofpidnames.append(val[1])
                else:
                    flagadd = False
                curr_bitscore = val[11]
                qseq = val[12]
                sseq = val[13]
                scopy = [i for num, i in enumerate(sseq) if qseq[num] != '-']
                head_spaces = int(val[6])-1
                trail_spaces = queries_len[query_name] - int(val[7])
                scopy_new = ['-']*head_spaces
                scopy_new.extend(scopy)
                scopy_new.extend(['-']*trail_spaces)
                listofseqs.append(scopy_new)
                listofnames.append(val[1])
                blast_lines.append(line.strip())
                counter += 1 

        #This is to take care of last query hits when the file ends
        tuplescore = execute(listofseqs, queries_len[query_name], raiseto, gamma_half_values, gamma_values, listofpidnames)
        write_summary(tuplescore, listofnames, output_file_outliers, output_file_summary, fblast, blast_lines, listofpidnames, str(ignore_flag))

    output_file_summary.close()
    output_file_outliers.close()
    fblast.close()

if __name__ == "__main__":
    main()




