"""
-------------------------------------------------------------------------------
This script generates relavant top BLAST hits, which we refer as outlier hits,
for each query sequence. 
Written by N. R. Shah [nidhi@cs.umd.edu] for WABI/AMB submission 2017.
-------------------------------------------------------------------------------
"""

import os
import sys
import math
import argparse
from scipy import special
from itertools import groupby
import numpy as np

def fasta_iter(fasta_name):
    queries = {}
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        hsplit = header.split(" ")
        seq = "".join(s.strip() for s in next(faiter))
        queries[hsplit[0]] = seq
    return queries

def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

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

def calc_col_score(c, ctotal, gamma_half_values, gamma_values):
    col_score = [gamma_half_values[x] for x in c]
    col_score_sub = [gamma_half_values[0] for x in c]
    return sum(col_score) - sum(col_score_sub) - gamma_values[ctotal]

def execute(listofseqs, width, raiseto, gamma_half_values, gamma_values):
    entropylist = []
    print ( len(listofseqs))
    wholescore = 0
    a = [[1 if x == 'A' else 0 for x in row ] for row in listofseqs]
    c = [[1 if x == 'C' else 0 for x in row] for row in listofseqs]
    t = [[1 if x == 'T' else 0 for x in row] for row in listofseqs]
    g = [[1 if x == 'G' else 0 for x in row] for row in listofseqs]
    asum = np.sum(a,axis=0)
    csum = np.sum(c,axis=0)
    tsum = np.sum(t,axis=0)
    gsum = np.sum(g,axis=0)
    ctotal = asum + csum + tsum + gsum  
    entropylist = [calc_entropy([asum[i],csum[i], tsum[i], gsum[i]],ctotal[i]) for i in range(0, width)]
    lscore_new = [0 if ctotal_i == 0 else sum([1 for x in [asum[i], csum[i], tsum[i], gsum[i]]]) for i, ctotal_i in enumerate(ctotal)]
    colscore_new = [0 if ctotal_i == 0 else (entropylist[i]**(raiseto))*calc_col_score([asum[i], csum[i], tsum[i], gsum[i]], ctotal_i, gamma_half_values, gamma_values) for i, ctotal_i in enumerate(ctotal)]
    wholescore = sum(colscore_new)
    scorearray = []
    scorearray.append(wholescore)
    for k in range(1,len(listofseqs)):
        score_x = 0
        score_y = 0
        xasum = np.sum(a[0:k],axis=0)
        xcsum = np.sum(c[0:k],axis=0)
        xtsum = np.sum(t[0:k],axis=0)
        xgsum = np.sum(g[0:k],axis=0)
        yasum = np.sum(a[k:],axis=0)
        ycsum = np.sum(c[k:],axis=0)
        ytsum = np.sum(t[k:],axis=0)
        ygsum = np.sum(g[k:],axis=0)
        c_total_x = xasum + xcsum + xtsum + xgsum
        c_total_y = ctotal - c_total_x
        colscore_x = sum([0 if ctotal_i == 0 else (entropylist[i]**(raiseto))*calc_col_score([xasum[i], xcsum[i], xtsum[i], xgsum[i]], ctotal_i, gamma_half_values, gamma_values) for i, ctotal_i in enumerate(c_total_x)])
        colscore_y = sum([0 if ctotal_i == 0 else (entropylist[i]**(raiseto))*calc_col_score([yasum[i], ycsum[i], ytsum[i], ygsum[i]], ctotal_i, gamma_half_values, gamma_values) for i, ctotal_i in enumerate(c_total_y)])
        score = colscore_x + colscore_y - wholescore
        scorearray.append(score)
        #This can even tell when query is sufficiently different than the rest i.e. when the peak in the scorearray is observed at index 1
        if len(scorearray) >= 2:
            if scorearray[-2] > scorearray[-1] and scorearray[-2] > scorearray[-3] and scorearray[-2] > 0:
            # This returns index in the scorearray with maximum score and the score array itself
                return len(scorearray)-2, scorearray
    return -1, scorearray

def write_summary(tuplescore, listofnames, fw, fblast, blast_lines):
    index, scorearray = tuplescore
    if index == 1:
        fw.write(''.join(['\t'.join([listofnames[0], 'NA', listofnames[0], str(scorearray[index])]),'\n']))
        return
    if index == -1:
        fw.write(''.join(['\t'.join([listofnames[0], 'NA', 'NA', 'NA']),'\n']))
        return
    outliers = listofnames[1:index]
    fw.write(''.join(['\t'.join([listofnames[0], ';'.join(outliers), 'NA', str(scorearray[index])]),'\n']))
    fblast.write('\n'.join(blast_lines[0:len(outliers)]))
    return

def main():
    parser = argparse.ArgumentParser(description="A tool to decide which of the top BLAST hits are related to the query sequence")
    parser.add_argument("-q","--query_file",help="A fasta file of query sequences",required=True)
    parser.add_argument("-b","--blast_file", help="BLAST output file with output format: -outfmt \" 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq \"",required=True)
    parser.add_argument("-a","--raiseto",help="A hyperparameter to control weight given to conserved coloumns in the multiple sequence alignment step (default value = 2.7)",default=2.7,required=False)
    parser.add_argument("-out","--output_file", help="output file name (default - outlier.txt)",default="outlier.txt",required=False)
    parser.add_argument("-blast","--blast_op", help="output file name for subsetting BLAST hits (default - subset_blast.txt)",default="subset_blast.txt",required=False)
    parser.add_argument("-max","--max_blast_hits", help="Maximum number of BLAST hits per query (default = 300)",default=300,required=False)
    args = parser.parse_args()
    raiseto = float(args.raiseto)

    #Pre-computes gamma values 
    max_blast_hits = int(args.max_blast_hits)
    gamma_half_values = [special.gammaln(i+0.5) for i in range(0,max_blast_hits+5)]
    gamma_values = [special.gammaln(i+2) for i in range(0,max_blast_hits+5)]

    #Open file handles for summary and subset blast output
    fw = open(str(args.output_file),'w')
    fw.write('#query_sequence\tcandidate_DB_seqs(outliers)\tquery_unrelated2DB\tscore_of_cut\n')
    fblast = open(str(args.blast_op), 'w')

    #Read BLAST file
    queries = fasta_iter(str(args.query_file))

    listofseqs = []
    listofnames = []
    blast_lines = []
    current = 'None'
    with open(str(args.blast_file)) as f:
        for line in f:
            val = line.strip().split('\t') 
            if current == 'None':
                current = val[0]
                query = queries[val[0]]
                ## New code
                listofseqs = []
                listofnames = []
                blast_lines = []
                listofseqs.append(list(query.upper()))
                listofnames.append(val[0])
                counter = 0
            if val[0] != current:
                print (current, len(listofseqs))
                tuplescore = execute(listofseqs, len(query), raiseto, gamma_half_values, gamma_values)
                write_summary(tuplescore, listofnames, fw, fblast, blast_lines)
                current = val[0]
                query = queries[val[0]]
                listofseqs = []
                listofnames = []
                blast_lines = []
                listofseqs.append(list(query.upper()))
                listofnames.append(val[0])
                counter = 0
            if val[0] == current and counter <= max_blast_hits:
                seqlen = abs(int(val[7])-int(val[6]))+1

                #Checking whether the BLAST hit at least covers 90% of the query length 
                if seqlen < 0.9*len(query):
                    continue

                qseq = val[12]
                sseq = val[13]
                scopy = [i for num, i in enumerate(sseq) if qseq[num] != '-']
                head_spaces = int(val[6])-1
                trail_spaces = len(query) - int(val[7])
                scopy_new = ['-']*head_spaces
                scopy_new.extend(scopy)
                scopy_new.extend(['-']*trail_spaces)
                listofseqs.append(scopy_new)
                listofnames.append(val[1])
                blast_lines.append(line.strip())
                counter += 1 

        #This is to take care of last query hits when the file ends
        tuplescore = execute(listofseqs, len(query), raiseto, gamma_half_values, gamma_values)
        write_summary(tuplescore, listofnames, fw, fblast, blast_lines)

    fw.close()
    fblast.close()

if __name__ == "__main__":
    main()




