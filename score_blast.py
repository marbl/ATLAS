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
gamma_values = []
gamma_half_values = []
queries = {}
raiseto = 2.7

def fasta_iter(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip()
        hsplit = header.split(" ")
        seq = "".join(s.strip() for s in faiter.next())
        queries[hsplit[0]] = seq

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

def execute(a,c,t,g,listdb,width,raiseto):
	entropylist = []
	wholescore = 0
	asum = np.sum(a,axis=0)
	csum = np.sum(c,axis=0)
	tsum = np.sum(t,axis=0)
	gsum = np.sum(g,axis=0)

	for i in xrange(0,width):
		col_score = 0
		ctotal = asum[i] + csum[i] + tsum[i] + gsum[i]
		c_x = [0,0,0,0]
		c_x[0] = asum[i]
		c_x[1] = csum[i]
		c_x[2] = tsum[i]
		c_x[3] = gsum[i]
		entropy = calc_entropy(c_x,ctotal)
		entropylist.append(entropy)
		if ctotal == 0:
			col_score  = 0
		else:
			for j in xrange(0,4):
				col_score += gamma_half_values[c_x[j]]
				col_score -= gamma_half_values[0]
			col_score -= gamma_values[ctotal] 
			#scale with entropy
			scaled_col_score = (entropy**(raiseto))*col_score
			wholescore += scaled_col_score
	scorearray = []
	scorearray.append(wholescore)
	for k in xrange(1,a.shape[0]):
		score_x = 0
		score_y = 0
		xasum = np.sum(a[0:k,:],axis=0)
		xcsum = np.sum(c[0:k,:],axis=0)
		xtsum = np.sum(t[0:k,:],axis=0)
		xgsum = np.sum(g[0:k,:],axis=0)
		yasum = np.sum(a[k:,:],axis=0)
		ycsum = np.sum(c[k:,:],axis=0)
		ytsum = np.sum(t[k:,:],axis=0)
		ygsum = np.sum(g[k:,:],axis=0)
		for i in xrange(0, width):
			c_x = [0,0,0,0]
			c_y = [0,0,0,0]
			c_x[0] = xasum[i]
			c_x[1] = xcsum[i]
			c_x[2] = xtsum[i]
			c_x[3] = xgsum[i]
			c_y[0] = yasum[i]
			c_y[1] = ycsum[i]
			c_y[2] = ytsum[i]
			c_y[3] = ygsum[i]
			c_total_x = sum(c_x)
			c_total_y = sum(c_y)
			col_score_x = 0
			col_score_y = 0
			if c_total_x != 0 and c_total_y != 0:
				for j in xrange(0,4):
					col_score_x += gamma_half_values[c_x[j]]
					col_score_x -= gamma_half_values[0]
					col_score_y += gamma_half_values[c_y[j]]
					col_score_y -= gamma_half_values[0]
				col_score_x -= gamma_values[c_total_x] 
				col_score_y -= gamma_values[c_total_y]
				#scale with entropy
				entropy = entropylist[i]
				scaled_col_score_x = (entropy**(raiseto))*col_score_x
				scaled_col_score_y = (entropy**(raiseto))*col_score_y
				score_x += scaled_col_score_x
				score_y += scaled_col_score_y

			elif c_total_x == 0 and c_total_y != 0:
				for j in xrange(0,4):
					col_score_y += gamma_half_values[c_y[j]]
					col_score_y -= gamma_half_values[0]
				col_score_y -= gamma_values[c_total_y]
				entropy = entropylist[i]
				scaled_col_score_y = (entropy**(raiseto))*col_score_y
				score_y += scaled_col_score_y
			elif c_total_y == 0 and c_total_x != 0:
				for j in xrange(0,4):
					col_score_x += gamma_half_values[c_x[j]]
					col_score_x -= gamma_half_values[0]
				col_score_x -= gamma_values[c_total_x]
				entropy = entropylist[i]
				scaled_col_score_x = (entropy**(raiseto))*col_score_x
				score_x += scaled_col_score_x
		score = score_x + score_y - wholescore
		scorearray.append(score)
		if len(scorearray) >= 3:
			if scorearray[-2] > scorearray[-1] and scorearray[-2] > scorearray[-3] and scorearray[-2] > 0:
			# if scorearray[-2] > scorearray[-1] and scorearray[-2] > scorearray[-3]:
				return listdb[0:len(scorearray)-3], scorearray

	return [],scorearray

def main():
	parser = argparse.ArgumentParser(description="A tool to decide which of the top BLAST hits are related to the query sequence")
	parser.add_argument("-q","--query_file",help="A fasta file of query sequences",required=True)
	parser.add_argument("-b","--blast_file", help="BLAST output file with output format: -outfmt \" 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq \"",required=True)
	parser.add_argument("-a","--raiseto",help="A hyperparameter to control weight given to conserved coloumns in the multiple sequence alignment step (default value = 2.7)",default=2.7,required=False)
	parser.add_argument("-out","--output_file", help="output file name (default - outlier.txt)",default="outlier.txt",required=False)
	args = parser.parse_args()
	raiseto = float(args.raiseto)
	#Pre-computes gamma values 
	for i in xrange(0,104):
		gamma_half_values.append(special.gammaln(i+ 0.5))
		gamma_values.append(special.gammaln(i+2))

	fw = open(str(args.output_file),'w')
	#Read BLAST file
	fasta_iter(str(args.query_file))
	current = 'None'
	with open(str(args.blast_file)) as f:
		for line in f:
			val = line.strip().split('\t') 
			if current == 'None':
				current = val[0]
				query = queries[val[0]]
				listofdb = []
				a = np.zeros((101, len(query)), dtype=bool)
				c = np.zeros((101, len(query)), dtype=bool)
				t = np.zeros((101, len(query)), dtype=bool)
				g = np.zeros((101, len(query)), dtype=bool)
				sq = np.array(list(query))
				a[0,:] =  ((sq=='A')*1)
				c[0,:] =  ((sq=='C')*1)
				t[0,:] =  ((sq=='T')*1)
				g[0,:] =  ((sq=='G')*1)
				counter = 1
			if val[0] != current:
				tuplescore = execute(a[0:counter,:], c[0:counter,:], t[0:counter,:], g[0:counter,:], listofdb, len(query), raiseto)
				fw.write(str(current) + '\t'+  ';'.join(tuplescore[0]) +'\n') 
				current = val[0]
				query = queries[val[0]]
				listofdb = []
				a = np.zeros((101, len(query)), dtype=bool)
				c = np.zeros((101, len(query)), dtype=bool)
				t = np.zeros((101, len(query)), dtype=bool)
				g = np.zeros((101, len(query)), dtype=bool)
				sq = np.array(list(query))
				a[0,:] =  ((sq=='A')*1)
				c[0,:] =  ((sq=='C')*1)
				t[0,:] =  ((sq=='T')*1)
				g[0,:] =  ((sq=='G')*1)
				counter = 1
			if val[0] == current and counter < 101:
				seqlen = abs(int(val[7])-int(val[6]))+1
				#Checking whether the BLAST hit at least covers 90% of the query length 
				if seqlen < 0.9*len(query):
					continue
				qseq = val[12]
				sseq = val[13]
				qcopy = bytearray(qseq)
				scopy = bytearray (sseq)
				for i, e in reversed(list(enumerate(qseq))):
					if e == '-':
						del qcopy[i]
						del scopy[i]
				qseq = str(qcopy)
				sseq = str(scopy)	
				
				for i in xrange(0,int(val[6])-1):
					sseq = '-'+ sseq
				for i in xrange(len(sseq), len(query)):
					sseq = sseq + '-'
				sq = np.array(list(sseq))
				a[counter,:] =  ((sq=='A')*1)
				c[counter,:] =  ((sq=='C')*1)
				t[counter,:] =  ((sq=='T')*1)
				g[counter,:] =  ((sq=='G')*1)
				counter += 1 
				listofdb.append(val[1])
		tuplescore = execute(a[0:counter,:], c[0:counter,:], t[0:counter,:], g[0:counter,:], listofdb, len(query), raiseto)
		fw.write(str(current) + '\t'+  ';'.join(tuplescore[0]) +'\n')
	fw.close()

if __name__ == "__main__":
	main()




