"""
-------------------------------------------------------------------------------
This script generates consensus taxonomic assignment based on selected top 
BLAST hits (outliers) for each query sequence.
Written by N. R. Shah [nidhi@cs.umd.edu] for WABI/AMB submission 2017.
-------------------------------------------------------------------------------
"""

import os
import sys
import numpy as np
import argparse
taxonomy = {} 
def main():

	parser = argparse.ArgumentParser(description="A script that finds consensus taxonomy for the top BLAST hits that are related to query sequence")
	parser.add_argument("-t","--db_file",help="file containing taxonomic information for database sequences; refer to the db_taxonomy file in example folder",required=True)
	parser.add_argument("-o","--outlier_file", help="file containing outlier BLAST hits for each query sequence (default - outlier.txt)",default="outlier.txt",required=False)
	parser.add_argument("-out","--output_file", help="output file name (default - consensus_taxonomy.txt) ",default="consensus_taxonomy.txt",required=False)
	args = parser.parse_args()

	fw = open(str(args.output_file),'w')

	#Read taxonomy
	with open(str(args.db_file)) as f:
		for line in f:
			l= line.strip()
			lsplit =  l.split('\t')
			taxonomy[lsplit[0]] = lsplit[1]
			
	#Read outlier file
	with open(str(args.outlier_file)) as f1:
		for line in f1:
			l = line.strip() 
			lsplit = l.split('\t')
			if len(lsplit) <= 1:
				continue
			else:
				outliers =  lsplit[1].split(';')
				names =[]
				for i in outliers:
					database_name = taxonomy[i].split(';')
					while len(database_name) < 7:
						database_name.append('')
					names.append(database_name)

				names1 = np.array(names)
				final_name = []
				if len(names1.shape) == 1:
					continue
				for i in xrange(0, names1.shape[1]):
					if len(np.unique(names1[:,i])) == 1:
						final_name.append(np.unique(names1[:,i])[0])	
					else:
						break
				fw.write(str(lsplit[0]) + '\t' + ';'.join(final_name)+'\n')
	fw.close()
				

if __name__ == "__main__":
	main()