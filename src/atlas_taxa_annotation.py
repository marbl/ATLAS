import os
import sys
import argparse

def main():

    parser = argparse.ArgumentParser(description="Reports for each read associated partition and taxonomy")
    parser.add_argument("-i","--partition_file", help="(*_partition_map_FINAL.txt)", required = True)
    parser.add_argument("-r","--read_partition", help="(*read_to_partition_assignment.txt)", required = True)
    parser.add_argument("-t","--database_taxonomy", help="database sequence to taxonomy map", required = True)
    parser.add_argument("-o","--output_taxa_file", help="The output file with sequence, assigned partition, and reduced set of taxa", required = True)
    args = parser.parse_args()

    partition_seqs = {}
    with open(args.partition_file) as f:
    	for line in f:
    		val = line.strip().split('\t')
    		if line.startswith('#'):
    			continue
    		if val[0] not in partition_seqs:
    			partition_seqs[val[0]] = {}
    		for dbseq in val[1].strip().split(';'):
    			partition_seqs[val[0]][dbseq.strip().split('.')[0]] = 1

    taxonomy = {}
    with open(args.database_taxonomy) as f:
    	for line in f:
    		val = line.strip().split('\t')
    		taxonomy[val[0]] = val[1]

    fw = open(args.output_taxa_file, 'w')
    fw.write('#sequence\tassigned_partition\tsmallest_set_of_taxa_in_partition\n')
    with open(args.read_partition) as f:
    	for line in f:
    		if line.startswith('#'):
    			continue
    		val = line.strip().split('\t')
    		read = val[0]
    		partition = val[2]
    		if partition in partition_seqs:
    			taxalist = [taxonomy[y] for y in partition_seqs[partition] if y in taxonomy]
    			taxalist_uniq = list(set(taxalist))
    			for t_iter in taxalist_uniq:
    				fw.write(read+'\t'+partition+'\t'+t_iter+'\n')

    fw.close()
    		# print (val)
    		# exit()

if __name__ == "__main__":
    main()
