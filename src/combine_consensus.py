import os
import sys
'''
This program gives one annotation for each read/OTU. 
It reports annotation from partition analysis when available. 
When partition analysis fails to annotate the read, it reports consensus taxonomy based on outliers.
If that fails too, it just writes NA
'''
#read consensus based on outliers
outliers = {}
# sys.argv[1] is consensus_taxonomy_based_on_outliers.txt
with open(sys.argv[1]) as f:
	for line in f:
		val = line.strip().split('\t')
		if len(val) == 2 and val[1].strip() != "":
			outliers[val[0]] = val[1].strip()

fw = open(sys.argv[3], 'w')
#read consensus based on partitions
#sys.argv[2] is consensus_taxonomy_based_on_partition.txt
with open (sys.argv[2]) as f:
	for line in f:
		val = line.strip().split('\t')
		taxa = ""
		if len(val) != 2 or val[1].strip() == "":
			if val[0] in outliers:
				taxa = outliers[val[0]]
			else:
				taxa = "NA"
		else:
			taxa = val[1].strip()
		fw.write(val[0]+'\t'+taxa+'\n')

fw.close()
		
			
