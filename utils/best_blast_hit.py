import os 
import sys
#uery coverage constant
QCOV_CONST=80
# read file taxmap_ncbi_ssu_ref_nr99_138.1.txt
taxonomy = {}
with open(sys.argv[1]) as f:
	for num, line in enumerate(f):
		if line.startswith("primary"):
			continue
		val = line.strip().split('\t')
		taxa = val[-2]
		taxonomy[val[0]] = taxa

# read the blast file
seen = {}
fw = open(sys.argv[3], 'w') 
with open (sys.argv[2]) as f:
	for line in f:
		val = line.strip().split('\t')
		qcov = abs(int(val[7]) - int(val[6]))*100/int(val[-1])
		if qcov >= QCOV_CONST:
			if val[0] not in seen:
				seen[val[0]] = 1
				taxa = val[1].strip().split('.')[0]
				if taxa in taxonomy:
					fw.write(val[0]+'\t'+val[1]+'\t'+taxonomy[taxa]+'\n')

fw.close()
