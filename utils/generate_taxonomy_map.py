import os
import sys
from ete3 import NCBITaxa
'''
This program creates the database_taxonomy.tsv file that ATLAS needs to assign consensus taxonomic labels
Input to this program is taxmap_ncbi_ssu_ref_nr99_138.1.txt, and the output is database_taxonomy.tsv
'''
ncbi = NCBITaxa()
levels= ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

seen = {}
fw = open('database_taxonomy.tsv', 'w')
# sys.argv[1] is taxmap_ncbi_ssu_ref_nr99_138.1.txt
with open(sys.argv[1]) as f:
	for num,line in enumerate(f):
		if num == 0:
			continue

		val = line.strip().split('\t')
		accession = val[0]
		if accession in seen:
			continue
		seen[accession] = 1
		taxmap = {}
		taxaval = val[3].strip().split(';')
		taxaval_new = [x.strip().split('(')[0].strip().split('<')[0].strip() for x in taxaval if x.strip() != ""]
		taxid = ncbi.get_name_translator(taxaval_new)
		rev_name = {}
		for k,v in taxid.items():
			for viter in v:
				rev_name[viter] = k
		taxid_rank = ncbi.get_rank([x for sublist in taxid.values() for x in sublist])
		rev_rank = dict((v,k) for k,v in taxid_rank.items())

		list2print = []
		for l in levels:
			if l in rev_rank:
				if l == "species":
					sp = '_'.join(rev_name[rev_rank[l]].strip().split()[0:2])
					list2print.append(sp)
				else:
					list2print.append(rev_name[rev_rank[l]].strip().replace(' ','_'))
			else:
				list2print.append('NA')
		fw.write(accession+'\t'+';'.join(list2print)+'\n')

fw.close()
