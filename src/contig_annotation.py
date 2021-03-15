import os
import sys
import argparse
import json

def get_LCA(taxidlist):
	
	answer = {}
	for i in taxidlist:
		val = i.strip().split(';')
		for num,name in enumerate(val):
			if num not in answer:
				answer[num] = {}
			if name == "NA":
				continue
			if name not in answer[num]:
				answer[num][name] = 0
			answer[num][name] += 1
	final = []
	for i in range(0, 7):
		if len(answer[i]) > 1:
			final.append("NA")
		elif len(answer[i]) == 0:
			final.append("NA")
		else:
			for x in answer[i]:
				final.append(x)
				break
	return ';'.join(final)

def get_annotation(tid, merged, sn, parent, level, levels):
	if tid not in parent:
		if tid in merged:
			tid = merged[tid]
			if tid not in parent:
				return None
		else:
			return None
	tidlist = [tid]
	while tid in parent:
		a = parent[tid]
		tidlist.append(a)
		tid = a
	result = {}
	for i in tidlist:
		if level[i] in levels:
			result[level[i]] = sn[i]
	final = []
	for i in ["superkingdom" , "phylum", "class" , "order" , "family", "genus" , "species"]:
		if i in result:
			final.append(result[i])
		else:
			final.append("NA")

	return ';'.join(final)


def main():
	parser = argparse.ArgumentParser(description="annotate longer sequences (contigs)")
	parser.add_argument("-s", "--outliers_file", help = "outliers.txt file for all genes in contigs", required=True)
	parser.add_argument("-b", "--blast_file", help="blast file to get subject id to taxid mapping", required=True)
	parser.add_argument("-m", "--merged_file", help = "merged.dmp file from NCBI taxonomy", required=True)
	parser.add_argument("-n", "--names_file", help="names.dmp file from NCBI taxonomy", required=True)
	parser.add_argument("-o", "--nodes_file", help="nodes.dmp file from NCBI taxonomy", required=True)
	parser.add_argument("-f", "--fetchmg_file", help="*all.marker_genes_scores.table file from fetchmg", default=None, required = False)
	parser.add_argument('-out', "--output_file", help="Output annotation json file", default="annotation.json", required=False)

	args = parser.parse_args()

	## First get taxonomy id for all subjects from BLAST file
	## TODO make this a user defined parameter
	TAXCOL = -1
	subject2id = {}
	with open(args.blast_file) as f:
		for line in f:
			val = line.strip().split('\t')
			if val[1] not in subject2id:
				subject2id[val[1]] = val[TAXCOL]

	outliers = {}
	#Read the outliers detected for the genes in contigs
	with open(args.outliers_file) as f:
		for line in f:
			val = line.strip().split('\t')
			if line.startswith('#'):
				continue
			if val[0] not in outliers:
				outliers[val[0]]= {}
				sID = val[1].strip().split(';')
				for x in sID:
					# print (x, "here")
					if x in subject2id:
						# print (x, "e")
						tid_x = subject2id[x]
						if tid_x not in outliers[val[0]]:
							outliers[val[0]][tid_x] = 0
						outliers[val[0]][tid_x] += 1


	
	### Read the merged.dmp names.dmp and nodes.dmp files now
	merged = {}
	# file loc - '/fs/cbcb-lab/mpop/projects/ilana_brito/annotation_taxonomy/quick_analysis_L82/ncbi_taxdump/merged.dmp'
	with open(args.merged_file) as f:
		for line in f:
			val = line.strip('\t|\n').split('\t|\t')
			merged[val[0]] = val[1].strip()
	sn = {} # scientific name
	#file loc - '/fs/cbcb-lab/mpop/projects/ilana_brito/annotation_taxonomy/quick_analysis_L82/ncbi_taxdump/names.dmp'
	with open(args.names_file) as f:
		for line in f:
			val = line.strip('\t|\n').split('\t|\t')
			if val[3].strip() == "scientific name":
				sn[val[0]] = val[1].strip()
	#read nodes.dmp file
	parent = {}
	level = {}
	#file loc- '/fs/cbcb-lab/mpop/projects/ilana_brito/annotation_taxonomy/quick_analysis_L82/ncbi_taxdump/nodes.dmp'
	with open(args.nodes_file) as f:
		for line in f:
			val = line.strip('\t|\n').split('\t|\t')
			if val[0].strip() == "1":
				level[val[0].strip()] = val[2].strip()
				continue
			parent[val[0].strip()] = val[1].strip()
			level[val[0].strip()] = val[2].strip()

	levels = {"superkingdom" : 1, "phylum": 1, "class":1 , "order": 1 , "family": 1, "genus":1 , "species":1}

	## Add marker gene annotation information to the output
	MARKER_GENE_FLAG = False
	if args.fetchmg_file != None:
		MARKER_GENE_FLAG = True
		marker_genes = {}
		with open(args.fetchmg_file) as f:
			for line in f:
				if line.startswith('#'):
					continue
				val = line.strip().split('\t')
				if val[0] not in marker_genes:
					marker_genes[val[0]] = []
				marker_genes[val[0]].append(val[2])



	## Start generating final json files now
	contigs = {}
	contig_taxidlist = {}
	for g in outliers:
		# splitting at last occurence of '_' to get the contig name
		contig_name = g.strip().rsplit('_', 1)[0]
		if contig_name not in contig_taxidlist:
			contig_taxidlist[contig_name] = {}
		if contig_name not in contigs:
			contigs[contig_name] = {}
			contigs[contig_name]["genes"]= {}
		if g not in contigs[contig_name]["genes"]:
			contigs[contig_name]["genes"][g] = {}
		
		taxidlist = {}
		for tid in outliers[g]:
			fulltaxa = get_annotation(tid, merged, sn, parent, level, levels)
			if fulltaxa not in taxidlist:
				taxidlist[fulltaxa] = 0
			if fulltaxa not in contig_taxidlist[contig_name]:
				contig_taxidlist[contig_name][fulltaxa] = 0
			contig_taxidlist[contig_name][fulltaxa]+= outliers[g][tid]
			taxidlist[fulltaxa] += outliers[g][tid]

		contigs[contig_name]["genes"][g]["LCA"] = get_LCA(taxidlist)
		if MARKER_GENE_FLAG:
			if g in marker_genes:
				contigs[contig_name]["genes"][g]["marker_genes"] = ';'.join(marker_genes[g])
			else:
				contigs[contig_name]["genes"][g]["marker_genes"] = None
		contigs[contig_name]["genes"][g]["taxa_annotation"] = {}
		for x in sorted(taxidlist,key=taxidlist.get, reverse=True):
			contigs[contig_name]["genes"][g]["taxa_annotation"][x] = taxidlist[x]
	key_order = ["LCA", "frequent_annotation", "taxa_annotation","genes"]
	final = {}
	final["sequences"] = {}
	for c in contigs:
		first=True
		contigs[c]["taxa_annotation"]= {}
		for x in sorted(contig_taxidlist[c],key=contig_taxidlist[c].get, reverse=True):
			if first==True:
				first = False
				contigs[c]["frequent_annotation"] = x
			contigs[c]["taxa_annotation"][x] = contig_taxidlist[c][x]

		contigs[c]["LCA"] = get_LCA(contig_taxidlist[c])
		final["sequences"][c] = {k: contigs[c][k] for k in key_order}

	with open(args.output_file, 'w') as outfile:
		json.dump(final, outfile, indent=2)



if __name__ == '__main__':
	main()