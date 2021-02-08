import json
import os
import sys
import argparse
import operator

def main():
    parser = argparse.ArgumentParser(description="writes all output in json object")
    parser.add_argument("-i","--partition_file", help="(*_partition_map_FINAL.txt)", required = True)
    parser.add_argument("-r","--read_partition", help="(*read_to_partition_assignment_summary.txt)", required = True)
    parser.add_argument("-b","--subset_blast", help="(*subset_blast.txt)", required = True)
    parser.add_argument("-co","--consensus_outliers", help="(*consensus_taxonomy_based_on_outliers.txt)", required = True)
    parser.add_argument("-cp","--consensus_partition", help="(*consensus_taxonomy_based_on_partition.txt)", required = True)
    parser.add_argument("-t","--database_taxonomy", help="database sequence to taxonomy map", required = True)
    parser.add_argument("-o","--output_json_file", help="The output file", required = True)
    args = parser.parse_args()
    database_taxa = {}
    with open(args.database_taxonomy) as f:
        for line in f:
            val = line.strip().split('\t')
            database_taxa[val[0]] = val[1].strip()
           

    partition2seqs = {}
    with open(args.partition_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            val = line.strip().split('\t')
            partition2seqs[val[0]] = val[1].strip().split(';')

    top_blast_hit = {}
    with open(args.subset_blast) as f:
        for line in f:
            val = line.strip().split('\t')
            if val[0] not in top_blast_hit:
                dbname = val[1].strip().split('.')[0]
                if dbname in database_taxa:
                    top_blast_hit[val[0]] = database_taxa[dbname]
                 
    lca_outliers = {}
    with open(args.consensus_outliers) as f:
        for line in f:
            val = line.strip().split('\t')
            if val[0] not in lca_outliers:
                if len(val) == 1:
                    lca_outliers[val[0]] = None
                    # print (val)
                else:
                    lca_outliers[val[0]] = val[1].strip()

    lca_partitions = {}
    with open(args.consensus_partition) as f:
        for line in f:
            val = line.strip().split('\t')
            if val[0] not in lca_partitions:
                if len(val) == 1:
                    lca_partitions[val[0]] = None
                else:
                    lca_partitions[val[0]] = val[1].strip()
    reads = {}                
    with open(args.read_partition) as f:
        for line in f:
            if line.startswith ('#'):
                # print(line)
                continue
            val = line.strip().split('\t')
            # print (val)
            if val[0] not in reads:
                reads[val[0]] = {}
                # Taxa based on top BLAST hit employing criteria for outlier detection
                if val[0] in top_blast_hit:
                    reads[val[0]]["taxa_top_hit"] = top_blast_hit[val[0]]
                else:
                    reads[val[0]]["taxa_top_hit"] = None
                ### Assign winning partition number to read and taxa based on partition sequences
                if val[2].strip().lower() != "na":
                    reads[val[0]]["assigned_partition"] = val[2].strip()
                    reads[val[0]]["taxa_partition"] = {}
                    p_taxa = val[1].strip().split(';')
                    curr = {}
                    for p in p_taxa:
                        p_label = database_taxa[p.strip().split('.')[0]]
                        if p_label not in curr:
                            curr[p_label] = 0
                        curr[p_label] += 1
                    for x in sorted(curr, key=curr.get, reverse=True):
                        reads[val[0]]["taxa_partition"][x] = curr[x]
                else:
                    reads[val[0]]["assigned_partition"] = None
                    reads[val[0]]["taxa_partition"] = None

                ### LCA of partition sequences
                if val[0] in lca_partitions:
                    reads[val[0]]["LCA_partition"] = lca_partitions[val[0]]
                else:
                    reads[val[0]]["LCA_partition"] = None
            #Now work with outliers and get 
                outliers = val[4].strip().split(';')
                curr= {}
                if val[4].strip().lower() != "na":
                    reads[val[0]]["taxa_outliers"] = {}
                    # reads[val[0]]["DB_outlier_sequences"] = outliers
                    for p in outliers:
                        p_label = database_taxa[p.strip().split('.')[0]]
                        if p_label not in curr:
                            curr[p_label] = 0
                        curr[p_label] += 1
                    for x in sorted(curr, key=curr.get, reverse=True):
                        reads[val[0]]["taxa_outliers"][x] = curr[x]
                else:
                    reads[val[0]]["taxa_outliers"] = None
                    # reads[val[0]]["DB_outlier_sequences"] = None
                ### LCA of outlier sequences
                if val[0] in lca_outliers:
                    reads[val[0]]["LCA_outliers"] = lca_outliers[val[0]]
                else:
                    reads[val[0]]["LCA_outliers"] = None
    ## Use partitions only if you want to output partitions in the final output. 
    # partitions = {}
    # for p in partition2seqs:
    #     if p not in partitions:
    #         partitions[p] = {}
    #         if partition2seqs[p] != None and len(partition2seqs[p]) != 0:
    #             partitions[p]["DB_sequences"] = partition2seqs[p]
    #             curr = {}
    #             for xp in partition2seqs[p]:
    #                 xp_label = database_taxa[xp.strip().split('.')[0]]
    #                 if xp_label not in curr:
    #                     curr[xp_label] = 0
    #                 curr[xp_label] += 1
    #             partitions[p]["taxa_partition"] = {}
    #             for xp in sorted(curr, key=curr.get, reverse=True):
    #                 partitions[p]["taxa_partition"][xp] = curr[xp]
    #         else:
    #             partitions[p]["DB_sequences"] = None
    #             partitions[p]["taxa_partition"] = None

    final = {}
    final["sequences"] = reads
    # final["partitions"] = partitions
    with open(args.output_json_file, 'w') as outfile:
        json.dump(final, outfile, indent=2)
if __name__ == '__main__':
    main()