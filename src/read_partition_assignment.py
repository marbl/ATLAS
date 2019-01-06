import os 
import argparse
import sys
from collections import defaultdict

def create_partition_map(dbtopartition_file):
    seq2partition = {}
    partition2seqs = {}
    high_num_unused = 0
    with open(dbtopartition_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            val = line.strip().split('\t')
            seq2partition[val[0]] = val[1]
            if val[1] in partition2seqs:
                partition2seqs[val[1]].append(val[0])
            else:
                partition2seqs[val[1]] = [val[0]]
            if int(val[1]) >= high_num_unused:
                high_num_unused = int(val[1]) + 1
    return seq2partition, partition2seqs, high_num_unused

def create_read_partition_map_first_flag(candidate_file, seq2partition, partition2seqs, high_num_unused, output_prefix):
    fw = open(output_prefix + '_read_to_partition_assignment.txt', 'w')
    fw.write('#partition assignment criteria: Assign read the parition to which the first candidate db sequence belongs in the outlier set (first_flag=True).\n')
    fw.write('#read\tdb_seqs_in_winning_partitions\tpartition_num(winner)\n')
    fs = open(output_prefix + '_read_to_partition_summary.txt', 'w')
    fs.write('#partition assignment criteria: Assign read the parition to which the first candidate db sequence belongs in the outlier set (first_flag=True).\n')
    fs.write('#read\tdb_seqs_in_winning_partitions\tpartition_num(winner)\tpartition_counts\tcandidate_db_seqs(outliers)\tsize_of_candidate_db_seqs_set\n')
    with open(candidate_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            val = line.strip().split('\t')
            candidates = val[1].split(';')
            if len(candidates) == 0:
                continue    # No outliers for this read
            counter = defaultdict(int)
            for c in candidates:    
                if c in seq2partition:
                    counter[seq2partition[c]] += 1
                else:
                    seq2partition[c] = str(high_num_unused)
                    partition2seqs[str(high_num_unused)] = [c]
                    high_num_unused += 1
                    counter[seq2partition[c]] += 1
            partition_counts = []
            for partition in counter:
                partition_counts.append(str(partition)+':'+str(counter[partition]))
            winner = seq2partition[candidates[0]]
            db_seqs_in_winner = ';'.join(partition2seqs[winner])
            fw.write(val[0] + '\t' + db_seqs_in_winner + '\t' + str(winner) + '\n')
            fs.write(val[0] + '\t' + db_seqs_in_winner + '\t' + str(winner) + '\t' + ';'.join(partition_counts) + '\t' + val[1] + '\t' + str(len(candidates)) + '\n')
    fw.close()
    fs.close()
    return seq2partition, partition2seqs

def create_read_partition_map_threshold(candidate_file, seq2partition, partition2seqs, high_num_unused, output_prefix, threshold):
    fw = open(output_prefix + '_read_to_partition_assignment.txt', 'w')
    fw.write('#partition assignment criteria: Assign read to the parition(s) that exceeds the ' + str(threshold) + ' proportion in the candidate (outlier) set (first_flag=False).\n')
    fw.write('#read\tdb_seqs_in_winning_partitions\tpartition_num(winner)\n')
    fs = open(output_prefix + '_read_to_partition_assignment_summary.txt', 'w')
    fs.write('#partition assignment criteria: Assign read to the parition(s) that exceeds the ' + str(threshold) + ' proportion in the candidate (outlier) set (first_flag=False).\n')
    fs.write('#read\tdb_seqs_in_winning_partitions\tpartition_num(winner)\tpartition_counts\tcandidate_db_seqs(outliers)\tsize_of_candidate_db_seqs_set\n')
    with open(candidate_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            val = line.strip().split('\t')
            candidates = val[1].strip().split(';')
            if len(candidates) == 0:
                continue    # No outliers for this read

            counter = defaultdict(int)
            for c in candidates:    
                if c in seq2partition:
                    counter[seq2partition[c]] += 1
                else:
                    seq2partition[c] = str(high_num_unused)
                    partition2seqs[str(high_num_unused)] = [c]
                    high_num_unused += 1
                    counter[seq2partition[c]] += 1
            
            winner = []
            len_candidates = len(candidates)
            partition_counts = []
            for partition in counter:
                partition_counts.append(str(partition)+':'+str(counter[partition]))
                if float(counter[partition])/len_candidates >= threshold:
                    winner.append(partition)
            db_seqs_in_winner = []
            for w in winner:
                db_seqs_in_winner += partition2seqs[w]
            db_seqs_in_winner_str = ';'.join(db_seqs_in_winner)
            winner_str = ';'.join(winner)
            partition_counts_str = ';'.join(partition_counts)
            if len(winner) > 0:
                fw.write(val[0] + '\t' + db_seqs_in_winner_str + '\t' + winner_str + '\n')
                fs.write(val[0] + '\t' + db_seqs_in_winner_str + '\t' + winner_str + '\t' + partition_counts_str + '\t' + val[1] + '\t' + str(len_candidates) + '\n')
            else:
                fs.write(val[0] + '\tNA\tNA\t' + partition_counts_str + '\t' + val[1] + '\t' + str(len_candidates) + '\n')
    fw.close()
    fs.close()
    return seq2partition, partition2seqs

def write_partition_map(partition2seqs, output_prefix):
    fw = open(output_prefix + '_partition_map_FINAL.txt', 'w')
    fw.write('#partition_num\tdb_seqs\tsize_of_partition\n')
    for partition in partition2seqs:
        fw.write(partition + '\t' + ';'.join(partition2seqs[partition]) + '\t' + str(len(partition2seqs[partition])) + '\n')
    fw.close()

def main():
    parser = argparse.ArgumentParser(description="Maps reads to the partition number")
    parser.add_argument("-c","--candidate_file", help="Reads to candidate db sequences i.e. outliers file (results_outliers.txt)", required = True)
    parser.add_argument("-p","--dbtopartition_file", help="Database sequences to partition assignment file", required = True)
    parser.add_argument("-t","--threshold", help="Consider partition only if the number of candidates in outlier set surpasses this threshold (between 0 and 1, default = 0.5). Any value above 0.5 will lead to majority partition. Can be interpreted as confidence threshold", default = 0.5, required = False)
    parser.add_argument("-first","--first_flag", help="If set to True, it will override the threshold selection criteria, and will output the partition of the best/first candidate db sequence", default = "False", required = False)
    parser.add_argument("-prefix","--output_prefix", help="Output prefix", default = "results", required = False)
    args = parser.parse_args()

    #process partition map
    seq2partition, partition2seqs, high_num_unused = create_partition_map(args.dbtopartition_file)

    #make partition assignment for reads
    if args.first_flag.startswith("T") or args.first_flag.startswith("t"):
        seq2partition, partition2seqs = create_read_partition_map_first_flag(args.candidate_file, seq2partition, partition2seqs, high_num_unused, args.output_prefix)
    else:
        seq2partition, partition2seqs = create_read_partition_map_threshold(args.candidate_file, seq2partition, partition2seqs, high_num_unused, args.output_prefix, float(args.threshold))
    #update the final partitions to DB sequence map
    write_partition_map(partition2seqs, args.output_prefix)
if __name__ == '__main__':
    main()