import os
import sys
import argparse
import time
import subprocess
from subprocess import Popen, PIPE
def main():
    bin_dir = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(description="BLAST relevant hits detection and partition assignment pipeline")
    parser.add_argument("-q","--query_file", help="A fasta file of query sequences",required=True)
    parser.add_argument("-b","--blast_file", help="BLAST output file with output format: -outfmt \" 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq \"", required=True)
    parser.add_argument("-t","--db_file",help="file containing taxonomic information for database sequences; refer to the db_taxonomy file in example folder", required=True)
    parser.add_argument('-o','--dir',help='Output directory to put results',required=False,default='OD_output')
    parser.add_argument("-a","--raiseto", help="A hyperparameter to control weight given to conserved coloumns in the multiple sequence alignment step (default value = 2.7)", default="2.7", required=False)
    parser.add_argument("-max","--max_blast_hits", help="Maximum number of BLAST hits per query (default = 500)", default= "500" , required=False)
    parser.add_argument("-qc","--qc_threshold", help="Minimum query coverage for the hit to qualify (value between 0 and 1, default = 0.9) ", default= "0.9", required=False)
    parser.add_argument("-pid","--pid_threshold", help="Select hits >= pid (percent identity) when we cannot find candidate hits through MSA technique (default=100)", default= "100" , required=False)
    args = parser.parse_args()
    try:
      import networkx
    except ImportError:
      raise ImportError('Looks like you do not have networkx. Please rerun with networkx python module installed.')
      sys.exit(1)
    try:
      import community
    except ImportError:
      raise ImportError('Looks like you do not have networkx. Please rerun with community python module installed.')
      sys.exit(1)
    if not os.path.exists(args.dir):
        os.makedirs(args.dir)

    print (time.strftime("%c")+': Starting outlier detection step..', file = sys.stderr) 
    try:
        p = subprocess.check_output('python '+bin_dir+'/score_blast.py  -q '+args.query_file+' -b ' + args.blast_file + ' -a ' + args.raiseto + ' -blast ' + args.dir + '/subset_blast.txt -max ' + args.max_blast_hits + ' -qc ' + args.qc_threshold + ' -pid ' + args.pid_threshold + ' -out ' + args.dir + '/results', shell = True)
    except subprocess.CalledProcessError as err:
        print (time.strftime("%c")+': Failed to find outliers (relevant BLAST hits), terminating the process....\n' + str(err.output), file =sys.stderr)
        sys.exit(1)

    print (time.strftime("%c")+':Starting LCA based taxonomy assignment to outliers...', file = sys.stderr)
    try:
        p = subprocess.check_output('python '+bin_dir+'/assign_taxon.py  -t '+ args.db_file +' -o ' + args.dir + '/results_outliers.txt -out ' + args.dir + '/consensus_taxonomy_based_on_outliers.txt ', shell = True)
    except subprocess.CalledProcessError as err:
            print (time.strftime("%c")+': Failed to assign the LCA of relevant hits as taxonomic annotation to reads, terminating the process....\n' + str(err.output), file = sys.stderr)
            sys.exit(1)
    
    print(time.strftime("%c")+': Starting Stage2: Partitioning/Clustering database to create \"often confused together\" groups', file=sys.stderr)
    try:
        p = subprocess.check_output('python '+bin_dir+'/make_edge_list.py  -q '+ args.query_file +' -b ' + args.dir + '/subset_blast.txt -o ' + args.dir + '/results_outliers.txt -out ' + args.dir + '/edges.list', shell = True)
    except subprocess.CalledProcessError as err:
            print (time.strftime("%c")+': Failed to make edges of the co-occurence graph, terminating the process....\n' + str(err.output), file = sys.stderr)
            sys.exit(1)
    print(time.strftime("%c")+': Making partitions/clusters', file=sys.stderr)
    try:
        p = subprocess.check_output('python '+bin_dir+'/make_partitions.py  -e '+ args.dir + '/edges.list -o ' + args.dir + '/database_partition_map.txt', shell = True)
    except subprocess.CalledProcessError as err:
            print (time.strftime("%c")+': Failed to make database groups, terminating the process....\n' + str(err.output), file = sys.stderr)
            sys.exit(1)

    print(time.strftime("%c")+': Read to partition assignment', file=sys.stderr)
    try:
        p = subprocess.check_output('python '+bin_dir+'/read_partition_assignment.py  -c '+ args.dir +'/results_outliers.txt -p ' + args.dir + '/database_partition_map.txt -prefix '+ args.dir +'/results', shell = True)
    except subprocess.CalledProcessError as err:
            print (time.strftime("%c")+': Failed to make assignment of reads to partition number, terminating the process....\n' + str(err.output), file = sys.stderr)
            sys.exit(1)

    print(time.strftime("%c")+': Get taxonomic annotation for the reads based on partition assignment', file=sys.stderr)
    try:
        p = subprocess.check_output('python '+bin_dir+'/assign_taxon.py  -t '+ args.db_file +' -o ' + args.dir + '/results_read_to_partition_assignment.txt -out ' + args.dir + '/consensus_taxonomy_based_on_partition.txt ', shell = True)
    except subprocess.CalledProcessError as err:
            print (time.strftime("%c")+': Failed to get taxonomic annotations for reads, terminating the process....\n' + str(err.output), file = sys.stderr)
            sys.exit(1)

    print(time.strftime("%c")+': Get taxonomic composition of database partition', file=sys.stderr)
    try:
        p = subprocess.check_output('python '+bin_dir+'/assign_taxon.py  -t '+ args.db_file +' -o ' + args.dir + '/results_partition_map_FINAL.txt -out ' + args.dir + '/results_partition_taxonomy.txt ', shell = True)
    except subprocess.CalledProcessError as err:
            print (time.strftime("%c")+': Failed to get taxonomic annotations for partitions, terminating the process....\n' + str(err.output), file = sys.stderr)
            sys.exit(1)
    print(time.strftime("%c")+': Done...', file=sys.stderr)
if __name__ == '__main__':
    main()