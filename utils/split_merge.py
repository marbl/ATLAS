import os
import sys
import argparse
from itertools import groupby

def fasta_iter(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq

def main():
    parser = argparse.ArgumentParser(description="Splits a BLAST file and splits query sequence file")
    parser.add_argument("-q","--query_file",help="A fasta file of query sequences",required=True)
    parser.add_argument("-b","--blast_file", help="BLAST output file with output format",required=True)
    parser.add_argument("-t","--threshold",help="Number of query sequences per BLAST chunk", default=1000, required=False)
    parser.add_argument("-p","--prefix", help="Output prefix",default = '',required=False)
    args = parser.parse_args()
    
    current = "None"
    threshold = int(args.threshold)
    counter = 0
    current_num = 0
    queries = {}
    file_name = str(args.prefix) +'myseq'
    with open(args.blast_file) as f:
        for line in f:
            valname = line.strip().split('\t')[0]
            if current == "None":
                fw = open(''.join([file_name,'_',str(current_num),'.blast.out']), 'w')
                current = valname
                queries[valname] = ''.join([file_name,'_',str(current_num)])
                counter += 1
            if valname == current:
                fw.write(''.join([line.strip(),'\n']))
            else:
                if counter < threshold:
                    current = valname
                    queries[valname] = ''.join([file_name,'_',str(current_num)])
                    fw.write(''.join([line.strip(),'\n']))
                    counter += 1
                else:
                    current_num += 1
                    fw.close()
                    fw = open(''.join([file_name,'_',str(current_num),'.blast.out']), 'w')
                    fw.write(''.join([line.strip(),'\n']))
                    counter = 1
                    current = valname
                    queries[valname] = ''.join([file_name,'_',str(current_num)])


    query_files = list(set([queries[key] for key in queries]))
    file_handles = [open(''.join([filename,'.fasta']),'w') for filename in query_files]

    handle_map = {query_files[n]: file_handles[n] for n in range(len(query_files))}

    if args.query_file.lower().endswith(('.fasta', '.fa', '.fna')):
        fiter = fasta_iter(args.query_file)
        for ff in fiter:
            if ff[0] not in queries:
                continue
            file_handler = handle_map[queries[ff[0]]]
            file_handler.write(''.join(['>',ff[0],'\n',ff[1],'\n']))
    else:
        print ("Need a fasta file for reads")

if __name__ == '__main__':
    main()