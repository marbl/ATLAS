#!/usr/bin/env python3
import networkx as nx
import community
import argparse
import numpy as np

def process_graph(edge_list_file, threshold, output_file):
    G = nx.Graph()
    with open(edge_list_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            val = line.strip().split('\t')
            valsplit = val[0].strip().split(' ')
            if float(valsplit[0]) > 0:
                G.add_edge(valsplit[1],val[1],weight = float(valsplit[0]))
    #Open the partition output file 
    fw = open(output_file, 'w', 1)
    fw.write('#db_sequence\tpartition_number\n')
    #Find partitions in each connected component
    counter = 0 
    current_highest = 0
    subgraphs = list(nx.connected_component_subgraphs(G))
    for subgraph in subgraphs:
        part = community.best_partition(subgraph)

        #check for low-weight connections
        current_comm_nodes = list(part.keys()) 
        current_comm = G.subgraph(current_comm_nodes)
        degree_val = [x[1] for x in current_comm.degree(current_comm_nodes)]
        deg_med = np.median(degree_val)
        deg_std = np.std(degree_val)
        for db in part:
            if current_comm.degree(db) < (deg_med - 2 * deg_std) or current_comm.degree(db) < 2:
                print (degree_val, deg_med, deg_std, "partition just before-" + str(part[db] + counter))
                continue
            partition_num = part[db] + counter
            if partition_num > current_highest:
                current_highest = partition_num
            fw.write(db + '\t' + str(partition_num) + '\n')
        counter = current_highest + 1

def main():

    parser = argparse.ArgumentParser(description="Creates a graph and subsequently partitions of DB sequences that are often confused together")
    parser.add_argument("-e","--edge_list_file", help="edge list with weights", required = True)
    parser.add_argument("-t","--weight_threshold", help="Consider edges with weight greater than this threshold (default = 0)", default = 0, required = False)
    parser.add_argument("-o","--partition_file", help="The output file for Database sequence grouping", required = True)
    args = parser.parse_args()
    
    process_graph(args.edge_list_file, float(args.weight_threshold), args.partition_file)

if __name__ == "__main__":
    main()
