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
    current_community_id = 0
    next_community_id = 0
    subgraphs = list(nx.connected_component_subgraphs(G))
    for subgraph in subgraphs:
        part = community.best_partition(subgraph)
        communities = list(set(part.values()))
        for c_ind in communities:
            current_community_id = next_community_id
            current_comm_nodes = [x for x in part if part[x] == c_ind]
            current_comm = nx.Graph(G.subgraph(current_comm_nodes))
            current_comm_edges = current_comm.edges(data=True)
            edge_wts = [x[2]['weight'] for x in current_comm_edges]
            edge_med = np.median(edge_wts)
            edge_std = np.std(edge_wts)
            to_remove_edges = []
            # print ("Before", current_community_id, len(list(nx.connected_component_subgraphs(current_comm))))
            for e in current_comm_edges:
                if e[2]['weight'] == 1 or e[2]['weight'] < (edge_med - 2 * edge_std):
                    to_remove_edges.append((e[0], e[1]))
            current_comm.remove_edges_from(to_remove_edges)
            sub_sub_communitites = list(nx.connected_component_subgraphs(current_comm))
            for sub_sub_ind in sub_sub_communitites:
                if len(sub_sub_ind.nodes()) == 1:
                    continue
                else:
                    current_community_id = next_community_id
                    for nodes_ind in sub_sub_ind.nodes():
                        fw.write(nodes_ind + '\t' + str(current_community_id) + '\n')
                        next_community_id = current_community_id + 1
            # print ("Before", current_community_id, [len(x.nodes()) for x in list(nx.connected_component_subgraphs(current_comm))])
            # for db in current_comm_nodes:
            #     if current_comm.degree(db, 'weight') == 0:
            #         continue
            #     fw.write(db + '\t' + str(current_community_id) + '\n')
            #     next_community_id = current_community_id + 1

def main():

    parser = argparse.ArgumentParser(description="Creates a graph and subsequently partitions of DB sequences that are often confused together")
    parser.add_argument("-e","--edge_list_file", help="edge list with weights", required = True)
    parser.add_argument("-t","--weight_threshold", help="Consider edges with weight greater than this threshold (default = 0)", default = 0, required = False)
    parser.add_argument("-o","--partition_file", help="The output file for Database sequence grouping", required = True)
    args = parser.parse_args()
    
    process_graph(args.edge_list_file, float(args.weight_threshold), args.partition_file)

if __name__ == "__main__":
    main()
