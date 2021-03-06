#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

# Importation des modules

import argparse
import sys
import os
import statistics
import random
from random import randint
import networkx as nx
from networkx import algorithms
import matplotlib.pyplot as plt

__author__ = "Hocine Meraouna & Hager Elharty"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Hocine Meraouna & Hager Elharty"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Hocine Meraouna & Hager Elharty"
__email__ = "hocine.meraouna@gmail.com & hager.elharty@outlook.fr"
__status__ = "Developpement"


# I) CREATION DU GRAPHE de DE BRUIJN
##  a) Identification des kmer unique

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()

def read_fastq(fastq):
    '''
    A function that takes a fastq file and return an iterator
    of a given sequence
    :Parameters:
          fastq : fastq file
    Returns: generator of a given sequence
    '''
    fastq_file = open(fastq)
    lines = iter(fastq_file.readlines())
    for line in lines :
        yield next(lines)
        next(lines)
        next(lines)

def cut_kmer(seq, kmer):
    '''
    a function that takes a sequence and the length of the k-mer
    :Parameters:
          seq : a sequence as str given by
          kmer : length of the k-mer
    :Returns: a generator of the k-mers
    '''
    seq = seq.strip('\n')
    for j in range(len(seq) - kmer + 1):
        yield seq[j:j+kmer]

def build_kmer_dict(fastq, km_len):
    """
    a function that takes a fastq file, the k-mer length and returns
    a dictionary having as keys the k-mer and as values the number of
    each k-mer iteration
    :Parameters:
          fastq : fastq file
          km_len : length of the k-mer
    Returns: A dictionary that contains the number of each k-mer iteration
    """
    dico = {}
    for seq in read_fastq(fastq):
        for k_mer in cut_kmer(seq, km_len):
            if k_mer not in dico:
                dico[k_mer] = 1
            else :
                dico[k_mer] += 1
    return dico

## b) Construction de l’arbre de de Bruijn

def build_graph(dico_kmer):
    """
    A function that builds an nx graph using a k-mer dictionary
    :Parameters:
          dico_kmer : dictionary that contains the number of each k-mer iteration
    Returns: the k-mers tree (graph) the nodes will be the k-mers and the weights
            will be the k-mer iterations
    """
    graph = nx.DiGraph()
    for k_mer in dico_kmer:
        node1 = k_mer[:-1]
        node2 = k_mer[1:]
        graph.add_edge(node1 , node2 , weight = dico_kmer[k_mer])
    return graph


def show_graph(graph):
    """
    A function that shows the nx graph
    :Parameters:
          graph : a networkx digraph graph
    """
    nx.draw(graph, pos=nx.spring_layout(graph))
    #nx.draw_networkx_labels(graph, nx.spring_layout(graph))
    plt.draw()
    plt.show()
 
# II)Parcours du graphe de de Bruijn

def get_starting_nodes(graph):
    '''
    A function that creates a list of entry nodes of a graph
    :Parameters:
          graph : nx graph
    Returns: a list of entry nodes
    '''
    list_entre = []
    for node in graph :
        pred = list(graph.predecessors(node))
        if not pred :
            #print("Pas de predecesseur\n")
            list_entre.append(node)
    return list_entre

def get_sink_nodes(graph):
    '''
    A function that creates a list of sink nodes of a graph
    :Parameters:
          graph : nx graph
    Returns: a list of sink nodes
    '''
    list_sink = []
    for node in graph :
        succ = list(graph.successors(node))
        if not succ :
            #print("Pas de predecesseur\n")
            list_sink.append(node)
    return list_sink

def get_contigs(graph, lst_start, lst_end):
    """
    A function that creates a tuple of contigs and their length
    :Parameters:
          graph : nx graph
          lst_start : list of entry nodes
          lst_end : list of sink nodes
    Returns: a tuple of contigs and their length
    """
    contigs = []
    for source in lst_start :
        for target in lst_end :
            if algorithms.has_path(graph, source, target):
                path = algorithms.shortest_path(graph, source, target)
                contig = path[0]
                for i in range(len(path)-1):
                    contig += path[i+1][-1]
                contigs.append((contig, len(contig)))
    return contigs

def fill(text, width=80):
    """
    Split text with a line return to respect fasta format
    """
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contigs, out_file):
    """
    A function that creates contig file followinng fasta format
    :Parameters:
          contigs : contig dictionary
          out_file : the name of the file to create
    """
    with open(out_file,"w") as f_out:
        for i in range(len(contigs)):
            f_out.write(">contig_{} len={}\n{}\n\n".format(i,contigs[i][1], fill(contigs[i][0])))

# III) Simplification du graphe de de Bruijn:
## a)Résolution des bulles
def std(list_values):
    """
    A function that gives the standard deviation of a list of values
    :Parameters:
          list_values : list of values
    Returns: the standard deviation value
    """
    return statistics.stdev(list_values)


def path_average_weight(graph, path):
    """
    A function that calculates the average weight of a given path in a graph
    :Parameters:
          graph : nx graph
          path : the path in the graph
    Returns: the average weight of the path
    """
    new_g = graph.subgraph(path)
    weight = new_g.degree(nbunch = new_g, weight = "weight")
    somme=0
    for i in weight:
        somme += i[1]
    mean_wei = somme/len(path)
    return mean_wei

def remove_paths(graph, list_paths, delete_entry_node, delete_sink_node):
    """
    A function that takes a graph, a list of paths and 2 boolean values to determinate if
    we keep or delete the entry and sink nodes and returns a cleaned graph
    :Parameters:
          graph : nx graph
          list_paths : list of paths
          delete_entry_node : boolean True/False
          delete_sink_node : boolean True/Flase
    Returns: a cleaned graph from unneeded paths
    """
    clean_graph = graph
    entry_node = 1
    sink_node = 1
    if delete_entry_node :
        entry_node -= 1
    if delete_sink_node :
        sink_node -= 1
    for path in list_paths:
        for i in range(entry_node, len(path)-sink_node):
            clean_graph.remove_node(path[i])
    return clean_graph


def select_best_path(graph, list_paths, lst_len_path, lst_mean_weight, delete_entry_node = False,
    delete_sink_node = False):
    """
    A function that takes a graph, a list of paths, a list of mean weights and 2 booleans
    and returns the best selected path of the graph
    we consider that the best path is :
        - highly frequented
        - has a big weight
    :Parameters:
          graph : nx graph
          lst_paths : list of paths
          lst_len_path : list of paths lenght
          lst_mean_weight : list of average weights
          delete_entry_node : boolean True/False
          delete_sink_node : boolean True/Flase
    Returns: the best selected path
    """
    max_weight = max(lst_mean_weight)
    ind_w = []
    for i, wei in enumerate(lst_mean_weight):
        if wei == max_weight:
            ind_w.append(i)
    if len(ind_w) > 1:
        max_len = max(lst_len_path)
        ind_l = []

        for i in ind_w:
            if lst_len_path[i] == max_len:
                ind_l.append(i)

        if len(ind_l) > 1:
            random.seed(9001)
            rnd_ind = randint(0,len(list_paths))
            best_path = ind_l[rnd_ind]
        else:
            best_path = ind_l[0]
    else:
        best_path = ind_w[0]
    bad_paths = list_paths[:best_path] + list_paths[best_path+1:]
    clean_graph = remove_paths(graph, bad_paths, delete_entry_node, delete_sink_node)
    return clean_graph

def solve_bubble(graph, ancestor_node, successor_node):
    """
    A function that takes a graph, an ancestor and a successor nodes and returns
    A graph cleaned from bubbles (it calls select_best_path to keep only one
    path from the bubble)
    :Parameters:
          graph : nx graph
          ancestor_node : an ancestor node
          successor_node : a successor node
    Returns: a graph cleaned from bubbles
    """
    all_paths = list(nx.algorithms.simple_paths.all_simple_paths(graph,ancestor_node,
        successor_node))
    graph_path_weights = []
    graph_path_lengths = []
    for i in all_paths:
        graph_path_weights.append(path_average_weight(graph, i))
        graph_path_lengths.append(len(i))
    graph_no_bull = select_best_path(graph, all_paths,graph_path_lengths, graph_path_weights)
    return graph_no_bull

def simplify_bubbles(graph):
    """
    A function that takes a graph and clean it from bubbles, it will find the bubbles
    and call solve_bubble to remove them
    :Parameters:
          graph : nx graph
    Returns: the graph without bubbles
    """
    list_nodes = graph.nodes()
    list_bubbles = []

    for node in list_nodes:
        node_predecessors = list(graph.predecessors(node))
        if len(node_predecessors) > 1:
            for i in range(len(node_predecessors)):
                for j in range(i,len(node_predecessors)):
                    lowest_predecessor = nx.lowest_common_ancestor(graph,node_predecessors[i],
                        node_predecessors[j])
                    if lowest_predecessor:
                        list_bubbles.append([lowest_predecessor, node])

    graph_no_bull = graph
    for bubble in list_bubbles:
        graph_no_bull = solve_bubble(graph_no_bull, bubble[0], bubble[1])

    return graph_no_bull

## b) Détection des pointes (tips)
def solve_entry_tips(graph, list_entre):
    """
    A function that take in arguments a graph and a list of entry nodes and
    return a graph without undesirable paths
        :Parameters:
            graph : nx graph
            list_entre : a list of entry nodes
        :Returns: 
            A graph with no entry tips
    """
    node_pred = []
    lst_path = []
    wei_path = []
    len_path = []
    for node in list_entre:
        for desc in nx.descendants(graph, node):
            pred = list(graph.predecessors(desc))
            if len(pred) > 1 and desc not in node_pred:
                for path in nx.all_simple_paths(graph, node, pred):
                    lst_path.append(path)
                    wei_path.append(len(path))
                    len_path.append(path_average_weight(graph, path))
        graph = select_best_path(graph, lst_path, len_path, wei_path, False, False)
    return graph

def solve_out_tips(graph, list_sink):
    """
    A funtction that take in arguments a graph un a list of out nodes and return
    a graph without undesirable paths
        :Parameters: 
            graph : nx graph 
            list_sink : a list of out nodes
        :Return:
        A graph with no out tips       
    """
    node_desc = []
    lst_path = []
    wei_path = []
    len_path = []
    for node in list_sink:
        for anc in nx.ancestors(graph, node):
            succ = list(graph.successors(anc))
            if len(succ) > 1 and anc not in node_desc:
                for path in nx.all_simple_paths(graph, node, anc):
                    lst_path.append(path)
                    wei_path.append(len(path))
                    len_path.append(path_average_weight(graph, path))
        if len(wei_path) > 0 :
            graph = select_best_path(graph, lst_path, len_path, wei_path, False, False)
    return graph

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    #args = get_arguments()

if __name__ == '__main__':
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("fasta_file", help="the fasta file", type=str)
    PARSER.add_argument("len_kmer", help="the length of the kmers", type=int)
    ARGS = PARSER.parse_args()
    FASTA_FILE = ARGS.fasta_file
    LEN_KMER = ARGS.len_kmer

    # 1. Reading the file and building the graph
    dic = build_kmer_dict(FASTA_FILE, LEN_KMER)
    print(dic)
    print("\n")
    G = build_graph(dic)
    #show_graph(G)

    starting_node = get_starting_nodes(G)
    print(starting_node)
    out_node = get_sink_nodes(G)
    print(out_node)



    # 2. Bubble resolution : remove all bubbles
    G = simplify_bubbles(G)

    # 3. Resolution of entry and out tips
    G = solve_entry_tips(G, starting_node)
    G = solve_out_tips(G, out_node)

    # 4. writting contig in fasta file:
    conti = get_contigs(G, starting_node, out_node)
    print(conti)
    save_contigs(conti,FASTA_FILE+".fna")
    
    







 
    
    
    
    