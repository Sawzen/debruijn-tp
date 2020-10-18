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
import sys , os
import pprint
import scipy
import networkx as nx
from networkx import algorithms
import matplotlib.pyplot as plt
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


# 1) CREATION DU GRAPHE de DE BRUIJN
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
    '''take a fastq file and return an iterator
    of a given sequence
    '''
    fastq_file = open(fastq)
    lines = iter(fastq_file.readlines())
    for line in lines :
        yield next(lines)
        next(lines)
        next(lines)


def cut_kmer(seq, kmer):
    ''' Fonction qui prend deux arugments 
            seq =  une liste sequence en str donnée par la fonction read_fastq
            kmer : une taille k et retourne un générateur de k_mer   '''
    seq = seq.strip('\n')
    for j in range(len(seq) - kmer + 1):
        yield seq[j:j+kmer]
        
        
def build_kmer_dict(fastq, km_len):
    """
    """
    dico = {}
    for seq in read_fastq(fastq):
        for k_mer in cut_kmer(seq, km_len):
            if k_mer not in dico:
                dico[k_mer] = 1
            else :
                dico[k_mer] += 1
    return dico
    

def build_graph(dico_kmer):
    """
    prendra en entrée un dictionnaire de k-mer et créera l’arbre de k-mers 
    préfixes et suffixes décrit précédemment. Les arcs auront pour paramètre 
    obligatoire un poids nommé “weight”
    """
    graph = nx.DiGraph()
    for k_mer in dico_kmer:
        node1 = k_mer[:-1]
        node2 = k_mer[1:]
        graph.add_edge(node1 , node2 , weight = dico_kmer[k_mer])
    return graph
    
def show_graph(graph):
    """
    """
    nx.draw(graph, pos=nx.spring_layout(graph))
    #pos = nx.spring_layout(graph)
    #nx.draw_networkx_labels(graph, pos)
    plt.draw()
    plt.show()

def get_starting_nodes(graph):
    '''take a graph and return a list of entry
    nodes
    '''
    list_entre = []
    for node in graph :
        pred = list(graph.predecessors(node))
        if (not pred) :
            #print("Pas de predecesseur\n")
            list_entre.append(node)
    return list_entre

def get_sink_nodes(graph):
    '''take a graph and return a list of entry
    nodes
    '''
    list_sink = []
    for node in graph :
        succ = list(graph.successors(node))
        if (not succ) :
            #print("Pas de predecesseur\n")
            list_sink.append(node)
    return list_sink

def get_contigs(graph, lst_start, lst_end):
    """
    prend un graphe, une liste de noeuds d’entrée et une liste de sortie et 
    retourne une liste de tuple(contig, taille du contig)
    """
    contigs = []
    for source in lst_start :
        for target in lst_end :
            if algorithms.has_path(graph, source, target) == True :
                path = algorithms.shortest_path(graph, source, target)
                contig = path[0]
                for i in range(len(path)-1):
                    contig += path[i+1][-1]
                contigs.append((contig, len(contig)))
    return contigs

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contigs, fichier_out):
    with open(fichier_out,"w") as f_out:
        for i in range(len(contigs)):
            f_out.write(">contig_{} len={}\n{}\n\n".format(i,contigs[i][1], fill(contigs[i][0])))

def std(list_values):
    return statistics.stdev(list_values)
    
# def path_average_weight(graph, contigs):
    graph_1path = graph.subgraph(contigs)
    # weight = 0
    # c = 0
    # for link in graph_1path.edges(data = True):
        # weight += link[2]["weight"]
        # c += 1
        # print(link)
    # mean_weight = weight/c  
    # return mean_weight
   
def path_average_weight(graph, path):
    """Take a graph and a path and return average weigth"""
    new_G = graph.subgraph(path)
    weight = new_G.degree(nbunch = new_G, weight = "weight")
    mean_wei = weight/len(path)
    return mean_wei
    
def remove_paths(graph, list_paths, delete_entry_node, delete_sink_node):
    clean_graph = graph 
    entry_node = 1
    sink_node = -2
    if delete_entry_node :
        entry_node = 0
    if delete_sink_node :
        entry_node = -1  
    for i in range(len(list_paths)): 
        clean_graph.remove_node_from(path[entry:sink])
    return clean_graph
    
def select_best_path(graph, list_paths, lst_len_path, lst_mean_weight, delete_entry_node = False, delete_sink_node = False):
    """
    qui prend un graphe, une liste de chemin, une liste donnant la longueur de chaque chemin, 
    une liste donnant le poids moyen de chaque chemin, delete_entry_node pour indiquer si les noeuds 
    d’entrée seront supprimés et delete_sink_node pour indiquer si les noeuds de sortie seront supprimés et 
    retourne un graphe nettoyé des chemins indésirables. Par défaut, delete_entry_node et delete_sink_node seront ici à False.
    Le meilleur chemin (liste de noeuds consécutif et acyclique) sera identifié par 3 critères:
    Un chemin est plus fréquent
    Un chemin est plus long
    Le hasard, vous imposerez une seed à 9001
    """
    max_weight = max(lst_mean_weight)
    ind_w = []
    for i, w in enumerate(lst_mean_weight):
        if w == max_weight:
            ind_w.append(i)
    if len(ind_w) > 1:
        max_len = max(path_lengths)
        ind_l = []
        
        for i, l in enumerate(lst_len_path):
            if l == count_max_len:
                ind_l.append(i)

        if len(ind_l) > 1:
            Random.seed(9001)
            best_path = random.choice[ind_l]
        else:
            best_path = ind_l[0]
    else:
        best_path = ind_w[0]
    bad_paths = list_paths[:best_path] + list_paths[best_path+1:]
    clean_graph = remove_paths(graph, bad_paths, delete_entry_node, delete_sink_node)
    return clean_graph    

def solve_bubble(graph, ancestor_node, successor_node):
    """
    """
    all_paths = list(nx.algorithms.simple_paths.all_simple_paths(graph, ancestor_node,successor_node))
    graph_path_weights = []
    graph_path_lengths = []
    for i in all_paths:
        graph_path_weights.append(path_average_weight(graph, i))
        graph_path_lengths.append(len(i))
    graph_no_bull = select_best_path(graph, all_paths,graph_path_lengths, graph_path_weights)
    return graph_no_bull

def simplify_bubbles(graph):
    """
    """
    list_nodes = graph.nodes()
    list_bubbles = []
    
    for node in list_nodes:
        node_predecessors = list(graph.predecessors(node))
        if len(node_predecessors) > 1:
            for i in range(len(node_predecessors)):
                for j in range(i,len(node_predecessors)):
                    lowest_predecessor = nx.lowest_common_ancestor(graph,node_predecessors[i], node_predecessors[j])
                    list_bubbles.append([lowest_common_ancestor, node])
                    
    graph_no_bull = graph                    
    for bubble in list_bubbles:
        graph_no_bull = solve_bubble(graph_no_bull, bubble[0], bubble[1])
    return graph_no_bull
    
#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

if __name__ == '__main__':
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("fasta_file", help="the fasta file", type=str)
    PARSER.add_argument("len_kmer", help="the length of the kmers", type=int)
    ARGS = PARSER.parse_args()
    FASTA_FILE = ARGS.fasta_file
    LEN_KMER = ARGS.len_kmer

    # ficher eva71_two_reads.fq
    dico = build_kmer_dict(FASTA_FILE, LEN_KMER)
    print(dico)
    print("\n")
    G = build_graph(dico)
    show_graph(G)
    starting_node = get_starting_nodes(G)
    print(starting_node)
    sink_node = get_sink_nodes(G)
    print(sink_node)

    print("\n")

    contig = get_contigs(G, starting_node, sink_node)
    print(contig)

    print("\n")
    save_contigs(contig,FASTA_FILE+".fna")
    
    print(contig[0][0])
    # print(nx.algorithms.simple_paths.all_simple_paths(G,"TCAGAGC", "AATTGTG"))
    # weight_G_1path = path_average_weight(G,nx.algorithms.simple_paths.all_simple_paths(G,"TCAGAGC", "AATTGTG")[0])
    # c = nx.path_graph(G,1)
    # weight_G_1path = path_average_weight(G, contig[0][0])
    # print(weight_G_1path)
    
    







 
    
    
    
    