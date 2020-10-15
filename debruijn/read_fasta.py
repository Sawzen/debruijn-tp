# Importation des modules

import argparse
import sys , os
import pprint
import networkx as nx
from networkx import algorithms
import matplotlib.pyplot as plt

# 1) CREATION DU GRAPHE de DE BRUIJN
##  a) Identification des kmer unique

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
    kmer_tree = nx.DiGraph()
    for k_mer in dico_kmer:
        node1 = k_mer[:-1]
        node2 = k_mer[1:]
        kmer_tree.add_edge(node1 , node2 , weight = dico_kmer[k_mer])
    return kmer_tree
    
def show_graph(graph):
    """
    """
    nx.draw(graph, pos=nx.spring_layout(graph))
    #pos = nx.spring_layout(graph)
    #nx.draw_networkx_labels(graph, pos)
    plt.draw()
    plt.show()  
    
#### Main
dico = build_kmer_dict("../data/eva71_two_reads.fq", 8)
print(dico)
G = build_graph(dico)
show_graph(G)


 
    
    
    
    