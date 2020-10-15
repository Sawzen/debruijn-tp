# Importation des modules

import argparse
import sys , os
import pprint
import networkx as nx
from networkx import algorithms

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
    
    
    
    
    