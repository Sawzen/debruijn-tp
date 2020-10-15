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

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Hager Elharty - Hocine Meraouna"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Elharty-Meraouna"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Elharty-Meraouna"
__email__ = "hocine.Meraouna@gmail.com"
__status__ = "Developpement"

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
	"""
	"""


def cut_kmer(seq, km_len):
	"""
	"""


def build_kmer_dict(fastq, km_len):
	"""
    prend un fichier fastq, une taille k- mer et retourne un dictionnaire 
    ayant pour clé le k-mer et pour valeur le nombre d’occurrence de ce k-mer
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
    main()
