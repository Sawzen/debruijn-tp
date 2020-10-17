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
    graph_clean = graph 
    entry_node = 1
    sink_node = -2
    if delete_entry_node :
        entry_node = 0
    if delete_sink_node :
        entry_node = -1  
    for i in range(len(list_paths)): 
        graph_clean.remove_node_from(path[entry:sink])
    return graph_clean
    
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
            best = random.choice[ind_l]
        else:
            best = ind_l[0]
    else:
        best = ind_w[0]
    paths.pop(best)
    return remove_paths(graph, list_paths, delete_entry_node, delete_sink_node)

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
    
    

