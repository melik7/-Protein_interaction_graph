import numpy as np
import random


### Chapter 1 ###

def read_interaction_file_dict(filename):

    '''
    Returns a dictionnary of all interactions for an undirected interaction graph

            Parameters:
                    filename : a file containing an undirected interaction graph 

            Returns:
                    dico (dict): a dictionnary with each proteins (keys) and proteins that interact with
    '''

    with open (filename, "r") as fh:
        lines = fh.readlines()
        dico = {}
        for line in lines[1:]: 
            info_ligne = line.split()
            vertex = info_ligne[0]
            neighbor = info_ligne[1]
            # Create key with the first value if it doesn't exist
            if vertex not in dico:   
                dico[vertex] = [neighbor]
            # Check that there are no duplicates before adding it 
            if neighbor not in dico [vertex]:  
                dico[vertex].append(neighbor) 
            
            # Create key for interactions undirected and not included in the file 
            if neighbor not in dico:  
                dico[neighbor] = [vertex]
            #  Check that there are no duplicates before adding it
            if vertex not in dico[neighbor]:   
                dico[neighbor].append(vertex)
        
        return dico


def read_interaction_file_list(filename):

    '''
    Returns a list of pairs of proteins interactions for an undirected interaction graph 

            Parameters:
                    filename : a file containing an undirected interaction graph 

            Returns:
                    interaction_list (list): a list of all interactions pairs  
    '''
    
    with open (filename, "r") as fh:
        lines = fh.readlines()
        interaction_list = []
        for line in lines[1:]:
            info_ligne = line.split()
            vertex = info_ligne[0]
            neighbor = info_ligne[1]
            # Check that the inverse couple doesn't exist to avoid duplicates of interactions pairs 
            if (neighbor,vertex) not in interaction_list and neighbor != vertex: 
                interaction_list.append((vertex, neighbor))

        return interaction_list


def read_interaction_file_mat(filename):

    '''
    Returns adjacency matrix representing an undirected interaction graph and a ordered list of vertices 

            Parameters:
                    filename : a file containing an undirected interaction graph 

            Returns:
                    mat (matrix): symetric matrix with 0 for non interaction and 1 for interaction
                    keys_dico (list): the ordered list of vertices
    '''

    dico = read_interaction_file_dict(filename)
    keys_dico = list(dico.keys())
    mat_size = len(keys_dico)
    mat = np.zeros((mat_size, mat_size))
    for vertex in dico:
        for neighbor in dico[vertex]:
            mat[keys_dico.index(vertex), keys_dico.index(neighbor)] = 1
              
    return mat, keys_dico


def read_interaction_file(filename):

    '''
    Returns a dictionnary, a list, a matrix and a ordered list of vertices of all interactions for an undirected interaction graph 

            Parameters:
                    filename : a file containing an undirected interaction graph 

            Returns:
                    d_int (dict): a dictionnary with each proteins (keys) and proteins that interact with
                    l_int (list): a list of all interactions pairs
                    m_int (matrix): symetric matrix with 0 for non interaction and 1 for interaction
                    l_som (list): the ordered list of vertices
    '''

    d_int = read_interaction_file_dict(filename)
    l_int = read_interaction_file_list(filename)
    m_int, l_som = read_interaction_file_mat(filename)

    return d_int, l_int, m_int, l_som


#### QUESTION STRUCTURE 5 #######

# Demander en argument de la fonction read_interaction_file ce que l'utilisateur veut récupérer
# Cela évitera de générer toutes les structures de données pour rien. 
# Ou sinon on pourrait sortir le dictionnaire des fonctions et en utiliser qu'un seul
# Cela permettrait de le générer qu'une seule fois et ca éviterait également d'aller
# lire plusieurs fois le fichier d'interactions.

def is_interaction_file(filename):

    '''
    Returns False if the file is empty, if it doesn't have the first line counting the number of interactions, 
    if it isn't the right number of interactions or if it doesn't have the right number of columns

            Parameters:
                    filename : a file containing an undirected interaction graph 

            Returns:
                    boolean : True or False 
    '''

    with open (filename, "r") as fh:
        lines = fh.readlines()
        if len(lines) ==  0:
            return False
        try:
            first_line = int(lines[0].split()[0])
        except:
            return False 
        else: 
            if first_line != len(lines)-1:
                return False 
        for line in lines[1:]: 
            char_numb = line.split()
            if len(char_numb) != 2:
                return False 

        return True 



### Chapter 2 ###

def count_vertices(file):

    '''
    Returns the number of vertices of an undirected graph.

    Parameters:
            file : a file containing an undirected interaction graph

    Returns: 
            integer : the number of vertices (which corresponds to the number 
            of the keys in the interaction dictionnnary)
    '''

    dico_int = read_interaction_file_dict(file)

    return len(dico_int.keys())


def count_edges(file):

    '''
    Returns the number of the edges of an undirected file. 

        Parameters:
                file : a file containing an undirected interaction graph

        Returns: 
                integer : the number of edges (which corresponds to the length of 
                the interaction list)
    '''

    list_int = read_interaction_file_list(file)

    return len(list_int)


def clean_interactome(filein, fileout):

    '''
    Returns a file containing an undirected interaction graph without redundant 
    interactions and homo-dimers.

        Parameters:
                filein : a file containing an undirected interaction graph
                fileout : a file which will contains the new undirected interaction
                graph

        Returns: 
                fileout : file containing the new graph 
    '''

    list_int_clean = read_interaction_file_list(filein)
    file_int_clean = open(fileout, "w")
    file_int_clean.write(str(len(list_int_clean)) + "\n")
    for prot in list_int_clean:
        file_int_clean.write(str(prot[0]) + str("\t") + str(prot[1]) + "\n")
    file_int_clean.close()

    print("The file has been generated")


def degree_prot(file):

    '''
    Returns a dictionnary where the keys correspond to the vertices and the values 
    their degree. 

        Parameters:
                file : a file containing an undirected interaction graph

        Returns: 
                dico_degree : the dictionnary containing the vertices and their degree
    '''

    dico_int = read_interaction_file_dict(file)
    dico_degree = {}
    for key in dico_int:
        dico_degree[key] = dico_degree.get(key, len(dico_int[key]))

    return dico_degree


def get_degree(file, prot):

    '''
    Returns the degree of a protein of interest and its name.

        Parameters:
                file : a file containing an undirected interaction graph
                prot (str) : name of the protein of interest

        Returns: 
                tuple : name of the protein of interest and its degree
    '''

    dico_degree = degree_prot(file)

    return dico_degree[prot]


def get_max_degree(file):

    '''
    Returns the name of the maximum degree protein and its degree.

        Parameters:
                file : a file containing an undirected interaction graph

        Returns: 
                tuple : name of the maximum degree protein and its degree
    '''

    dico_degree = degree_prot(file)
    k = list(dico_degree.keys())
    v = list(dico_degree.values())

    return k[v.index(max(v))], v[v.index(max(v))]


def get_ave_degree(file):

    '''
    Returns the average degree of the proteins in the interaction graph.

        Parameters:
                file : a file containing an undirected interaction graph

        Returns: 
                mean (float) : mean degree
    '''

    dico_degree = degree_prot(file)
    list_values = list(dico_degree.values())
    mean = round(sum(list_values) / len(list_values), 2)

    return mean 


def count_degree(file, deg):

    '''
    Returns the number of proteins in the interaction graph whose degree is exactly 
    equal to the degree of interest.

        Parameters:
                file : a file containing an undirected interaction graph
                deg (int) : degree of interest

        Returns:
                nber_prot (int) : number of proteins whose degree is equal to deg     
    '''

    dico_degree = degree_prot(file)
    nber_prot = sum(dico_degree[key] == deg for key in dico_degree)

    return nber_prot


def histogram_degree(file, dmin, dmax):

    '''
    Returns the number of the proteins having a degree equal to d, where d is a 
    degree between a minumum (dmin) and a maximum (dmax) degree. The result is 
    dispalyed as a histogram.

        Parameters:
                file : a file containing an undirected interaction graph
                dmin (int) : minimum degree
                dmax (int)  : maximum degree

        Returns:
                the histogram         
    '''

    list_deg = [i for i in range(dmin, dmax+1)]
    list_nber_prot = []
    for deg in list_deg:
        list_nber_prot.append(count_degree(file, deg))
    for i in range(0, len(list_nber_prot)):
        print(str(list_deg[i]) + " " + list_nber_prot[i] * '*')


# Nous pouvons constater que le nombre de protéines décroit à mesure que le degré augmente. 
# Ainsi, il semblerait que cette distribution corresponde une loi de puissance se caractérisant
# par la diminution de la fréquence d'un évènement (ici le nombre de protéines), au fur et 
# à mesure que la taille d'un autre évènement augmente (ici le degré des protéines). 
# La faible proportion de protéines ayant un fort degré et donc de nombreux voisins, suggère
# la présénce de "hubs" protéiques au sein du réseau d'interactions.   

