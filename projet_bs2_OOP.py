import numpy as np
import random
import copy

### Chapter 3 ###


#########  QUESTION 3.1.1  ##########
# On lit le fichier d'intercaction autant de fois que l'on parcoure de degrès
# soit dmax + 1 - dmin


class Interactome:
    def __init__(self, filename):
        
        
        self.int_list = []
        self.int_dict = {}
        with open (filename, "r") as fh:
            lines = fh.readlines()
            for line in lines[1:]:
                
        
                #Dictionnary 
                info_ligne = line.split()
                vertex = info_ligne[0]
                neighbor = info_ligne[1]
                # Create key with the first value if it doesn't exist
                if vertex not in self.int_dict:   
                    self.int_dict[vertex] = list(neighbor)
                # Check that there are no duplicates before adding it 
                if neighbor not in self.int_dict [vertex]:  
                    self.int_dict[vertex].append(neighbor) 
                
                # Create key for interactions undirected and not included in the file 
                if neighbor not in self.int_dict:  
                    self.int_dict[neighbor] = list(vertex)
                #  Check that there are no duplicates before adding it
                if vertex not in self.int_dict[neighbor]:   
                    self.int_dict[neighbor].append(vertex)
            
                #Interaction list                                
                info_ligne = line.split()
                vertex = info_ligne[0]
                neighbor = info_ligne[1]
                if (neighbor,vertex) not in self.int_list and (vertex,neighboor) not in self.int_list and neighbor != vertex:  
                    self.int_list.append((vertex, neighbor))
        
        self.proteins = list(self.int_dict.keys())
        self.filename = filename 
    
        keys_dico = list(self.int_dict.keys())
        mat_size = len(keys_dico)
        self.mat = np.zeros((mat_size, mat_size))
        for vertex in self.int_dict:
            for neighbor in self.int_dict[vertex]:
                self.mat[keys_dico.index(vertex), keys_dico.index(neighbor)] = 1
    
    # Define accessors
    def get_int_dict(self):
        return self.int_dict
    def get_int_list(self):
        return self.int_list
    def get_proteins(self):
        return self.proteins


    # Define mutators
    def set_int_dict(self, new_dict):
        self.int_dict = new_dict
    def set_int_list(self, new_list):
        self.int_list = new_list
    def set_proteins(self, new_prot):
        self.proteins = new_prot
    

    def is_interaction_file(self):

        '''
        Returns False if the file is empty, if it doesn't have the first line counting the number of interactions, 
        if it isn't the right number of interactions or if it doesn't have the right number of colomn
    
                Parameters:
                        filename : a file containing an undirected interaction graph 
    
                Returns:
                        boolean : True or False 
        '''
    
        with open (self.filename, "r") as fh:
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
    
    
    def count_vertices(self):
    
        '''
        Returns the number of vertices of an undirected graph.
    
        Parameters:
    
        Returns: 
                integer : the number of vertices (which corresponds to the number 
                of the keys in the interaction dictionnnary)
        '''
    
        return len(self.get_int_dict().keys())
    
    
    def count_edges(self):
    
        '''
        Returns the number of the edges of an undirected file. 
    
            Parameters:
    
            Returns: 
                    integer : the number of edges (which corresponds to the length of 
                    the interaction list)
        '''
    
        return len(self.get_int_list())
    
    
    def clean_interactome(self, fileout):
    
        '''
        Returns a file containing an undirected interaction graph without redundant 
        interactions and homo-dimers.
    
            Parameters:
                    fileout : a file which will contains the new undirected interaction graph
    
            Returns: 
                    fileout : file containing the new graph 
        '''
    
        file_int_clean = open(fileout, "w")
        file_int_clean.write(str(len(l_int_clean)) + "\n") ### ?????
        for prot in self.get_int_list():
            file_int_clean.write(str(prot[0]) + str("\t") + str(prot[1]) + "\n")
        file_int_clean.close()

        print("The file has been generated")
    

    def degree_prot(self):
    
        '''
        Returns a dictionnary where the keys correspond to the vertices and the values 
        their degree. 
    
            Parameters:
    
            Returns: 
                    dict_degree : the dictionnary containing the vertices and their degree
        '''
    
        dict_degree = {}
        for key in self.get_int_dict():
            dict_degree[key] = dict_degree.get(key, len(self.get_int_dict()[key]))

        return dict_degree
    
    
    def get_degree(self, prot):
    
        '''
        Returns the degree of a protein of interest and its name.
    
            Parameters:
                    prot (str) : name of the protein of interest
    
            Returns: 
                    tuple : name of the protein of interest and its degree
        '''
    
        dict_degree = self.degree_prot()

        return dict_degree[prot]
    
   
    def get_max_degree(self):
    
        '''
        Returns the name of the maximum degree protein and its degree.
    
            Parameters:
    
            Returns: 
                    tuple : name of the maximum degree protein and its degree
        '''
    
        dict_degree = self.degree_prot()
        k = list(dict_degree.keys())
        v = list(dict_degree.values())

        return k[v.index(max(v))], v[v.index(max(v))]  

    
    def get_ave_degree(self):
    
        '''
        Returns the average degree of the proteins in the interaction graph.
    
            Parameters:
    
            Returns: 
                    mean (float) : mean degree
        '''
    
        dict_degree = self.degree_prot()
        list_values = list(dict_degree.values())
        mean = round(sum(list_values) / len(list_values), 2)

        return mean 
    
    
    def count_degree(self, deg):
    
        '''
        Returns the number of proteins in the interaction graph whose degree is exactly 
        equal to the degree of interest.
    
            Parameters:
                    deg (int) : degree of interest
    
            Returns:
                    nber_prot (int) : number of proteins whose degree is equal to deg     
        '''
    
        dict_degree = self.degree_prot()
        nber_prot = sum(dict_degree[key] == deg for key in dict_degree)

        return nber_prot
    
    
    def histogram_degree(self, dmin, dmax):
    
        '''
        Returns the number of the proteins having a degree equal to d, where d is a 
        degree between a minumum (dmin) and a maximum (dmax) degree. The result is 
        dispalyed as a histogram.
    
            Parameters:
                    dmin (int) : minimum degree
                    dmax (int)  : maximum degree
    
            Returns:
                    the histogram         
        '''
    
        list_deg = [i for i in range(dmin, dmax+1)]
        list_nber_prot = []
        for deg in list_deg:
            list_nber_prot.append(self.count_degree(deg))
        for i in range(0, len(list_nber_prot)):
            print(str(list_deg[i]) + " " + list_nber_prot[i] * '*')
    

    def density(self):

        '''
        Returns the graph's density which is the retio between the number of edges 
        present and the total number of possible edges
    
            Parameters:

            Returns:
                    int : graph's density       
        '''

        count_edge = len(self.get_int_list())
        count_vertex = len(self.get_int_dict())
        count_tot_vertex = (count_vertex*(count_vertex-1))/2

        return count_edge/count_tot_vertex


    def clustering(self, prot):

        '''
        Returns the clustering coefficient of a protein. This coefficient corresponds
        to the propability to have two connected nodes knowing they have one neighbor 
        in common 
    
            Parameters:
                    prot (str) : name of the protein of interest

            Returns:
                    int : clustering coefficient of a protein       
        '''

        neighbors = self.get_int_dict()[prot]
        count_vertex = len(neighbors)
        neighbors2 = copy.deepcopy(neighbors)
        count = 0 
        # All neighbors of prot
        for protein in neighbors:  
            del neighbors2[0]
            # All neighbors of prot except neighbors already done   
            for i in neighbors2:  
                # Neighbors's neighbors of prot 
                for k in self.get_int_dict()[protein]:
                    if i == k: 
                        count+=1
        
        count_tot_vertex = (count_vertex*(count_vertex-1))/2

        return count/count_tot_vertex if count_tot_vertex != 0 else 0


    def dict_CC(self):

        '''
        Returns a dictionnary containing all the connected components of the graph. 
        The keys correspond to a numbering of the connected components, so the first 
        connected component get the number 0, the second get the number 1 etc... The values
        correspond to the list of the proteins in the connected component.  
    
            Parameters:

            Returns:
                    relat_comp_dict (dict) : connected component dictionary        
        '''
        
        dico = copy.deepcopy(self.get_int_dict()) 
        # First element from which the component will be built.
        random.seed(1)
        start = random.choice(list(dico))
        # List of all elements added to the component. This list avoid to neighbors to be recounted.  
        relat_comp_list = [start] 
        relat_comp_d = {}
        # Component number which will be incremented.
        num_comp = 0
        
        # Temporary list of neighbors line already added to the component and which allow to extend the component. 
        temporary_list1 = [start]
        # Temporary list of news neighbors of tempory_list
        temporary_list2 = []
        
        while True:
            # Line of neighbors
            for i in temporary_list1: 
                for j in dico[i]:
                    if j not in relat_comp_list and j not in temporary_list2:
                        temporary_list2.append(j)
            if len(temporary_list2) == 0:
                # If temporary_list2 = 0, it's mean that there are no more neighbors to add. 
                # So the component is added to the final dictionary 
                relat_comp_dict[num_comp] = relat_comp_list
                num_comp += 1
            
                # We remove the elements of relat_comp_list
                for i in relat_comp_list:
                    del dico[i]
                
                if len(dico) != 0:
                    # If there are still elements in the dictionary, a new element is selected randomly, 
                    # to build a new connected component. 
                    relat_comp_list = []
                    start = random.choice(list(dico))
                    temporary_list2 = [start]
                else:
                    break   
                
            relat_comp_list += temporary_list2
            # We go to the next line of neighbors 
            temporary_list1 = temporary_list2 
            temporary_list2 = []
     
        return relat_comp_dict
    

    def count_CC(self):

        '''
        Returns the number of connected components and the number of proteins for each of them.  

        A new dictionary is created : dictio, its keys correspond to the numbering of the connected components,
        and its values to the number of proteins or the len of the compenent. 
        The total number of connected components corresponds to the len of the dictionary.
        
            Parameters:

            Returns:
                    dictio (dict) : connected components dictionary  
                            
        '''

        dict_ref = copy.deepcopy(self.dict_CC())
        dictio = {}
        for component in dict_ref:
            dictio[component] = len(dict_ref[component])

        return len(dictio), dictio 
    

    def write_CC(self, fileout):

        '''
        Returns a file with all connected components.  

        The file is written as before, the number of proteins for each component 
        and then a list of all protein couple. 
        
            Parameters:

            Returns:
                    file : connected components file  
                            
        '''

        dictio = copy.deepcopy(self.dict_CC())
        file_int_clean = open(fileout, "w")
        for comp in dictio:
            file_int_clean.write(str(len(dictio[comp])) + " \t" + str(dictio[comp]) + "\n")
        file_int_clean.close()

        print("The file has been generated")
    
    
    def extract_CC(self, prot):
        
        '''
        Returns all vertex of the connected component of a protein.   
        
            Parameters:
                    prot (str) : name of the protein of interest

            Returns:
                    list : all vertex  
                            
        '''

        dictio = copy.deepcopy(self.dict_CC())
        for comp in dictio:
            if prot in dictio[comp]:

                return dictio[comp]

        # return None
    
                
    def compute_CC(self):

        '''
        Returns a list (lcc) of the connected component number in which each protein 
        in our protein list is found.  

        The break is prevents the dictionary from being browsed if the protein has been found.
        
            Parameters:

            Returns:
                    lcc (list) : connected component number  
                            
        '''

        dictio = copy.deepcopy(self.dict_CC())
        lcc = []
        for prot in self.get_proteins():
            a = 0
            for comp in dictio:
                if prot in dictio[comp]:
                    lcc.append(comp)
                    a += 1
                    break
            if a == 0:    
                lcc.append(None)
        return lcc
        
