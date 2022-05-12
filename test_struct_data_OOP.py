# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 13:15:17 2021

@author: melik, marie
"""

import unittest
from projet_bs2_OOP import * 
import random


lien = "toy_2comp.txt"

# To check the function "is_interaction_file"
empty_file = "empty_file.txt"
bad_file = "bad_file.txt"

interactome1 = Interactome(lien)
dict1 = interactome1.get_int_dict()


class Test_assembly(unittest.TestCase):


    def test_format_dict(self):   
        
        draw1 = random.choice(list(dict1.keys()))
        draw2 = random.choice(dict1[draw1])
        answer = draw1 in dict1[draw2]
        
        size_list = []
        for i in dict1:
            size_list.append(len(dict1[i]))
        list_verif = [5, 3, 3, 4, 3, 3, 1, 2, 3, 2, 3, 1, 1, 1, 1]
        
        # Check that the duplicate exists in the dictionary 
        self.assertEqual(answer, True)
        # Check that the result is an instance of list
        self.assertIsInstance(dict1, dict)
        # Check the length of the generated dict
        self.assertEqual(len(dict1), 15)
        # Check dictionary structure
        self.assertEqual(size_list, list_verif)
        
        
    
    def test_format_list(self):  
        list1 = interactome1.get_int_list()
        draw1 = random.choice(list1)
        answer = (draw1[1], draw1[0]) in list1

        # Check that the reverse duplicate does not exist in the list
        self.assertEqual(answer, False)
        # Check that the result is an instance of list
        self.assertIsInstance(list1, list)
        # Check the length of the generated list
        self.assertEqual(len(list1), 18)
    
    def test_format_mat(self):  
        mat1 = interactome1.mat
        def is_sym_mat(mat1):
            a = mat1.shape[0]
            for i in range(a):
                for j in range(a): 
                    if mat1[i][j] != mat1[j][i]:
                        return False
            return True
        answer = is_sym_mat(mat1)
        first_line_mat_l = [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        # Check that the matrix is ​​symmetrical
        self.assertEqual(answer, True)
        # Check that the size of the matrix is ​​equivalent to the size of the dictionary
        self.assertEqual(mat1.shape[0], len(dict1))
        # Check the first row of the matrix
        self.assertEqual(list(mat1[0]), first_line_mat_l)
        
        
    def test_format_interactionfile(self):
        interactome_empty = Interactome(empty_file)
        interactome_bad = Interactome(bad_file)
        
        self.assertEqual(interactome_empty.is_interaction_file(), False)
        self.assertEqual(interactome_bad.is_interaction_file(), False)
        self.assertEqual(interactome1.is_interaction_file(), True)
        
    
    def test_count_edge(self):  
        self.assertEqual(interactome1.count_edges(), 18)

        
    def test_count_vertices(self):   
        self.assertEqual(interactome1.count_vertices(), 15)
    
    
    def test_get_degree(self):
        self.assertEqual(interactome1.get_degree("B"), 3)
        self.assertEqual(interactome1.get_degree("D"), 3)
        self.assertEqual(interactome1.get_degree("F"), 1)
        self.assertEqual(interactome1.get_degree("H"), 2)
        self.assertEqual(interactome1.get_degree("J"), 2)
    
    
    def test_get_max_degree(self): # check that the result is a string
        self.assertEqual(interactome1.get_max_degree()[0], "A")
        self.assertEqual(interactome1.get_max_degree()[1], 5)
    
    
    def test_get_max_degree1(self): # check that the maximum degree of the protein corresponds well
        result = interactome1.get_max_degree()
        self.assertEqual(interactome1.get_degree(result[0]), result[1])
    
    def test_get_ave_degree(self):  #  check that the result with the example file is 2.4
        result = interactome1.get_ave_degree()
        self.assertEqual(result, 2.4)
    
    def test_count_degree(self):   #  check the result
        self.assertEqual(interactome1.count_degree(1), 5)
        self.assertEqual(interactome1.count_degree(2), 2)
        self.assertEqual(interactome1.count_degree(3), 6)
        self.assertEqual(interactome1.count_degree(4), 1)
        
    def test_density(self):
        self.assertEqual(round(interactome1.density(), 3), 0.171)
        
    def test_clustering(self):
        self.assertEqual(interactome1.clustering("A"), 0.5)
        self.assertEqual(interactome1.clustering("B"), 1)
        self.assertEqual(round(interactome1.clustering("C"),2), 0.67)
        self.assertEqual(round(interactome1.clustering("D"),2), 0.67)
        self.assertEqual(interactome1.clustering("F"), 0)
        self.assertEqual(interactome1.clustering("H"), 1)
        
    def test_count_CC(self):
        self.assertEqual(interactome1.count_CC()[0], 3)
        
        dict2 = interactome1.count_CC()[1]
        size_comp = list(dict2.values())
        size_comp.sort(reverse = True)
        self.assertEqual(interactome1.count_CC()[1][0], size_comp[0])
        self.assertEqual(interactome1.count_CC()[1][1], size_comp[1])
        self.assertEqual(interactome1.count_CC()[1][2], size_comp[2])
        

    def test_extract_CC(self):

        comp1_l = ['B', 'A', 'G', 'C', 'D', 'E', 'F']
        comp2_l = ['I', 'H', 'J', 'K', 'L', 'M']
        comp3_l = ['Z', 'Y']
        
        self.assertEqual(interactome1.extract_CC("A"), comp1_l)
        self.assertEqual(interactome1.extract_CC("H"), comp2_l)
        self.assertEqual(interactome1.extract_CC("Z"), comp3_l)
     
        
    def test_compute_CC(self):
        
        result_l = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2]
        self.assertEqual(interactome1.compute_CC(), result_l)
        
  
        
    
unittest.main()
