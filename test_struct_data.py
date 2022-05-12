# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 13:15:17 2021

@author: melik, marie 
"""

import unittest
from projet_bs2 import * 
import random
import io
import sys

lien = "toy_2comp.txt"

# To check the function "is_interaction_file"
empty_file = "empty_file.txt"
bad_file = "bad_file.txt"


dict1 = read_interaction_file_dict(lien)


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
        list1 = read_interaction_file_list(lien)
        draw1 = random.choice(list1)
        answer = (draw1[1], draw1[0]) in list1

        # Check that the reverse duplicate does not exist in the list
        self.assertEqual(answer, False)
        # Check that the result is an instance of list
        self.assertIsInstance(list1, list)
        # Check the length of the generated list
        self.assertEqual(len(list1), 18)
    
    def test_format_mat(self):  
        mat1 = read_interaction_file_mat(lien)[0]
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
        
        self.assertEqual(is_interaction_file(empty_file), False)
        self.assertEqual(is_interaction_file(bad_file), False)
        self.assertEqual(is_interaction_file(lien), True)
        
    
    def test_count_edge(self): 
        result = count_edges(lien)
        self.assertEqual(result, 18)

        
    def test_count_vertices(self):    
        result = count_vertices(lien)
        self.assertEqual(result, 15)
    
    
    def test_get_degree(self):
        self.assertEqual(get_degree(lien, "B"), 3)
        self.assertEqual(get_degree(lien, "D"), 3)
        self.assertEqual(get_degree(lien, "F"), 1)
        self.assertEqual(get_degree(lien, "H"), 2)
        self.assertEqual(get_degree(lien, "J"), 2)
    
    
    def test_get_max_degree(self): # check that the result is a string
        self.assertEqual(get_max_degree(lien)[0], "A")
        self.assertEqual(get_max_degree(lien)[1], 5)
    
    
    def test_get_max_degree1(self): # check that the maximum degree of the protein corresponds well
        result = get_max_degree(lien)
        self.assertEqual(get_degree(lien, result[0]), result[1])
    
    def test_get_ave_degree(self):  #  check that the result with the example file is 2.4
        result = get_ave_degree(lien)
        self.assertEqual(result, 2.4)
    
    def test_count_degree(self):   #  check the result 
        self.assertEqual(count_degree(lien, 1), 5)
        self.assertEqual(count_degree(lien, 2), 2)
        self.assertEqual(count_degree(lien, 3), 6)
        self.assertEqual(count_degree(lien, 4), 1)
        
    def test_hist_degree(self):
        def test_hist():
            capturedOutput = io.StringIO()                  
            sys.stdout = capturedOutput                     
            histogram_degree(lien, 1,1)                                    
            sys.stdout = sys.__stdout__
            return capturedOutput.getvalue()
        self.assertEqual(test_hist()[:-1], "1 *****")
        

        
 
  
        
    
unittest.main()
