"""
test_setup_methods.py
J. Sumabat, N. Lai, S. Peck, Y. Huang, 3/12/2024

test_setup_methods.py contains the test cases for setup_methods.py.
"""

import unittest
import sys
sys.path.append('../lib/consensus_formation')
import setup_methods as sm
import os

em='../data/empty_folder/'
ghst = '../data/GHOST/'
wdir = '../working_directory/'

class TestCheckMakeDir(unittest.TestCase):
    
    def test_dir_existing(self):
        method_call = sm.check_make_dir(em)
        expected = str(os.path.abspath(em))+"/"
        self.assertEqual(method_call, expected)
        
    def test_dir_nonexisting(self):
        os.rmdir(ghst)
        method_call = sm.check_make_dir(ghst)
        expected = str(os.path.abspath(ghst))+"/"
        self.assertEqual(method_call, expected)
        

class TestSetupDirectorySystem(unittest.TestCase):
    
    def test_output_length(self):
        method_call = len(sm.setup_directory_system(wdir))
        expected = 6
        self.assertEqual(method_call, expected)
        

if __name__ == '__main__':
    unittest.main()