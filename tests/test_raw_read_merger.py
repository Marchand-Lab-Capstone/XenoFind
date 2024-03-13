"""
test_raw_read_merger.py
J. Sumabat, N. Lai, S. Peck, Y. Huang, 3/12/2024

test_raw_read_merger.py contains testing classes and cases for the methods within
raw_read_merger.py.
"""

import unittest
import sys
sys.path.append('../lib/consensus_formation')
import raw_read_merger as rrm
import os

f5='../data/fast5_folder/'
p5='../data/pod5_folder/'
em='../data/empty_folder/'
mix='../data/mixed_folder/'
tar='../data/working_directory/merged_pod5/'
fname='merged'

class TestMergeReadsCommand(unittest.TestCase):

    def test_fast5_win(self):
        method_call = rrm.merge_reads_command(f5,
                                              'fast5',
                                               tar,
                                               fname)
        expected = "pod5 convert fast5 --recursive --force-overwrite {} -o {}{}.pod5".format(f5, tar, fname)
        self.assertEqual(method_call, expected)

        
    def test_pod5_win(self):
        method_call = rrm.merge_reads_command(p5,
                                              'pd5',
                                               tar,
                                               fname)
        expected = "pod5 --force-overwrite --recursive {} -o {}{}.pod5".format(p5, tar, fname)
        self.assertEqual(method_call, expected)  
    

class TestValidateReadDirectory(unittest.TestCase):
    def test_fast5_only(self):
        self.assertEqual(rrm.validate_read_directory(f5), 'fast5')
        
    def test_pod5_only(self):
        self.assertEqual(rrm.validate_read_directory(p5), 'pod5')
    
    def test_empty(self):
        self.assertEqual(rrm.validate_read_directory(em), '')
        
    def test_mixed(self):
        self.assertEqual(rrm.validate_read_directory(mix), '')
        

class TestValidateTargetDirectory(unittest.TestCase):
    def test_valid_file(self):
        self.assertEqual(rrm.validate_target_directory(em), True)
    
    def test_invalid_file(self):
        self.assertEqual(rrm.validate_target_directory("aslkdjf"), False)
        

class TestGenerateMergedPod5(unittest.TestCase):
    def test_fast5_only(self):
        self.assertEqual(str(rrm.generate_merged_pod5(f5, tar, fname)), "(0, 'fast5', True)")
    
    def test_pod5_only(self):
        self.assertEqual(str(rrm.generate_merged_pod5(p5, tar, fname)), "(0, 'pod5', True)")
        
    def test_mixed(self):
        self.assertEqual(str(rrm.generate_merged_pod5('mix', tar, fname)), "(1, '', True)")
    
    def test_bad_output(self):
        self.assertEqual(str(rrm.generate_merged_pod5(f5, 'asdfef', fname)), "(1, 'fast5', False)")
        
        
if __name__ == '__main__':
    unittest.main()