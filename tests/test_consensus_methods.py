"""
test_consensus_methods.py
J. Sumabat, N. Lai, S. Peck, Y. Huang, 3/12/2024

test_consensus_methods.py contains testing classes and cases for methods within consensus_methods.py
"""

import unittest
import sys
sys.path.append('../lib/consensus_formation')
import consensus_methods as cm
import os
import pysam


class TestBasecallCommand(unittest.TestCase):
    
    def test_structure(self):
        # Tests the output structure of the method. 
        method_call = cm.basecall_command('dorado', 'p5path/pod5.pod5', 'outpath/', 'outname') 
        expected = 'dorado basecaller hac --no-trim --emit-fastq p5path/pod5.pod5 > outpath/outname.fq'
        self.assertEqual(method_call, expected)
        
    # Further testing requires error handling within the method.
    # assumes dorado is in PATH when passing dorado


class TestMapToReference(unittest.TestCase):
    
    def test_structure(self):
        # Tests the output structure of the method. 
        method_call = cm.map_to_reference('minimap2', 'ref_path/ref.fasta', 'basecall_path/basecall.sam', 'outpath/', 'outname' )
        expected = 'minimap2 -ax map-ont --score-N 0 --MD --min-dp-score 10 ref_path/ref.fasta basecall_path/basecall.sam > outpath/outname.sam'
        self.assertEqual(method_call, expected)

     # as with basecall, requires minimap2 to be in path as well, and needs more error handling


class TestVsearchCommand(unittest.TestCase):
    
    def test_structure(self):
        # Tests the output structure of the method. 
        method_call = cm.vsearch_command('vsearch', 'fasta_path/fasta.fasta', 'outpath/', 'outname', .9)
        expected = "vsearch --cluster_fast fasta_path/fasta.fasta --id 0.9 --clusterout_sort --consout outpath/outname.fasta"
        self.assertEqual(method_call, expected)

    # Assumes files do actually exist in order to function, does not have error handling within method. 


class TestMedakaConsensusCommand(unittest.TestCase):
    
    def test_structure(self):
        # Tests the output structure of the method. 
        method_call = cm.medaka_consensus_command('medaka', 'fasta_path/trimmed_fasta.fasta', 'fasta_path/filtered_fasta.fasta', 'outpath/')
        expected = "medaka -i fasta_path/trimmed_fasta.fasta -d fasta_path/filtered_fasta.fasta -o outpath/ -m r1041_e82_400bps_hac_v4.2.0 -f -b 300"
        self.assertEqual(method_call, expected)


    # As above, so below. No error handling in here, assumes files are extant.


class TestReadTrim(unittest.TestCase):
    
    def test_output_type(self):
        # Tests the output to see if its a list.
        method_call = type(cm.read_trim('../data/working_dir/basecall_directory/bc_aligned.sam'))
        expected = type(["a", "a", "a"])
        self.assertEqual(method_call, expected)
        
    def test_sub_type(self):
        # Tests the output to see if the list contains a string.
        method_call = type(cm.read_trim('../data/working_dir/basecall_directory/bc_aligned.sam')[0])
        expected = type("this is a string")
        self.assertEqual(method_call, expected)
        

class TestSortFasta(unittest.TestCase):
    
    def test_output_type(self):
        # Tests the output to see if its a list
        method_call = type(cm.sort_fasta('../data/working_dir/fasta_directory/trimmed.fasta'))
        expected = type(["liiiiist"])
        self.assertEqual(method_call, expected)
        
    def test_sub_type(self):
        # Test the first value of the output to see if it's a string.
        method_call = type(cm.sort_fasta('../data/working_dir/fasta_directory/trimmed.fasta')[0])
        expected = type("Lorem ipsum dolor")
        self.assertEqual(method_call, expected)

class TestWriteToFasta(unittest.TestCase):
    
    def test_output_directory_exists(self):

        # Test if the directory is extant.
        list_data = cm.sort_fasta('../data/working_dir/fasta_directory/trimmed.fasta')
        method_call = os.path.exists(cm.write_to_fasta('../data/working_dir/fasta_directory/','trimmed', list_data))
        expected = True
        self.assertEqual(method_call, expected)


class TestFirstConsensus(unittest.TestCase):
    
    def test_overall(self):

        # Test if the final output directory matches what it should be
        wdir = "../data/working_dir/"
        rdir ="../data/reads/"
        refdir = "../data/reference_fasta/xBSn_90mer_fake_randomer.fa"
        method_call = cm.first_consensus(wdir, rdir, refdir)
        expected = str(os.path.abspath('../data/working_dir/rough_consensus_output/consensus.fasta'))
        self.assertEqual(method_call, expected)


class TestExtractNIndexes(self):

    def test_outtype(self):
        fastafile = "../data/reference_fasta/xBSn_90mer_fake_randomer.fa"
        method_call = type(cm.extract_n_indexes(fastafile))
        expected = type([(1, 2)])
        self.assertEqual(method_call, expected)


class TestRenameConsensusHeaders(self):

    def test_output_exists(self):
        method_call = cm.rename_consensus_headers('../data/working_dir/rough_consensus_output/consensus.fasta', 1, 2, '../data/working_dir/rough_consensus_output/labeled_consensus.fasta')
        expected = str(os.path.abspath('../data/working_dir/rough_consensus_output/labeled_consensus.fasta'))
        self.assertEqual(method_call, expected)

if __name__ == '__main__':
    unittest.main()
