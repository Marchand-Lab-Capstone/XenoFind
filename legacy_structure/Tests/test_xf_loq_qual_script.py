import sys
sys.path.append('../lib/')
import xf_low_qual_script as xlfq
import unittest

# These are valid directories
working_dir = 'working_directory/'
raw_dir = '../Data/sample_fast5/subset_fast5' # currently contains both fast5 and pod5 conversions
ref_dir = '../Data/sample_fast5/fasta_metadata/fastafile.fa'
basecall_dir = '../Data/sample_fast5/merged_folder/calls.bam'


class TestGetDirectory(unittest.TestCase):
    def test_working_exists(self):
        # pass a valid directory, type working
        method_call = xlfq.get_directory(working_dir, 'working')
        expected_output = working_dir
        self.assertEqual(method_call, expected_output)

    def test_working_none(self):
        xlfq.input = lambda: 'working_dir'
        method_call = xlfq.get_directory(None, 'working')
        expected_output = working_dir
        # pass none, type working
    
    def test_raw_exists(self):
        # pass valid directory, type raw
    
    def test_raw_none(self):
        # pass none, type raw
    
    def test_ref_exists(self):
        # pass valid directory, type ref
    
    def test_ref_none(self):
        # pass none, type ref
    
    def test_base_exists(self):
        # pass valid directory, type base
    
    def test_base_none(self):
        # pass none, type ref
            
            
class TestGenerateWorkingFolders(unittest.TestCase):
    def test1(self):
         self.assertEqual(xflq.generate_working_folders(#params),#response)
         
class TestValidateDirType(unittest.TestCase):
    def test1(self):
         self.assertEqual(xflq.validate_dirype(#params),#response)
         
class TestExtractReadInfo(unittest.TestCase):
    def test1(self):
          self.assertEqual(xflq.extract_read_info(#params),#response)
    
class TestRunSript(unittest.TestCase):
    def test1(self):
         self.assertEqual(xflq.run_script(#params),#response)
