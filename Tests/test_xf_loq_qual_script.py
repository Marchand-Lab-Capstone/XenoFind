import sys
sys.path.append('../lib/')
import xf_low_qual_script as xflq
import unittest

class TestGetDirectory(unittest.TestCase):
    def test1(self):
        self.assertEqual(xflq.get_directory(#params),#response)
            
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
