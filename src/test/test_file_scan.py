import unittest
from input.file_scan import FileScan

class TestFileScan(unittest.TestCase):
    def setUp(self):
        self.file_scan = FileScan()

    def test_scan(self):
        lines = self.file_scan.scan("test_file.txt")
        self.assertIsInstance(lines, list)

if __name__ == '__main__':
    unittest.main()
