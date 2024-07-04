import unittest
from input.input import Input

class TestInput(unittest.TestCase):
    def setUp(self):
        self.input = Input()

    def test_init(self):
        result = self.input.init("file1.fa", "file2.fa")
        self.assertTrue(result)

if __name__ == '__main__':
    unittest.main()
