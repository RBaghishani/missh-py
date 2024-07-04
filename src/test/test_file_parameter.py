import unittest
from file_parameter import FileParameter

class TestFileParameter(unittest.TestCase):
    def setUp(self):
        self.param = FileParameter()

    def test_init_single_end(self):
        result = self.param.init("file1.fa")
        self.assertTrue(result)

    def test_add_spaced_qmer(self):
        self.param.add_spaced_qmer("1**1", "1**1")
        self.assertEqual(len(self.param.spaced_qmers), 1)
        self.assertEqual(self.param.spaced_qmers[0], ("1**1", "1**1"))

if __name__ == '__main__':
    unittest.main()
