import unittest
from parameter.parameter import Parameter

class TestParameter(unittest.TestCase):
    def setUp(self):
        self.param = Parameter()

    def test_init(self):
        result = self.param.init("param_file.txt")
        self.assertTrue(result)

    def test_add_spaced_qmer(self):
        self.param.add_spaced_qmer("1**1", "1**1")
        self.assertEqual(len(self.param.spaced_qmers), 1)
        self.assertEqual(self.param.spaced_qmers[0], ("1**1", "1**1"))

if __name__ == '__main__':
    unittest.main()
