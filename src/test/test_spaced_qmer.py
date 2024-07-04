import unittest
from spaced.spaced_qmer import SpacedQmer

class TestSpacedQmer(unittest.TestCase):
    def setUp(self):
        self.qmer = SpacedQmer()

    def test_initialize(self):
        self.qmer.initialize("1**1")
        self.assertEqual(self.qmer.get_pattern(), "1**1")

if __name__ == '__main__':
    unittest.main()
