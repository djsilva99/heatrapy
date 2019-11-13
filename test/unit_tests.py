"""unit_tests.

Contains the unit tests for the heatrapy modules
"""

import unittest
from .. import heatrapy as ht


class SingleObjects(unittest.TestCase):
    """Test all single_objects components."""

    def test_implicit_k(self):
        """Test single_object with the implicit_k(k) solver."""
        solution = 278
        example = ht.single_object(300, boundaries=(0, 200))
        example.compute(10, 5)
        self.assertEqual(int(example.temperature[5][0]), solution)


if __name__ == '__main__':
    unittest.main()
