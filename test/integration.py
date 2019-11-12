"""integration_tests.

Contains the integration tests for the heatrapy modules
"""

import unittest
import sys
sys.path.append('../')
import heatrapy as ht


class SolidActiveRegenerator(unittest.TestCase):
    """Test solid_active_regenerator."""

    def main(self):
        """Test solid_active_regenerator with the implicit_k(k) solver."""
        solution = 278
        example = ht.single_object(boundaries=[0, 200])
        self.assertEqual(int(round(example.compute(10)[5])), solution)


if __name__ == '__main__':
    unittest.main()
