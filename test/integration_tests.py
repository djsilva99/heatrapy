"""integration_tests.

Contains the integration tests for the heatrapy modules
"""

import time
import unittest
from .. import heatrapy as htp


class SolidActiveRegenerator(unittest.TestCase):
    """Test solid_active_regenerator."""

    def test_main(self):
        """Test solid_active_regenerator with the implicit_k(x) solver."""
        solution = 751
        file_name = str(time.time())[:10]
        example = htp.solid_active_regenerator_1D(file_name, max_cycle_number=2)
        self.assertEqual(int(example[0]), solution)


class FluidActiveRegenerator(unittest.TestCase):
    """Test fluid_active_regenerator."""

    def test_main(self):
        """Test fluid_active_regenerator with the implicit_k(x) solver."""
        solution = 88455
        file_name = str(time.time())[:10]
        example = htp.fluid_active_regenerator_1D(
            file_name, max_cycle_number=2, type_study='fixed_temperature_span',
            h_mcm_fluid=1.8e6, h_leftreservoir_fluid=1.8e6,
            h_rightreservoir_fluid=1.8e6
        )
        self.assertEqual(int(example[1]), solution)


if __name__ == '__main__':
    unittest.main()
