"""unit_tests.

Contains the unit tests for the heatrapy modules
"""

import unittest
from .. import heatrapy as ht


class SingleObjects(unittest.TestCase):
    """Test all single_objects components."""

    def test_implicit_general(self):
        """Test singleObject with the implicit_general solver."""
        solution = 245
        example = ht.SingleObject(300, boundaries=(0, 200))
        example.compute(30, 5, solver='implicit_general')
        self.assertEqual(int(example.temperature[5][0]), solution)

    def test_explicit_general(self):
        """Test singleObject with the explicit_general solver."""
        solution = 245
        example = ht.SingleObject(300, boundaries=(0, 200))
        example.compute(30, 5, solver='explicit_general')
        self.assertEqual(int(example.temperature[5][0]), solution)

    def test_explicit_k(self):
        """Test singleObject with the explicit_k(k) solver."""
        solution = 270
        example = ht.SingleObject(
            300, boundaries=(0, 200), materials=('Cu', 'Gd'),
            borders=(1, 11, 22), materials_order=(0, 1)
        )
        example.compute(2000, 5, solver='explicit_k(x)')
        self.assertEqual(int(example.temperature[5][0]), solution)

    def test_implicit_k(self):
        """Test singleObject with the implicit_k(k) solver."""
        solution = 270
        example = ht.SingleObject(
            300, boundaries=(0, 200), materials=('Cu', 'Gd'),
            borders=(1, 11, 22), materials_order=(0, 1)
        )
        example.compute(2000, 5, solver='implicit_k(x)')
        self.assertEqual(int(example.temperature[5][0]), solution)


class SystemObjects(unittest.TestCase):
    """Test all system_objects components."""

    def test_implicit_general(self):
        """Test systemObjects with the implicit_general solver."""
        solution = 246
        example = ht.SystemObjects()
        example.objects[0].temperature = [[200., 200.] for i in range(12)]
        example.contact_add(((0, 3), (1, 3), 1e6))
        example.compute(200, 5, solver='implicit_general')
        self.assertEqual(int(example.objects[0].temperature[5][0]), solution)

    def test_explicit_general(self):
        """Test systemObjects with the explicit_general solver."""
        solution = 246
        example = ht.SystemObjects()
        example.objects[0].temperature = [[200., 200.] for i in range(12)]
        example.contact_add(((0, 3), (1, 3), 1e6))
        example.compute(200, 5, solver='explicit_general')
        self.assertEqual(int(example.objects[0].temperature[5][0]), solution)

    def test_implicit_k(self):
        """Test systemObjects with the implicit_k(x) solver."""
        solution = 246
        example = ht.SystemObjects()
        example.objects[0].temperature = [[200., 200.] for i in range(12)]
        example.contact_add(((0, 3), (1, 3), 1e6))
        example.compute(200, 5, solver='implicit_k(x)')
        self.assertEqual(int(example.objects[0].temperature[5][0]), solution)

    def test_explicit_k(self):
        """Test systemObjects with the explicit_k(x) solver."""
        solution = 246
        example = ht.SystemObjects()
        example.objects[0].temperature = [[200., 200.] for i in range(12)]
        example.contact_add(((0, 3), (1, 3), 1e6))
        example.compute(200, 5, solver='explicit_k(x)')
        self.assertEqual(int(example.objects[0].temperature[5][0]), solution)


if __name__ == '__main__':
    unittest.main()
