# Import libraries
import pytest
# from pathlib import Path

# Path to the Object class (modify if your file is located elsewhere)
from heatrapy.dimension_2.objects import Object as Object2D


class TestObject:

    def test_init_valid_inputs(self):
        # Arrange
        amb_temperature = 273.15  # Room temperature in Kelvin
        material = 'Cu'
        dx = 0.01
        dy = 0.01
        dt = 0.1
        size = (10, 10)
        file_name = None
        boundaries = (0, 0, 0, 0)
        Q = []
        Q0 = []
        initial_state = False
        materials_path = False

        # Act
        obj = Object2D(amb_temperature, material, dx, dy, dt, size, file_name,
                       boundaries, Q, Q0, initial_state, materials_path)

        # Assert
        assert obj.amb_temperature == amb_temperature
        assert obj.size == size
        assert obj.dx == dx
        assert obj.dy == dy
        assert obj.dt == dt
        assert len(obj.temperature) == size[0]
        assert len(obj.temperature[0]) == size[1]

    def test_init_invalid_amb_temperature(self):
        # Arrange
        amb_temperature = "invalid"
        material = 'Cu'
        dx = 0.01
        dy = 0.01
        dt = 0.1
        size = (10, 10)
        file_name = None
        boundaries = (0, 0, 0, 0)
        Q = []
        Q0 = []
        initial_state = False
        materials_path = False

        # Act with expected error
        with pytest.raises(ValueError):
            Object2D(amb_temperature, material, dx, dy, dt, size, file_name,
                     boundaries, Q, Q0, initial_state, materials_path)

    def test_init_invalid_dx(self):
        # Arrange
        amb_temperature = 273.15
        material = 'Cu'
        dx = "invalid"
        dy = 0.01
        dt = 0.1
        size = (10, 10)
        file_name = None
        boundaries = (0, 0, 0, 0)
        Q = []
        Q0 = []
        initial_state = False
        materials_path = False

        # Act with expected error
        with pytest.raises(ValueError):
            Object2D(amb_temperature, material, dx, dy, dt, size, file_name,
                     boundaries, Q, Q0, initial_state, materials_path)

    def test_init_invalid_dt(self):
        # Arrange
        amb_temperature = 273.15
        material = 'Cu'
        dx = 0.01
        dy = 0.01
        dt = "invalid"
        size = (10, 10)
        file_name = None
        boundaries = (0, 0, 0, 0)
        Q = []
        Q0 = []
        initial_state = False
        materials_path = False

        # Act with expected error
        with pytest.raises(ValueError):
            Object2D(amb_temperature, material, dx, dy, dt, size, file_name,
                     boundaries, Q, Q0, initial_state, materials_path)
