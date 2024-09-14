import pytest
# from unittest import mock
# from pathlib import Path
from heatrapy.dimension_1.objects import Object as Object1D


# Mocking mats.CalMatPro class
class MockCalMatPro:
    def __init__(
        self, tadi, tadd, cpa, cp0, k0, ka, rho0, rhoa, lheat0, lheata
    ):
        pass

    def tadi(self, temp):
        return 10.0

    def tadd(self, temp):
        return 5.0

    def cpa(self, temp):
        return 1.0

    def cp0(self, temp):
        return 0.5

    def ka(self, temp):
        return 0.05

    def k0(self, temp):
        return 0.02

    def rho0(self, temp):
        return 1000

    def rhoa(self, temp):
        return 1050

    def lheat0(self):
        return [[273, 100]]

    def lheata(self):
        return [[273, 100]]


# Mocking the mats import
@pytest.fixture(autouse=True)
def mock_mats(monkeypatch):
    monkeypatch.setattr("heatrapy.mats.CalMatPro", MockCalMatPro)


# Tests
def test_object_initialization():
    """Test the Object class initialization with valid input."""
    obj = Object1D(
        amb_temperature=300,
        materials=('Cu',),
        borders=(1, 5),
        materials_order=(0,),
        dx=0.01,
        dt=0.1,
    )
    assert obj.amb_temperature == 300
    assert obj.materials_index == [None, 0, 0, 0, 0, None]
    assert obj.temperature == [
        [300, 300],
        [300, 300],
        [300, 300],
        [300, 300],
        [300, 300],
        [300, 300]
    ]


def test_object_invalid_initialization():
    """Test the Object class initialization with invalid input types."""
    with pytest.raises(ValueError):
        Object1D(amb_temperature="300")  # Invalid ambient temperature type


def test_object_activation():
    """Test the activation method."""
    obj = Object1D(
        amb_temperature=300,
        materials=('Cu',),
        borders=(1, 5),
        materials_order=(0,),
        dx=0.01,
        dt=0.1,
    )

    obj.activate(1, 3)

    assert obj.state[1] is True
    assert obj.state[2] is True
    assert obj.temperature[1][0] == 310  # tadi adds 10
    assert obj.temperature[2][0] == 310  # tadi adds 10


def test_object_deactivation():
    """Test the deactivation method."""
    obj = Object1D(
        amb_temperature=300,
        materials=('Cu',),
        borders=(1, 5),
        materials_order=(0,),
        dx=0.01,
        dt=0.1,
    )

    # First activate some points
    obj.activate(1, 3)

    # Then deactivate them
    obj.deactivate(1, 3)

    assert obj.state[1] is False
    assert obj.state[2] is False
    assert obj.temperature[1][0] == 305  # tadd subtracts 5
    assert obj.temperature[2][0] == 305  # tadd subtracts 5


def test_object_boundaries():
    """Test boundary conditions."""
    obj = Object1D(
        amb_temperature=300,
        materials=('Cu',),
        borders=(1, 5),
        materials_order=(0,),
        boundaries=(1, 0)  # Testing boundaries other than default (0, 0)
    )

    assert obj.boundaries == (1, 0)
