import pytest
from heatrapy.dimension_1.objects import SystemObjects as SystemObjects1D


@pytest.fixture
def example_system_object_1d():
    example = SystemObjects1D()
    example.objects[0].temperature = [[200., 200.] for i in range(12)]
    example.contact_add(((0, 3), (1, 3), 1e6))
    return example


def test_implicit_general_system_1d(example_system_object_1d):
    solution = 246
    example_system_object_1d.compute(200, 5, solver='implicit_general')
    assert int(
        example_system_object_1d.objects[0].temperature[5][0]
    ) == solution


def test_explicit_general_system_1d(example_system_object_1d):
    solution = 246
    example_system_object_1d.compute(200, 5, solver='explicit_general')
    assert int(
        example_system_object_1d.objects[0].temperature[5][0]
    ) == solution


def test_implicit_k_system_1d(example_system_object_1d):
    solution = 246
    example_system_object_1d.compute(200, 5, solver='implicit_k(x)')
    assert int(
        example_system_object_1d.objects[0].temperature[5][0]
    ) == solution


def test_explicit_k_system_1d(example_system_object_1d):
    solution = 246
    example_system_object_1d.compute(200, 5, solver='explicit_k(x)')
    assert int(
        example_system_object_1d.objects[0].temperature[5][0]
    ) == solution
