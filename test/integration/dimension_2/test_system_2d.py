import pytest
from heatrapy.dimension_2.objects import SystemObjects as SystemObjects2D


@pytest.fixture
def example_system_object_2d():
    example = SystemObjects2D(boundaries=((0, 0, 0, 0), (250, 0, 0, 0)))
    example.contact_add(((0, (4, 4)), (1, (7, 7)), 100000))
    return example


def test_explicit_general_system_2d(example_system_object_2d):
    solution = 275
    example_system_object_2d.compute(10, 5, solver='explicit_general')
    assert int(
        example_system_object_2d.objects[1].temperature[3][3][0]
    ) == solution


def test_explicit_k_system_2d(example_system_object_2d):
    solution = 275
    example_system_object_2d.compute(10, 5, solver='explicit_k(x)')
    assert int(
        example_system_object_2d.objects[1].temperature[3][3][0]
    ) == solution
