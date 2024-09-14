import pytest
from heatrapy.dimension_2.objects import SingleObject as SingleObject2D


@pytest.fixture
def example_single_object_2d():
    return SingleObject2D(300, boundaries=(200, 0, 0, 0), draw=[])


def test_explicit_general_single_2d(example_single_object_2d):
    solution = 245
    example_single_object_2d.compute(30, 5, solver='explicit_general')
    assert int(
        example_single_object_2d.object.temperature[5][5][0]
    ) == solution


def test_explicit_k_single_2d(example_single_object_2d):
    solution = 245
    example_single_object_2d.compute(30, 5, solver='explicit_k(x)')
    assert int(
        example_single_object_2d.object.temperature[5][5][0]
    ) == solution


def test_explicit_k_single_2d_with_activation(example_single_object_2d):
    solution = 245
    example_single_object_2d.activate((2, 2), (3, 3))
    example_single_object_2d.deactivate((2, 2), (3, 3))
    example_single_object_2d.compute(30, 5, solver='explicit_k(x)')
    assert int(
        example_single_object_2d.object.temperature[5][5][0]
    ) == solution


def test_explicit_k_single_2d_with_bondary_change(example_single_object_2d):
    solution = 245
    example_single_object_2d.change_boundaries((200, 0, 0, 0))
    example_single_object_2d.compute(30, 5, solver='explicit_k(x)')
    assert int(
        example_single_object_2d.object.temperature[5][5][0]
    ) == solution


def test_explicit_k_single_2d_with_material_change(example_single_object_2d):
    solution = 245
    example_single_object_2d.change_material("square", "Cu", (1, 1), (2, 2))
    example_single_object_2d.compute(30, 5, solver='explicit_k(x)')
    assert int(
        example_single_object_2d.object.temperature[5][5][0]
    ) == solution


def test_explicit_k_single_2d_with_power_change(example_single_object_2d):
    solution = 245
    example_single_object_2d.change_power(
        "square", "Q0", (1, 1), (2, 2), 0
    )
    example_single_object_2d.compute(30, 5, solver='explicit_k(x)')
    assert int(
        example_single_object_2d.object.temperature[5][5][0]
    ) == solution
