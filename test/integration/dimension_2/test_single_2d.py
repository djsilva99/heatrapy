import pytest
from heatrapy.dimension_2.objects import SingleObject as SingleObject2D


@pytest.fixture
def example_single_object_2d():
    return SingleObject2D(
        300,
        boundaries=(200, 0, 0, 0),
        draw=[
            "materials", "temperature", "state", "Q", "Q0"
        ]
    )


@pytest.fixture
def example_single_object_2d_activated():
    return SingleObject2D(
        300,
        boundaries=(200, 0, 0, 0),
        draw=[
            "materials", "temperature", "state", "Q", "Q0"
        ],
        initial_state=True
    )


def test_explicit_general_single_2d(example_single_object_2d):
    solution = 245
    example_single_object_2d.compute(30, 5, solver='explicit_general')
    example_single_object_2d.show_figure("temperature")
    example_single_object_2d.show_figure("materials")
    example_single_object_2d.show_figure("state")
    example_single_object_2d.show_figure("Q")
    example_single_object_2d.show_figure("Q0")
    assert int(
        example_single_object_2d.object.temperature[5][5][0]
    ) == solution


def test_explicit_k_single_2d(example_single_object_2d_activated):
    solution = 245
    example_single_object_2d_activated.activate((4, 4), 2, shape="circle")
    example_single_object_2d_activated.deactivate((4, 4), 2, shape="circle")
    example_single_object_2d_activated.activate((2, 2), (3, 3))
    example_single_object_2d_activated.deactivate((2, 2), (3, 3))
    example_single_object_2d_activated.compute(30, 5, solver='explicit_k(x)')

    assert int(
        example_single_object_2d_activated.object.temperature[5][5][0]
    ) == solution


def test_explicit_k_single_2d_with_activation(example_single_object_2d):
    solution = 245
    example_single_object_2d.activate((4, 4), 2, shape="circle")
    example_single_object_2d.deactivate((4, 4), 2, shape="circle")
    example_single_object_2d.activate((2, 2), (3, 3))
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
    example_single_object_2d.change_material("circle", "Cu", (4, 4), 2)
    example_single_object_2d.compute(30, 5, solver='explicit_k(x)')
    assert int(
        example_single_object_2d.object.temperature[5][5][0]
    ) == solution


def test_explicit_k_single_2d_with_power_change(example_single_object_2d):
    solution = 245
    example_single_object_2d.change_power(
        "square", "Q", (1, 1), (2, 2), 0
    )
    example_single_object_2d.change_power(
        "circle", "Q", (4, 4), 2, 0
    )
    example_single_object_2d.change_power(
        "square", "Q0", (1, 1), (2, 2), 0
    )
    example_single_object_2d.change_power(
        "circle", "Q0", (4, 4), 2, 0
    )
    example_single_object_2d.compute(30, 5, solver='explicit_k(x)')
    assert int(
        example_single_object_2d.object.temperature[5][5][0]
    ) == solution
