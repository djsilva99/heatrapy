import pytest
from heatrapy.dimension_1.objects.single import SingleObject as SingleObject1D


@pytest.fixture
def example_single_object_implicit():
    return SingleObject1D(300, boundaries=(0, 200))


@pytest.fixture
def example_single_object_explicit():
    return SingleObject1D(300, boundaries=(0, 200), draw=[])


@pytest.fixture
def example_single_object_with_materials():
    return SingleObject1D(
        300, boundaries=(0, 200), materials=('Cu', 'Gd_mag'),
        borders=(1, 11, 22), materials_order=(0, 1), draw=[]
    )


def test_implicit_general_solver(example_single_object_implicit):
    """Test singleObject with the implicit_general solver."""
    # given
    solution = 245

    # when
    example_single_object_implicit.activate(2, 3)
    example_single_object_implicit.deactivate(2, 3)
    example_single_object_implicit.compute(
        30,
        5,
        solver='implicit_general'
    )

    # then
    assert int(
        example_single_object_implicit.object.temperature[5][0]
    ) == solution


def test_explicit_general_solver(example_single_object_explicit):
    """Test singleObject with the explicit_general solver."""
    # given
    solution = 245

    # when
    example_single_object_explicit.compute(
        30,
        5,
        solver='explicit_general'
    )
    example_single_object_explicit.show_figure("temperature")

    # then
    assert int(
        example_single_object_explicit.object.temperature[5][0]
    ) == solution


def test_explicit_k_solver(example_single_object_with_materials):
    """Test singleObject with the explicit_k(k) solver."""
    # given
    solution = 270

    # when
    example_single_object_with_materials.compute(
        2000,
        5,
        solver='explicit_k(x)'
    )

    # then
    assert int(
        example_single_object_with_materials.object.temperature[5][0]
    ) == solution


def test_implicit_k_solver(example_single_object_with_materials):
    """Test singleObject with the implicit_k(k) solver."""
    # given
    solution = 270

    # when
    example_single_object_with_materials.compute(
        2000,
        5,
        solver='implicit_k(x)'
    )
    example_single_object_with_materials.change_boundaries((0, 0))

    # then
    assert int(
        example_single_object_with_materials.object.temperature[5][0]
    ) == solution
