"""Unit tests for 1D solver type checking.

Verifies that solvers raise TypeError for invalid inputs.
"""

import pytest

from heatrapy.dimension_1.solvers._latent_heat import (
    apply_latent_heat,
)
from heatrapy.dimension_1.solvers.explicit_general import (
    explicit_general,
)
from heatrapy.dimension_1.solvers.explicit_k import explicit_k
from heatrapy.dimension_1.solvers.implicit_general import (
    implicit_general,
)
from heatrapy.dimension_1.solvers.implicit_k import implicit_k


class TestApplyLatentHeatTypeErrors:
    """TypeError checks for apply_latent_heat."""

    def test_nx_not_list_raises(self):
        """nx must be a list, not a tuple or ndarray."""

        class FakeObj:
            num_points = 3

        with pytest.raises(TypeError, match="nx must be a list"):
            apply_latent_heat((1.0, 2.0, 3.0), FakeObj())

    def test_obj_not_object_raises(self):
        """obj must have num_points attribute."""
        with pytest.raises(
            TypeError, match="obj must be a thermal Object"
        ):
            apply_latent_heat([1.0, 2.0], "not_an_object")


class TestSolverTypeErrors:
    """TypeError checks for all four solver functions."""

    @pytest.mark.parametrize("solver", [
        explicit_general,
        explicit_k,
        implicit_general,
        implicit_k,
    ])
    def test_string_raises(self, solver):
        with pytest.raises(
            TypeError, match="obj must be a thermal Object"
        ):
            solver("not_an_object")

    @pytest.mark.parametrize("solver", [
        explicit_general,
        explicit_k,
        implicit_general,
        implicit_k,
    ])
    def test_none_raises(self, solver):
        with pytest.raises(
            TypeError, match="obj must be a thermal Object"
        ):
            solver(None)

    @pytest.mark.parametrize("solver", [
        explicit_general,
        explicit_k,
        implicit_general,
        implicit_k,
    ])
    def test_int_raises(self, solver):
        with pytest.raises(
            TypeError, match="obj must be a thermal Object"
        ):
            solver(42)
