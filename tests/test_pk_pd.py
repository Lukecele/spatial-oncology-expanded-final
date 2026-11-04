"""Unit tests for src/pk_pd.py"""

import numpy as np
import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.pk_pd import simulate_one_compartment


class TestSimulateOneCompartment:
    def test_initial_concentration(self):
        """At t=0, C(0) = dose / V."""
        c = simulate_one_compartment(dose=100, k_elim=0.15, V=10, times=[0])
        assert np.isclose(c[0], 10.0), f"Expected 10.0, got {c[0]}"

    def test_decay(self):
        """Concentration should decrease over time."""
        times = [0, 1, 2, 5, 10]
        c = simulate_one_compartment(dose=100, k_elim=0.15, V=10, times=times)
        for i in range(len(c) - 1):
            assert c[i] > c[i + 1], "Concentration should decrease monotonically"

    def test_half_life(self):
        """Verify half-life: t_half = ln(2) / k_elim."""
        k_elim = 0.15
        t_half = np.log(2) / k_elim
        c = simulate_one_compartment(dose=100, k_elim=k_elim, V=10, times=[0, t_half])
        assert np.isclose(c[1], c[0] / 2, rtol=1e-6), \
            f"At t_half={t_half:.2f}, C should be half of C(0)"

    def test_zero_dose(self):
        """Zero dose gives zero concentration."""
        c = simulate_one_compartment(dose=0, k_elim=0.15, V=10, times=[0, 1, 5])
        assert all(np.isclose(c, 0)), "Zero dose should give zero concentration"

    def test_output_length(self):
        """Output length matches input times."""
        times = np.linspace(0, 24, 50)
        c = simulate_one_compartment(dose=100, k_elim=0.15, V=10, times=times)
        assert len(c) == len(times)

    def test_non_negative(self):
        """Concentrations should always be non-negative."""
        times = np.linspace(0, 100, 200)
        c = simulate_one_compartment(dose=100, k_elim=0.5, V=5, times=times)
        assert np.all(c >= 0), "Concentrations must be non-negative"
