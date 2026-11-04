"""Unit tests for src/multiscale.py"""

import numpy as np
import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.multiscale import compute_cri, compute_cri_3d


class TestComputeCRI:
    def test_output_length(self):
        """CRI output should match number of time points."""
        times = np.linspace(0.5, 12, 12)
        cri = compute_cri(dose=100, k_elim=0.15, V=10, D=0.01, k_decay=0.05, times=times)
        assert len(cri) == len(times)

    def test_non_negative(self):
        """CRI values must be non-negative."""
        times = np.linspace(0.5, 24, 20)
        cri = compute_cri(dose=100, k_elim=0.15, V=10, D=0.01, k_decay=0.05, times=times)
        assert np.all(cri >= 0)

    def test_higher_dose_higher_cri(self):
        """Higher dose should produce higher CRI."""
        times = [1.0, 2.0, 4.0]
        cri_low = compute_cri(dose=10, k_elim=0.15, V=10, D=0.01, k_decay=0.05, times=times)
        cri_high = compute_cri(dose=1000, k_elim=0.15, V=10, D=0.01, k_decay=0.05, times=times)
        assert np.all(cri_high >= cri_low), "Higher dose should yield higher or equal CRI"

    def test_cri_decreases_over_time(self):
        """With elimination, CRI should generally decrease at later times."""
        times = np.linspace(1, 48, 20)
        cri = compute_cri(dose=100, k_elim=0.3, V=10, D=0.01, k_decay=0.05, times=times)
        # Last value should be less than first
        assert cri[-1] < cri[0], "CRI should decrease as drug is eliminated"


class TestComputeCRI3D:
    def test_non_negative(self):
        """3-D CRI must be non-negative."""
        cri = compute_cri_3d(dose=100, k_elim=0.15, V=10, D=0.01, k_decay=0.05, t_eval=2.0)
        assert cri >= 0

    def test_scalar_output(self):
        """3-D CRI returns a scalar."""
        cri = compute_cri_3d(dose=100, k_elim=0.15, V=10, D=0.01, k_decay=0.05, t_eval=1.0)
        assert np.isscalar(cri)
