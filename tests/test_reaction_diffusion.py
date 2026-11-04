"""Unit tests for src/reaction_diffusion.py"""

import numpy as np
import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.reaction_diffusion import reaction_diffusion_1d, reaction_diffusion_3d


class TestReactionDiffusion1D:
    def test_boundary_conditions(self):
        """c(0) = c0, c(L) = 0."""
        x, c = reaction_diffusion_1d(D=0.01, k_decay=0.05, c0=1.0)
        assert np.isclose(c[0], 1.0), f"Left BC violated: {c[0]}"
        assert np.isclose(c[-1], 0.0), f"Right BC violated: {c[-1]}"

    def test_monotonic_decay(self):
        """With decay and no internal source, concentration should decrease from left."""
        x, c = reaction_diffusion_1d(D=0.01, k_decay=0.05, c0=1.0, T=2.0, nt=5000)
        # Allow small numerical noise
        diffs = np.diff(c)
        assert np.all(diffs <= 1e-10), "Concentration should decrease monotonically from x=0"

    def test_higher_diffusion_spreads_more(self):
        """Higher D should yield higher concentration in the middle."""
        _, c_low = reaction_diffusion_1d(D=0.001, k_decay=0.01, c0=1.0, T=2.0, nt=5000)
        _, c_high = reaction_diffusion_1d(D=0.1, k_decay=0.01, c0=1.0, T=2.0, nt=5000)
        mid = len(c_low) // 2
        assert c_high[mid] > c_low[mid], "Higher diffusion should spread drug further"

    def test_no_decay_preserves_more(self):
        """k_decay=0 should preserve more drug than k_decay>0."""
        _, c_nodecay = reaction_diffusion_1d(D=0.01, k_decay=0.0, c0=1.0, T=1.0, nt=2000)
        _, c_decay = reaction_diffusion_1d(D=0.01, k_decay=0.1, c0=1.0, T=1.0, nt=2000)
        assert np.sum(c_nodecay) > np.sum(c_decay)

    def test_output_shape(self):
        """Output arrays should have correct shape."""
        nx = 51
        x, c = reaction_diffusion_1d(D=0.01, k_decay=0.05, c0=1.0, nx=nx)
        assert len(x) == nx
        assert len(c) == nx

    def test_non_negative(self):
        """Concentrations must be non-negative."""
        x, c = reaction_diffusion_1d(D=0.01, k_decay=0.05, c0=1.0, T=5.0, nt=10000)
        assert np.all(c >= -1e-12), "Concentrations must be non-negative"


class TestReactionDiffusion3D:
    def test_output_shape(self):
        """3-D output should have shape (n, n, n)."""
        n = 10
        c = reaction_diffusion_3d(D=0.005, k_decay=0.02,
                                   source_coords=[(5,5,5)],
                                   source_strength=1.0, grid_size=n)
        assert c.shape == (n, n, n)

    def test_source_has_highest_concentration(self):
        """Source location should have the highest concentration."""
        n = 10
        c = reaction_diffusion_3d(D=0.005, k_decay=0.02,
                                   source_coords=[(5,5,5)],
                                   source_strength=1.0, grid_size=n)
        assert c[5,5,5] == c.max()

    def test_non_negative(self):
        """All concentrations non-negative."""
        c = reaction_diffusion_3d(D=0.005, k_decay=0.02,
                                   source_coords=[(5,5,5)],
                                   source_strength=1.0, grid_size=10)
        assert np.all(c >= -1e-12)
