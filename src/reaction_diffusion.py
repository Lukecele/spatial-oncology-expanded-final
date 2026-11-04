# src/reaction_diffusion.py

"""Reaction-Diffusion solver for drug transport in tumour tissue.
Uses a simple finite-difference scheme (no FEniCS dependency for portability).
"""

import numpy as np

def reaction_diffusion_1d(D, k_decay, c0, L=1.0, nx=101, T=1.0, nt=1000):
    """Solve 1-D reaction-diffusion: dc/dt = D * d²c/dx² - k_decay * c
    Parameters
    ----------
    D : float – diffusion coefficient (mm²/h)
    k_decay : float – first-order decay rate (1/h)
    c0 : float – boundary concentration at x=0
    L : float – domain length (mm)
    nx : int – number of spatial grid points
    T : float – total simulation time (h)
    nt : int – number of time steps
    Returns
    -------
    x : np.ndarray – spatial grid
    c : np.ndarray – concentration profile at time T
    """
    dx = L / (nx - 1)
    dt = T / nt
    x = np.linspace(0, L, nx)
    c = np.zeros(nx)
    c[0] = c0

    for _ in range(nt):
        c_new = c.copy()
        for i in range(1, nx - 1):
            c_new[i] = c[i] + dt * (
                D * (c[i+1] - 2*c[i] + c[i-1]) / dx**2 - k_decay * c[i]
            )
        c_new[0] = c0      # Dirichlet BC at x=0
        c_new[-1] = 0.0    # Dirichlet BC at x=L
        c = c_new
    return x, c


def reaction_diffusion_3d(D, k_decay, source_coords, source_strength,
                          grid_size=20, L=1.0, T=0.5, nt=500):
    """Solve 3-D reaction-diffusion on a uniform grid with point sources.
    Parameters
    ----------
    D : float – diffusion coefficient
    k_decay : float – decay rate
    source_coords : list of (i,j,k) grid indices for drug sources
    source_strength : float – constant source term
    grid_size : int – number of points per axis
    L : float – domain side length
    T : float – simulation time
    nt : int – time steps
    Returns
    -------
    c : np.ndarray of shape (grid_size, grid_size, grid_size)
    """
    n = grid_size
    dx = L / (n - 1)
    dt = T / nt
    c = np.zeros((n, n, n))

    for _ in range(nt):
        laplacian = np.zeros_like(c)
        laplacian[1:-1,1:-1,1:-1] = (
            (c[2:,1:-1,1:-1] + c[:-2,1:-1,1:-1]
             + c[1:-1,2:,1:-1] + c[1:-1,:-2,1:-1]
             + c[1:-1,1:-1,2:] + c[1:-1,1:-1,:-2]
             - 6*c[1:-1,1:-1,1:-1]) / dx**2
        )
        c += dt * (D * laplacian - k_decay * c)
        for (si, sj, sk) in source_coords:
            c[si, sj, sk] += dt * source_strength
    return c


if __name__ == "__main__":
    x, c = reaction_diffusion_1d(D=0.01, k_decay=0.05, c0=1.0)
    print(f"1-D solution: max={c.max():.4f}, min={c.min():.4f}")
    c3 = reaction_diffusion_3d(D=0.005, k_decay=0.02,
                                source_coords=[(10,10,10)],
                                source_strength=1.0)
    print(f"3-D solution: max={c3.max():.4f}")
