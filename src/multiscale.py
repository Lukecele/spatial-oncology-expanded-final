# src/multiscale.py

"""Multiscale integration module.
Couples the PK/PD compartmental model with the reaction-diffusion PDE solver
to produce a Cellular Response Index (CRI).
"""

import numpy as np
from . import pk_pd
from . import reaction_diffusion as rd


def compute_cri(dose, k_elim, V, D, k_decay, times, L=1.0, nx=101):
    """Compute the Cellular Response Index (CRI) over time.

    CRI(t) integrates the spatial drug concentration profile weighted by
    a cell-kill function (Hill model) at each time point.

    Parameters
    ----------
    dose : float – administered dose
    k_elim : float – elimination rate constant
    V : float – volume of distribution
    D : float – tissue diffusion coefficient
    k_decay : float – tissue decay rate
    times : array-like – time points (hours)
    L : float – tissue domain length (mm)
    nx : int – spatial resolution

    Returns
    -------
    cri_values : np.ndarray – CRI at each time point
    """
    times = np.asarray(times, dtype=float)
    cri_values = np.zeros(len(times))

    for i, t in enumerate(times):
        # Plasma concentration at time t (boundary condition for PDE)
        c_plasma = pk_pd.simulate_one_compartment(dose, k_elim, V, [t])[0]

        # Solve spatial drug distribution at time t
        x, c_tissue = rd.reaction_diffusion_1d(
            D=D, k_decay=k_decay, c0=c_plasma, L=L, nx=nx, T=max(t, 0.01)
        )

        # Hill-model cell-kill: E(c) = c^2 / (IC50^2 + c^2)
        IC50 = 0.5  # default IC50 in µM
        effect = c_tissue**2 / (IC50**2 + c_tissue**2)

        # CRI = spatial integral of effect over tissue
        cri_values[i] = np.trapezoid(effect, x)

    return cri_values


def compute_cri_3d(dose, k_elim, V, D, k_decay, t_eval,
                    grid_size=20, L=1.0, IC50=0.5):
    """Compute CRI using 3-D reaction-diffusion."""
    c_plasma = pk_pd.simulate_one_compartment(dose, k_elim, V, [t_eval])[0]
    centre = grid_size // 2
    c3d = rd.reaction_diffusion_3d(
        D=D, k_decay=k_decay,
        source_coords=[(centre, centre, centre)],
        source_strength=c_plasma,
        grid_size=grid_size, L=L, T=max(t_eval, 0.01)
    )
    effect = c3d**2 / (IC50**2 + c3d**2)
    dx = L / (grid_size - 1)
    cri = np.sum(effect) * dx**3
    return cri


if __name__ == "__main__":
    times = np.linspace(0.5, 24, 24)
    cri = compute_cri(dose=100, k_elim=0.15, V=10,
                       D=0.01, k_decay=0.05, times=times)
    print("CRI over time:")
    for t, v in zip(times, cri):
        print(f"  t={t:5.1f}h  CRI={v:.4f}")
