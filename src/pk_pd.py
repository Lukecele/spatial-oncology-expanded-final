# src/pk_pd.py

"""PK/PD modelling utilities.
Implements basic compartmental model fitting using Stan via CmdStanPy.
"""

import numpy as np
import pandas as pd


def simulate_one_compartment(dose, k_elim, V, times):
    """Simple one-compartment IV bolus model.
    C(t) = (dose/V) * exp(-k_elim * t)
    """
    dose = float(dose)
    k_elim = float(k_elim)
    V = float(V)
    times = np.asarray(times, dtype=float)
    conc = (dose / V) * np.exp(-k_elim * times)
    return conc


def fit_one_compartment(data_df, dose, V_initial=1.0):
    """Fit a one-compartment model to concentration-time data.
    Parameters
    ----------
    data_df: pandas.DataFrame with columns ['time','conc']
    dose: administered dose (same units as concentration*V)
    V_initial: initial guess for volume of distribution
    Returns
    -------
    fit_result: CmdStanModel fit object
    """
    try:
        from cmdstanpy import CmdStanModel
    except ImportError:
        raise ImportError(
            "cmdstanpy is required for Bayesian fitting. "
            "Install it with: pip install cmdstanpy"
        )

    stan_code = """
    data {
        int<lower=1> N;               // number of observations
        vector[N] time;                // sampling times
        vector[N] conc;                // observed concentrations
        real<lower=0> dose;            // dose administered
        real<lower=0> V_guess;         // initial V guess
    }
    parameters {
        real<lower=0> k_elim;         // elimination rate constant
        real<lower=0> V;              // volume of distribution
        real<lower=0> sigma;          // observation noise
    }
    model {
        // priors (weakly informative)
        k_elim ~ normal(0.1, 0.5);
        V ~ normal(V_guess, 5.0);
        sigma ~ normal(0, 1);
        // likelihood
        conc ~ normal((dose / V) * exp(-k_elim * time), sigma);
    }
    """
    import tempfile, os
    with tempfile.TemporaryDirectory() as tmpdir:
        stan_path = os.path.join(tmpdir, "pk_pd_model.stan")
        with open(stan_path, "w") as f:
            f.write(stan_code)
        model = CmdStanModel(stan_file=stan_path)
        data = {
            "N": len(data_df),
            "time": data_df["time"].values,
            "conc": data_df["conc"].values,
            "dose": dose,
            "V_guess": V_initial,
        }
        fit = model.sample(data=data, chains=4, parallel_chains=4,
                           iter_sampling=1000, iter_warmup=500)
    return fit


if __name__ == "__main__":
    times = np.linspace(0, 24, 13)
    conc = simulate_one_compartment(dose=100, k_elim=0.15, V=10, times=times)
    df = pd.DataFrame({"time": times, "conc": conc})
    fit = fit_one_compartment(df, dose=100, V_initial=10)
    print(fit.summary())
