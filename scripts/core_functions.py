# core_functions.py

import numpy as np
import pandas as pd
import warnings
from scipy.interpolate import interp1d
from scipy.integrate import quad

def load_and_prepare_potentials(ap_path, dp_path):
    """
    Loads atomic and dipole potential data and calculates the total phase shift.
    
    Returns:
        tuple: A tuple containing z-values (pd.Series) and total phase shift (pd.Series).
    """
    AP = pd.read_csv(ap_path, header=None, sep=r'\s+')
    DP = pd.read_csv(dp_path, header=None, sep=r'\s+')

    element_potential = {"H": 25, "C": 130, "N": 108, "O": 97, "P": 267}

    def element_value(symbol):
        return element_potential.get(symbol, None)

    zVals = AP.iloc[:, 0][1:].astype(float)
    first_row = AP.iloc[0]
    first_elements = [str(item)[0] for item in first_row]
    atomic_vals = [element_value(element) for element in first_elements[1:]]

    AP2 = AP.drop(AP.columns[0], axis=1)
    AP2.iloc[0] = atomic_vals
    AP2 = AP2.astype(float)
    AP3 = AP2.mul(AP2.iloc[0]).iloc[1:]
    total_AP = AP3.sum(axis=1)

    total_phase_shift = 0.65 * (total_AP.reset_index(drop=True) + DP.iloc[:, 1])
    
    return zVals, total_AP, DP, total_phase_shift

def integratorV(zVals, total_phase_shift, radius=200, rEnd=650, rstep=0.5):
    """
    Integrates the phase shift profile to project it onto a spherical vesicle.
    """
    zProfileWidth = np.round(zVals.iloc[-1] - zVals.iloc[0])
    PSPstepSize = round(zVals[2] - zVals[1], 2)

    PSPNorm = np.array(total_phase_shift - total_phase_shift[0])
    PSPVes = np.array(zVals + radius)
    PSPVesicle = np.array([[x, y] for x, y in zip(PSPVes, PSPNorm)])

    duplet_list = np.array([[x, 0] for x in np.arange(0, PSPVes[0], PSPstepSize)])
    duplet_list2 = np.array([[x, 0] for x in np.arange(PSPVes[-1] + PSPstepSize, rEnd + PSPstepSize * 0.01, PSPstepSize)])
    PSP2 = np.concatenate((duplet_list, PSPVesicle, duplet_list2), axis=0)

    x2 = np.array([t[0] for t in PSP2])
    y2 = np.array([t[1] for t in PSP2])
    radial_intensity_function = interp1d(x2, y2, kind='quadratic', fill_value='extrapolate')

    xt = np.concatenate((
        np.arange(0, radius - zProfileWidth, 5 * rstep),
        np.arange(radius - zProfileWidth, radius + zProfileWidth, rstep),
        np.arange(radius + zProfileWidth, rEnd, 5 * rstep)
    ))

    results_list = []
    for n in xt:
        def integrand(z):
            try:
                val = 2 * radial_intensity_function(z) / np.sqrt(1 - (n / z) ** 2)
                return val if np.isfinite(val) else 0
            except (ZeroDivisionError, ValueError):
                return 0

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            result, _ = quad(integrand, n, rEnd, epsabs=1e-6, epsrel=1e-6, limit=100)

        results_list.append((n, result / 600))

    return list(zip(*results_list))