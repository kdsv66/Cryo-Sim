# analysis.py

import pandas as pd
import numpy as np
import re
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema

def analyze_profiles_from_csv(profile_csv_path, interp_points=5000):
    """
    Performs a detailed analysis of centered profiles from a CSV file.
    """
    df_profiles = pd.read_csv(profile_csv_path)
    distances = df_profiles.iloc[:, 0].to_numpy()
    profile_columns = df_profiles.columns[1:]
    analysis_results = []

    for obj_id in profile_columns:
        profile = df_profiles[obj_id].to_numpy()
        match = re.match(r"(\d{2,3})nm", obj_id)
        size_nm = int(match.group(1)) if match else np.nan

        # Min/Max analysis
        left_mask = distances < 0
        right_mask = distances > 0
        if not (np.any(left_mask) and np.any(right_mask)):
            continue
        
        ymin_left = np.min(profile[left_mask])
        ymin_right = np.min(profile[right_mask])
        ratio_minima = ymin_right / ymin_left if ymin_left != 0 else np.nan

        # Interpolated profile for finding extrema
        interp_func = interp1d(distances, profile, kind='cubic', fill_value="extrapolate")
        fine_radii = np.linspace(distances.min(), distances.max(), interp_points)
        fine_intensities = interp_func(fine_radii)

        minima_indices = argrelextrema(fine_intensities, np.less)[0]
        minima_radii = fine_radii[minima_indices]
        minima_values = fine_intensities[minima_indices]

        # Find key minima
        x_min_left, x_min_right = np.nan, np.nan
        if np.any(minima_radii < 0):
            left_radii = minima_radii[minima_radii < 0]
            left_values = minima_values[minima_radii < 0]
            x_min_left = left_radii[np.argmin(left_values)]
        if np.any(minima_radii > 0):
            right_radii = minima_radii[minima_radii > 0]
            right_values = minima_values[minima_radii > 0]
            x_min_right = right_radii[np.argmin(right_values)]

        min_distance = abs(x_min_right - x_min_left) if not (np.isnan(x_min_left) or np.isnan(x_min_right)) else np.nan

        analysis_results.append({
            "filename": obj_id, "size_nm": size_nm,
            "ymin_left": ymin_left, "ymin_right": ymin_right, "ratio_minima": ratio_minima,
            "x_min_left": x_min_left, "x_min_right": x_min_right, "min_distance": min_distance
        })

    return pd.DataFrame(analysis_results)