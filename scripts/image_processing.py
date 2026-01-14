# image_processing.py

import os
import re
import numpy as np
from tifffile import imread
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d  # <<< THIS LINE IS THE FIX

def process_image(filename, image_directory, angstrom_per_pixel):
    """Extracts metadata and image data for a single image file."""
    img_files = sorted([f for f in os.listdir(image_directory) if f.endswith(".tiff")])
    vesicle_diameter = int(re.search(r'(\d+)\s*nm', filename).group(1))
    pixel_radius = vesicle_diameter * 10 / (angstrom_per_pixel * 2)

    img = imread(os.path.join(image_directory, filename))
    h, w = img.shape
    cx, cy = w / 2, h / 2

    offset = 75 / angstrom_per_pixel
    return img, cx, cy, pixel_radius - offset, pixel_radius + offset, pixel_radius

def radial_intensity_profile(img, cx, cy, r1, r2, angstrom_per_pixel, smooth_sigma=2):
    """Calculates the radial intensity profile of an image within an annulus."""
    y, x = np.indices(img.shape)
    r = np.hypot(x - cx, y - cy)
    mask = (r >= r1) & (r <= r2)
    
    r_vals = r[mask]
    intensities = img[mask]
    
    step = 1 / angstrom_per_pixel
    bins = np.arange(r1, r2 + step, step)
    bin_idx = np.digitize(r_vals, bins)

    radial_means = np.array([
        intensities[bin_idx == i].mean()
        for i in range(1, len(bins))
        if np.any(bin_idx == i)
    ])
    radial_r = np.array([
        (bins[i] + bins[i - 1]) / 2
        for i in range(1, len(bins))
        if np.any(bin_idx == i)
    ])

    if smooth_sigma > 0:
        radial_means = gaussian_filter1d(radial_means, sigma=smooth_sigma)

    return radial_r, radial_means

def center_and_pad_profile(radii_nm, intensities, target_range=(-75, 75), step=1.0):
    """Centers a profile on its local maximum and pads it to a standard range."""
    center_idx = np.argmin(np.abs(radii_nm))
    search_range = 10
    start = max(0, center_idx - search_range)
    end = min(len(intensities), center_idx + search_range + 1)

    # Handle cases where the profile is empty or too short
    if len(intensities[start:end]) == 0:
        return np.array([]), np.array([])

    local_max_idx = start + np.argmax(intensities[start:end])
    x_shift = radii_nm[local_max_idx]
    y_shift = intensities[local_max_idx]

    shifted_radii = radii_nm - x_shift
    shifted_intensities = intensities - y_shift

    mask = (shifted_radii >= target_range[0]) & (shifted_radii <= target_range[1])
    shifted_radii = shifted_radii[mask]
    shifted_intensities = shifted_intensities[mask]
    
    # Handle cases where the mask results in an empty array
    if len(shifted_radii) == 0:
        return np.array([]), np.array([])
    
    padded_radii = np.arange(target_range[0], target_range[1] + step, step)
    interp_func = interp1d(shifted_radii, shifted_intensities, kind='linear', fill_value=(shifted_intensities[0], shifted_intensities[-1]), bounds_error=False)
    padded_intensities = interp_func(padded_radii)

    return padded_radii, padded_intensities