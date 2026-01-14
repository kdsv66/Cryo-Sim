# simulation.py

import numpy as np
import math
from scipy.interpolate import interp1d

class VesicleSimulator:
    """
    A class to encapsulate the vesicle image simulation process.
    """
    # 1. UPDATED INIT TO ACCEPT STRUCTURAL NOISE PARAMETERS
    def __init__(self, scaled_data_entry, params, angstrom_per_pixel, add_structural_noise=False, structural_noise_std=0.1):
        self.params = params
        self.scaled_data_entry = scaled_data_entry
        self.angstrom_per_pixel = angstrom_per_pixel
        # Store the structural noise settings
        self.add_structural_noise = add_structural_noise
        self.structural_noise_std = structural_noise_std 
        
        self.lambda_angstrom = self._electron_wavelength(params['kev'])
        self.projection_pixels, self.ft_matrix, self.ctf_matrix = None, None, None
        self.xco, self.yco = None, None

    def _electron_wavelength(self, kev):
        h, m_e, c, e = 6.626e-34, 9.109e-31, 2.997e8, 1.602e-19
        V = kev * 1e3
        return (h / np.sqrt(2 * m_e * e * V * (1 + (e * V) / (2 * m_e * c**2)))) * 1e10

    def _ctf_func(self, svec):
        delta_z, Q, B = self.params['delta_z'], self.params['Q'], self.params['B']
        s_norm_sq = np.linalg.norm(svec) ** 2
        chi_val = -np.pi * self.lambda_angstrom * delta_z * 1e4 * s_norm_sq
        return (np.sin(chi_val) - Q * np.cos(chi_val)) * np.exp(-B * s_norm_sq)

    def _calculate_projection(self):
        if self.projection_pixels is None:
            x_vals_scaled = np.array(self.scaled_data_entry[0])
            y_vals = np.array(self.scaled_data_entry[1])
            rEnd = x_vals_scaled[-1] * self.angstrom_per_pixel + 1.5
            pixels = math.ceil(2 * rEnd / self.angstrom_per_pixel)
            self.xco = self.yco = np.arange(-pixels / 2, pixels / 2 + 1)
            xx, yy = np.meshgrid(self.xco, self.yco)
            rValCo = np.sqrt(xx**2 + yy**2)
            intensityFunc = interp1d(x_vals_scaled, y_vals, kind='quadratic', fill_value='extrapolate')
            
            # Calculate base projection
            self.projection_pixels = intensityFunc(rValCo)

            # 2. APPLY STRUCTURAL NOISE IF ENABLED
            if self.add_structural_noise:
                noise = np.random.normal(0, self.structural_noise_std, self.projection_pixels.shape)
                self.projection_pixels = self.projection_pixels + noise

    def _calculate_fft(self):
        if self.ft_matrix is None:
            self._calculate_projection()
            self.ft_matrix = np.fft.fft2(self.projection_pixels)

    def _calculate_ctf_matrix(self):
        if self.ctf_matrix is None:
            self._calculate_projection()
            y_dim, x_dim = self.projection_pixels.shape
            s_coord_x = (np.floor(x_dim / 2) - np.abs(np.arange(-(x_dim - 1) / 2, (x_dim - 1) / 2 + 1))) / (self.angstrom_per_pixel * x_dim)
            s_coord_y = (np.floor(y_dim / 2) - np.abs(np.arange(-(y_dim - 1) / 2, (y_dim - 1) / 2 + 1))) / (self.angstrom_per_pixel * y_dim)
            s_coords_mat_x, s_coords_mat_y = np.meshgrid(s_coord_x, s_coord_y)
            s_coords = np.stack((s_coords_mat_x, s_coords_mat_y), axis=-1)
            self.ctf_matrix = np.array([[self._ctf_func(svec) for svec in row] for row in s_coords])
            
    def get_projection_image(self):
        self._calculate_projection()
        return self.projection_pixels

    def get_fft_image(self):
        self._calculate_fft()
        return np.abs(np.fft.fftshift(self.ft_matrix))

    def get_ctf_image(self):
        self._calculate_ctf_matrix()
        return np.abs(np.fft.fftshift(self.ctf_matrix))

    def get_noisefree_image(self, scaling_factor=2.5):
        """Generates the final image without any added noise."""
        self._calculate_fft()
        self._calculate_ctf_matrix()
        I0, m = self.params['I0'], self.params['m']
        
        img_conv_mat = self.ft_matrix * self.ctf_matrix
        IS_mat = 1 + m * img_conv_mat
        IS_abs_mat, IS_arg_mat = np.abs(IS_mat), np.angle(IS_mat)
        img_conv_ivt_mat_t = I0 + np.real(np.fft.ifft2(IS_abs_mat * np.exp(1j * IS_arg_mat)))
        img_conv_ivt_mat = img_conv_ivt_mat_t.copy()
        img_conv_ivt_mat[0, 0] = I0
        
        return scaling_factor * (img_conv_ivt_mat - np.mean(img_conv_ivt_mat))

    def get_noisy_image(self, scaling_factor=2.5, stddev=217):
        """Generates the final image with added Gaussian noise and offset."""
        noise_free_component = self.get_noisefree_image(scaling_factor)
        
        random_offset = np.random.uniform(-10, 10)
        gaussian_noise = np.random.normal(0, stddev, noise_free_component.shape)
        
        return noise_free_component + random_offset + gaussian_noise