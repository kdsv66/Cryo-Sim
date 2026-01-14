import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os
import tifffile as tiff  # Your original function requires this

def style_axes(ax):
    """Applies a consistent style to a matplotlib axis."""
    ax.axvline(x=0, color='gray', linestyle='--', linewidth=1)
    ax.tick_params(axis='both', which='major', direction='in', top=True, right=True, length=6, width=1, color='gray')
    ax.tick_params(axis='both', which='minor', direction='in', top=True, right=True, length=3, width=0.5, color='gray')
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))

def plot_phase_shifts(zVals, total_AP, DP, total_phase_shift, show_plot=True):
    """Plots the atomic, dipole, and total phase shifts."""
    fig, axs = plt.subplots(3, 1, figsize=(6, 9), sharex=True)
    axs[0].plot(zVals, total_AP, color='tab:blue')
    axs[0].set_ylabel(r'Phase shift [mrad $\mathrm{\AA}^{-1}$]')
    axs[0].set_title('Atomic Phase Shift')
    axs[1].plot(DP.iloc[:, 0], DP.iloc[:, 1], color='tab:green')
    axs[1].set_ylabel(r'Phase shift [mrad $\mathrm{\AA}^{-1}$]')
    axs[1].set_title('Dipole Phase Shift')
    axs[2].plot(zVals, total_phase_shift, color='tab:red')
    axs[2].set_xlabel(r'Distance from bilayer center [$\mathrm{\AA}$]')
    axs[2].set_ylabel(r'Phase shift [mrad $\mathrm{\AA}^{-1}$]')
    axs[2].set_title('Total Phase Shift')
    for ax in axs:
        style_axes(ax)
    plt.tight_layout()
    if show_plot:
        plt.show()
    else:
        plt.close(fig) # Frees up memory if plot isn't shown

def plot_integration_profiles(integrations, vesicle_sizes, show_plot=True):
    """Plots the projected intensity profiles for different vesicle sizes."""
    fig = plt.figure(figsize=(8, 5))
    labels = [f"{w} nm" for w in vesicle_sizes]
    for i, vesicle_data in enumerate(integrations):
        plt.plot(vesicle_data[0], vesicle_data[1], label=labels[i], alpha=0.85)
    plt.title('Projected Intensity Profiles for Vesicle Sizes')
    plt.xlabel('Radial Distance from Center [Å]')
    plt.ylabel('Integrated Phase Intensity')
    plt.xlim(0, vesicle_sizes[-1] * 6)
    ax = plt.gca()
    style_axes(ax)
    plt.legend(title='Vesicle Diameter', loc='center left', bbox_to_anchor=(1.01, 0.5))
    plt.tight_layout()
    if show_plot:
        plt.show()
    else:
        plt.close(fig) # Frees up memory if plot isn't shown

def save_tiff_and_plot(image_data, file_path, title, show_plot=True):
    """
    Saves image data to a TIFF file and optionally displays it in a plot.
    """
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    
    tiff.imwrite(file_path, image_data.astype(np.float32))
    print(f"✅ Image saved to '{file_path}'")
    
    if show_plot:
        plt.imshow(image_data, cmap="gray")
        plt.title(title)
        plt.axis('off')
        plt.show()
