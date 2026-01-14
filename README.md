# Cryo-Sim
Cryo-SIM converts molecular dynamics data—number density and dipole potentials—into synthetic cryo-EM images of liposomes. Users can customize voltage, defocus, resolution, and noise to mimic real microscopy. It outputs noisy and noise-free images plus radial profiles, helping researchers validate simulations against experimental data.
# This project uses Conda to manage its environment and dependencies. You must have Anaconda or Miniconda installed on your system.
Download the cryo-sim package to the directory of your choice.

In terminal
# Make sure you are in the same directory as your environment.yml and cryo_sim.ipynb file
conda env create -f environment.yml
conda activate cryo-sim
jupyter notebook

Jupyter notebook will open in browser
open cryo_sim.ipynb



/your-project-directory
|
|-- README.md                     # This helper file
|-- environment.yml               # Conda environment definition file
|
|-- cryo_sim.ipynb        	  # Main Jupyter Notebook for running the workflow
|
|-- /scripts/
|-- core_functions.py             # Functions for loading data and integration
|-- simulation.py                 # Contains the VesicleSimulator class
|-- image_processing.py           # Functions for radial profiling and centering
|-- analysis.py                   # Final analysis function for profiles
|-- plotting.py                   # plotting and visualization functions
|
|-- ndp_cut.dat                   # INPUT: Atomic potential data (required)
|-- dipole.dat                    # INPUT: Dipole potential data (required)
|
|-- /noisy_images/                # OUTPUT: where the final images with noise are produced and also directory for analysis
|-- /related_images/              # OUTPUT: Diagnostic TIFF images
|-- /noise_free/                  # OUTPUT: Final noise-free TIFF images
|
|-- /results/
|-- centered_profiles_wide.csv    # OUTPUT: Processed radial profiles
|-- profile_analysis_result.csv   # OUTPUT: Final analysis results
|-- bilayer_projection.json       # OUTPUT: integrated projection file name, change name if generating new file or keep same to use previously genrated file


cryo_sim.ipynb contains comments as well for help

changeable parameters
ANGSTROM_PER_PIXEL = 2.7
VESICLE_SIZES = np.arange(30, 130, 10) # (start diameter, end diameter - interval, interval) in nm
defocus = 2 # in μm
accelarating_voltage = 300 # in kV
ADD_STRUCTURAL_NOISE = False  # Set to False to disable
STRUCTURAL_NOISE_STD = 0.1
NOISY_NUM = 1 # number of noisy images to generate for each size
Scaling_Factor_img = 2.5 # control pixel intensity
Stddev_Noise_img = 200 # control noise deviation
intensity_profile_smoothing = 3 # control smoothing of intensity profile
For hiding additional plots
SHOW_PLOTS = False # will hide most of the additional plots in the notebook
