Cryo-Sim
Cryo-SIM converts molecular dynamics data—number density and dipole potentials—into synthetic cryo-EM images of liposomes. Users can customize voltage, defocus, resolution, and noise to mimic real microscopy. It outputs noisy and noise-free images plus radial profiles, helping researchers validate simulations against experimental data.


\# This project uses Conda to manage its environment and dependencies. You must have Anaconda or Miniconda installed on your system.

Download the cryo-sim package to the directory of your choice.

**In terminal**

\# Make sure you are in the same directory as your **environment.yml and cryo\_sim.ipynb** file

conda env create -f **environment.yml**

conda activate **cryo-sim**

jupyter notebook



**Jupyter notebook will open in browser**

open cryo\_sim.ipynb







/your-project-directory

|

|-- README.md                     # This helper file

|-- environment.yml               # Conda environment definition file

|

|-- cryo\_sim.ipynb        	  # Main Jupyter Notebook for running the workflow

|

|-- /scripts/

|-- core\_functions.py             # Functions for loading data and integration

|-- simulation.py                 # Contains the VesicleSimulator class

|-- image\_processing.py           # Functions for radial profiling and centering

|-- analysis.py                   # Final analysis function for profiles

|-- plotting.py                   # plotting and visualization functions

|

|-- ndp\_cut.dat                   # INPUT: Atomic potential data (required)

|-- dipole.dat                    # INPUT: Dipole potential data (required)

|

|-- /noisy\_images/                # OUTPUT: where the final images with noise are produced and also directory for analysis

|-- /related\_images/              # OUTPUT: Diagnostic TIFF images

|-- /noise\_free/                  # OUTPUT: Final noise-free TIFF images

|

|-- /results/

|-- centered\_profiles\_wide.csv    # OUTPUT: Processed radial profiles

|-- profile\_analysis\_result.csv   # OUTPUT: Final analysis results

|-- bilayer\_projection.json       # OUTPUT: integrated projection file name, change name if generating new file or keep same to use previously genrated file





**cryo\_sim.ipynb** contains comments as well for help



changeable parameters

ANGSTROM\_PER\_PIXEL = 2.7

VESICLE\_SIZES = np.arange(30, 130, 10) # (start diameter, end diameter - interval, interval) in nm

defocus = 2 # in μm

accelarating\_voltage = 300 # in kV

ADD\_STRUCTURAL\_NOISE = False  # Set to False to disable

STRUCTURAL\_NOISE\_STD = 0.1

NOISY\_NUM = 1 # number of noisy images to generate for each size

Scaling\_Factor\_img = 2.5 # control pixel intensity

Stddev\_Noise\_img = 200 # control noise deviation

intensity\_profile\_smoothing = 3 # control smoothing of intensity profile

For hiding additional plots

SHOW\_PLOTS = False # will hide most of the additional plots in the notebook
