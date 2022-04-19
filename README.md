# PEPICO_Elettra_VG_python
Python code for analyzing electron-ion coincidence experiments conducted at the GasPhase beamline of Elettra Sincrotrone Trieste using the VG electron analyzer.

Coincidence events of the raw data files corresponding to the same kinetic energy are analyzed with "main.py" which uses the four classes "VG_tof", "VG_Tof_only", "VG_Movie_Tof", and "VG_figures".

The final PEPICO matrix can then be obtained using the script "PEPICO_complete_sp.py".

A mini-manual is provided in "PEPICO_Elettra_VG.pdf"

In order to guarantee package compatibilities, it is recommended to work in a virtual environment that can be created from the PEPICO_Elettra_VG_env.yml file:
If you use Anaconda you can create the environment using the terminal or Anaconda Prompt:
1. "conda env create -f PEPICO_Elettra_VG_env.yml"
2. activate the new environment with "conda activate PEPICO_Elettra_VG_env"
