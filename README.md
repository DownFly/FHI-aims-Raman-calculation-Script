# FHI-aims Raman spectra calculation scripts

We used FHI-aims to calculate the Raman spectra data of 3504 two-dimensional structures. This repository stores the scripts we used in the calculation process.

### plot.py

A script for plotting Raman spectra.

### pos2geom.py

Convert POSCAR to geometry.in

### check.sh

Check if the structure is calculated successfully

### 01_prepare_calc.sh & 01_calc_aims_harmonic_raman.pl

Prepare for calculation

### 02_submit_calc.sh

Batch submit calculations, resubmit if errors are detected

### 03_read_calc.sh & 02_read_aims_harmonic_raman.pl

Read the calculation results
