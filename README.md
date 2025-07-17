# mos2-defects

Data and analysis in support of "Defect density quantification in monolayer MoS2 using helium atom micro-diffraction" by A Radic, NA von Jeinsen, K Wang, Y Zhu, I Sami, V Perez, DJ Ward, AP Jardine, M Chhowalla, SM Lambrick.

A pre-print of the manuscript can be found at arXiv: [https://doi.org/10.48550/arXiv.2409.18637](https://doi.org/10.48550/arXiv.2409.18637). 

The plotting and analysis scripts make use of the [SHeM diffraction analysis](https://github.com/slambrick/SHeM-diffraction-analysis) package. This package is expected to be in the directory the scripts are called. This can be achieved with 
```bash
git clone https://github.com/slambrick/SHeM-diffraction-analysis.git
```
(or downloading and extracting the `.zip`) then renaming the folder as `SHeM_diffraction_analysis` as python doesn't like `-` in directory names.

## Scripts

`sample1a_lattice_constant.ipynb` calculates the lattice constant of the pristine monolayer and produces plots to demonstrate the procedure. The data is plotted in figure S1 of the manuscript.

`samples1a_3_figures.ipynb` creates two 2D plots of the reciprocal lattice for a high defect sample and a 'pristine' sample. These are plotted in both figure 1 (small) and 3 (full sized).

`Defect_rate_comparison.ipynp` analyses the diffraction intensity of the samples with different defect density. It results in plots of the diffraction intensity as a function of the defect rate at 2 different samples.

`substrate_comparison.ipynb` compares "z-scans" from the substrate, from bulk MoS2 and monolayer MoS2 on two different substrates (directly on SiO2 or on hBN-on-SiO2).

`plot_fits.ipynb` plots the fits to the lattice gas model used to generate  figures S3, S4 and 4.

## Additional uses of data

Some of the data and analysis in this repository is also used in another publication: "Helium atom micro-diffraction as a
characterization tool for 2D materials". Specifically `samples1a_3_figures.ipynb` creates the 2D plot with identified lattice positions for figure 4. A [preprint of this manuscript may be found on arXiv](https://doi.org/10.48550/arXiv.2409.20461).

## Images

The directory "figure_2" contains the SHeM and optical images used for figure 2. Optical is provided as a `.png` file. The SHeM image is provided in a number of formats, first the `.png`, as well as the raw data as `.mat` and `.txt` files. Plus a [Gwyddion](https://gwyddion.net/) file for easy image processing.

## Supplementary information figures

The input potential files for [MultiScat](https://github.com/Cambridge-Atom-Scattering-Centre/Multiscat) ([Method](https://doi.org/10.1039/FT9908601641)) representing the potential used to generate the results in figure 6, and displayed in figure 5 is given in the "figure_S5" folder. 

The diffraction amplitudes resulting from the MultiScat simulation, which are displayed in figure S6 are in the "figure_S6" folder.

## Raw data

Raw data used in the various scripts to produce the figures is in the `data` folder. In general SHeM data is stored in `.mat` files.