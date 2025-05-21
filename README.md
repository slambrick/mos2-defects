# mos2-defects

Data and analysis in support of "Defect density quantification in monolayer MoS2 using helium atom micro-diffraction" by A Radic, NA von Jeinsen, K Wang, Y Zhu, I Sami, V Perez, DJ Ward, AP Jardine, M Chhowalla, SM Lambrick.

A pre-print of the manuscript can be found at arXiv: [https://doi.org/10.48550/arXiv.2409.18637](https://doi.org/10.48550/arXiv.2409.18637). 

The plotting and analysis scripts make use of the [SHeM diffraction analysis](https://github.com/slambrick/SHeM-diffraction-analysis) package. This package is expected to be in the directory the scripts are called. This can be achieved with 
```bash
git clone https://github.com/slambrick/SHeM-diffraction-analysis.git
```
(or downloading and extracting the `.zip`) then renaming the folder as `SHeM_diffraction_analysis` as python doesn't like `-` in directory names.

## Scripts

`sample1a_lattice_constant.ipynb` calculates the lattice constant of the pristine monolayer and produces plots to demonstrate the procedure.

`samples1a_3_figures.ipynb` creates two 2D plots of the reciprocal lattice for a high defect sample and a 'pristine' sample.