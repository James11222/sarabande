[![codecov](https://codecov.io/gh/James11222/sarabande/branch/main/graph/badge.svg?token=47GPJCFZLE)](https://codecov.io/gh/James11222/sarabande) 
![PyPI](https://img.shields.io/pypi/v/sarabande?color=blue%20&label=PyPi%20)

<p align="center">
  <img src="logo/logo_text_dm.png#gh-dark-mode-only" width="100%">
  <img src="logo/logo_text.png#gh-light-mode-only" width="100%">
</p>

A useful `python` package to measure the 3/4 PCFs of discrete periodic data in NlogN time. This is done using Fast Fourier Transforms. 


## Installation: 
The package is available on PyPi via the command `pip install sarabande`. To check if the code is working properly after pip installation, run 

```python
import sarabande

sarabande.check_install()
```
which will display a message stating if the package was properly installed.

## Basic Usage:

```python
import sarabande

NPCF_obj = sarabande.measure(**kwargs)
sarabande.calc_zeta(NPCF_obj)
zeta = NPCF_obj.zeta
```

Where `**kwargs` can be any of the arguments to the measure constructor function. The possible arguments are:

`Args:`

* `nPCF` ([`int`]): Must be either 3 or 4. Determines how many points we use in our nPCF.
* `projected` ([`bool`]): Flag to determine whether the user wants a projected 3/4 PCF or the Full. Defaults to False.
    - if `projected`:
        - `m_max` ([`int`]): If user chooses projected, we set an m_max (similar to the `ell_max` in 3D)
    - if not `projected`:
        - `ell_max` ([`int`]): If user choosees not projected (full nPCF) then ell_max is the highest order for calculation.

* `density_field_data` ([`ndarray`]): A square ndarray of data that is periodic. Must be 2D for projected and 3D for full.
* `save_dir` ([`string`]): A string to tell the algorithm where to save and store files. All temporary files will be stored here.
* `save_name` ([`string`]): A string to tell the algorithm what to name the files.
* `nbins` ([`int`]): Number of bins to be used in nPCF calculation.
* `bin_spacing` ([`string`]): A string to determine the spacing of bins. Options are `'LIN'`, `'INV'`, or `'LOG'`
* `bin_min` ([`int`]): The lower bound of the inner most bin. Default is 1. Optional.
* `physical_boxsize` ([`float`]): An optional parameter if using a physical scale. The length of one side of the data.
* `rmin` ([`float`]): minimum calculation distance (determins `bin_min`)
* `rmax` ([`float`]): maximum calculation distance (determins `bin_max`)
* `normalize` ([`bool`]): A boolean flag to normalize the 3/4 PCFs. Defaults to True. Can't use normalize without giving a `physical_boxsize`, `rmin`, and `rmax` first.
* `particles_on_grid` ([`bool`]): An optional boolean flag to modify the normalization scheme slightly. This is recommended if you are working with particles on the grid mesh where a given cell corresponds to a particle. 

We note that the `calc_zeta` method has an optional boolean argument `verbose_flag` which can be turned on and off depending on if the user wants to see the steps of the code printed. We also add an optional boolean argument `parallelized` which can be turned on and off if the user wishes to compute the Full 4PCF serially. This is added due to the instability of `concurrent.futures` and parallel processing in python across different machines. 

For an example, please visit the demo notebook in the analysis notebooks folder: `notebooks/Application_Example.ipynb`

 ## Workflow:    
The map of SARABANDE is as follows:

<p align="center">
  <img src="notebooks/paper_figures/workflow_dm.png#gh-dark-mode-only" width="100%">
  <img src="notebooks/paper_figures/workflow.png#gh-light-mode-only" width="100%">
</p>

For more information about each algorithm, please read [Sunseri et al. 2022](https://ui.adsabs.harvard.edu/abs/2022arXiv221010206S/abstract)

## Coverage [![codecov](https://codecov.io/gh/James11222/sarabande/branch/main/graph/badge.svg?token=47GPJCFZLE)](https://codecov.io/gh/James11222/sarabande) 
We provide a sunburst plot of the code coverage for sarabande below provided by codecov.io. The inner-most circle is the entire project, moving away from the center are folders then, finally, a single file. The size and color of each slice is representing the number of statements and the coverage, respectively.

<p align="center">
  <img src="https://codecov.io/gh/James11222/sarabande/branch/main/graphs/sunburst.svg?token=47GPJCFZLE">
</p>
