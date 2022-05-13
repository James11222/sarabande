[![codecov](https://codecov.io/gh/James11222/sarabande/branch/main/graph/badge.svg?token=47GPJCFZLE)](https://codecov.io/gh/James11222/sarabande) 
![PyPI](https://img.shields.io/pypi/v/sarabande?color=blue%20&label=PyPi%20)

<img src="https://github.com/James11222/sarabande/blob/main/Images/logo_text.png?raw=true" /> 

A useful `python` package to measure the 3/4 PCFs of discrete periodic data in NlogN time. This is done using Fast Fourier Transforms.


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
* `nbins` ([int]): Number of bins to be used in nPCF calculation.
* `bin_spacing` ([`string`]): A string to determine the spacing of bins. Options are `'LIN'`, `'INV'`, or `'LOG'`
* `bin_min` ([`int`]): The lower bound of the inner most bin. Default is 1. Optional.
* `physical_boxsize` ([`float`]): An optional parameter if using a physical scale. The length of one side of the data.
* `rmin` ([`float`]): minimum calculation distance (determins `bin_min`)
* `rmax` ([`float`]): maximum calculation distance (determins `bin_max`)

 ## Workflow:    
The map of SARABANDE is as follows:

<img src="https://github.com/James11222/sarabande/blob/main/Images/workflow.png?raw=true" width="100%" />

For more information about each algorithm, please read (Sunseri et al. 2022 in prep)
