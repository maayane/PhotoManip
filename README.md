# Photomanip
This package is for doing all kind of manipulations conversions, and modeling of photometric data. E.g. calculating the time tref at which the flux in a given band crosses zero (this time if ofter taken as the explosion date).

[![PyPI](https://img.shields.io/pypi/v/Photomanip.svg?style=flat-square)](https://pypi.python.org/pypi/Photomanip)

## Documenation

`Photomanip` is a package to manipulate, convert and model SNe Photometry.

## Credit

## How to install `Photomanip`?

### pip

`pip install Photomanip`

### Python version
* `python 3`

### Required python packages
* `math`
* `numpy`
* `scipy`
* `pylab`
* `emcee`
* `pylab`
* `pandas`
* `csv`

## How to run `Photomanip`?
Before running anything, you first need to define the path of a directory where all your transmission curves are stored. We provide such a directory with the package.
```python
>>> import Photomanip
>>> from Photomanip import Photomanip_fun
>>> filters_directory='./Filters' 
```

### Convert a mag (AB) value in a given flux into a flux value.

```python
>>> mag_value=20.55 #this is the magAB  value to convert
>>> flux=Photomanip_fun.magAB_in_filter_to_flux_in_filter(mag_value,Filter_vector=np.array([['swift','UVW2']]),filters_directory=filters_directory,verbose=False)
>>> print('the flux is',flux)
flux is 1.55285808219e-16
```

### Convert data downloaded from the Marshall into the format 'jd','mag','magerr','flux','fluxerr','absmag','absmagerr','filter','instr'.

First define the path of the photometry file downloaded from the Marshall
```python
data_path='./data_test.txt'# a file downloaded from the Marshall
```
Then define the path of the output file

```python
output_path='./data_formated.txt' #path to the file with the right format
```
Then run
```python
Photomanip_fun.read_data_Marshall_simple(data_path,no99=True,filters_directory=filters_directory,output_path=output_path)
```

<p align="center">
  <img src="./test/result_fit_sed_mat/day_1.359/SED_date_1.359.png" width="350">
  <img src="./test/result_fit_sed_mat/day_1.359/fit_result_FLux.png" width="350">
</p>

## Give it a try with the test data!

All the figures above were obtained by running `Photomanip` on the multiple-bands light curve of the Supernova PTF13dqy ([Yaron et al 2017](https://ui.adsabs.harvard.edu/#abs/2017NatPh..13..510Y/abstract)). The data is available in the `test` directory (including the output of the time-consuming interpolation step). You can reproduce all these results and figures by running `Photomanip` with the parameters file `params_test.py` as it is.


