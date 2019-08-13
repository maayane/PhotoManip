# PhotoManip
`PhotoManip` is a package to manipulate, convert and model SNe Photometry.

[![PyPI](https://img.shields.io/pypi/v/PhotoManip.svg?style=flat-square)](https://pypi.python.org/pypi/PhotoManip)

## Documentation

## Credit

## How to install `PhotoManip`?

### pip

`pip install PhotoManip`

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

## How to run `PhotoManip`?
Before running anything, you first need to define the path of a directory where all your transmission curves are stored. We provide such a directory with the package.
```python
>>> import PhotoManip
>>> from PhotoManip import PhotoManip_fun
>>> filters_directory='./Filters' 
```

### Convert a mag (AB) value in a given flux into a flux value.

```python
>>> mag_value=20.55 #this is the magAB  value to convert
>>> flux=PhotoManip_fun.magAB_in_filter_to_flux_in_filter(mag_value,Filter_vector=np.array([['swift','UVW2']]),filters_directory=filters_directory,verbose=False)
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
PhotoManip_fun.read_data_Marshall_simple(data_path,no99=True,filters_directory=filters_directory,output_path=output_path)
```

### Calculate the time "tref" at which the flux in a given band crosses zero

The explosion date is often approximated by the extrapolated time at which the flux (e.g. in the r-band) crosses zero. `PhotoManip` allows you to calculate this time by modeling the light curve with
1. a concave exponent
2. a power law
and then chosing the best model.

First, specify the band you would like to use (r or g), the path to the data (the data should be in the format 'jd','mag','magerr','flux','fluxerr','absmag','absmagerr','filter','instr'; see previous section how to produce such format from Marshall-downloaded data) and the number of points you would like to include in the fit:
```python
band='r'
days_rising=10
```

Then run
```python
tref=PhotoManip_fun.tref_from_P48('./data_formated_aatqzim.txt',band=band,days_rising=days_rising)[0]
```

The code will first show you the light curve and nmber of points to include in the fit, and warn you to adjust the `days_rising` parameter.
<p align="center">
  <img src="./test/result_fit_sed_mat/day_1.359/SED_date_1.359.png" width="350">
</p>

It will the run mcmc to model the rising chunk with both a concave exponent and a power law, output plots and summarize the results for you:


All plots can be found in the `` directory.
 

## Give it a try with the test data!

All the figures above were obtained by running `PhotoManip` on the multiple-bands light curve of the Supernova PTF13dqy ([Yaron et al 2017](https://ui.adsabs.harvard.edu/#abs/2017NatPh..13..510Y/abstract)). The data is available in the `test` directory (including the output of the time-consuming interpolation step). You can reproduce all these results and figures by running `PhotoManip` with the parameters file `params_test.py` as it is.


