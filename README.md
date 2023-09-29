# MF_lensing
MF_lensing provides a python package for analyzing full-sky maps with Minkowski Functionals(MF). 

## Characteristics
* Measure **MFs** from masked maps or alm.
* Measure **non-Gaussianity parameters** (Skewness & Kurtosis) either locally or globally.
* Analysis maps with **needlet filters**.
* Generate corrections for MFs using non-Gaussianity parameters measured in weak non-Gaussian maps, based on the MFs differences between Gaussian and non-Gaussian maps.

## Installation
To install MF_lensing, navigate to the directory where you downloaded the package and install it using `pip`:
```python
cd MF_lensing
pip install .
```

## Usage
Several examples are provided in the `example` folder.
To use the function `Planck_data.py`, which is based on Planck data, you should first download the related Planck data. It can be obtained, for instance, from [Planck Legacy Archive](https://pla.esac.esa.int/#home).

## Developer
* Yuqi Kang
