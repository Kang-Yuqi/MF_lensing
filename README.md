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

## Citing
If you use MF_lensing in your research, please cite it as follows:

```
@article{Hamann:2023tdu,
    author = "Hamann, Jan and Kang, Yuqi",
    title = "{A Minkowski functional analysis of the Cosmic Microwave Background weak lensing convergence}",
    eprint = "2310.14618",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1088/1475-7516/2024/05/076",
    journal = "JCAP",
    volume = "05",
    pages = "076",
    year = "2024"
}
```
## Developer
* Yuqi Kang
