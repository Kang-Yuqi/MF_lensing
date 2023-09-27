from setuptools import setup, find_packages

setup(
    name="MFLensingTools",
    version="0.1",
    author="Yuqi Kang",
    packages=find_packages(),
    python_requires='>3.0.0',
    install_requires=[
        "numpy>=1.17.0",
        "healpy>=1.16.0",
        "matplotlib>=3.5.0",
        "mtneedlet"
    ]
)