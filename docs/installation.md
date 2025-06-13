---
title: Installation and Dependencies
sidebar_position: 3
---

It is highly recommendable that `TidyScreen` is installed in an separate Conda virtual environment to prevent inconsistencies among the libraries currently installed in your system and those specifically required by the package itself:

```python
# Create an isolated environment with the required libraries
$ conda create -n tidyscreen python=3.10 chemicalite

# Activate the corresponding environment
$ conda activate tidyscreen

# Install TidyScreen from the corresponding GitHub repository
$ pip install git+https://github.com/alfredoq/TidyScreen_v2
```

As part of compounds parameterization, `TidyScreen` implements a charge assignment method by taking profit of the accuracy and speed of a neural network-based approximation ([*Wang, Y; et.al.*](https://pubs.acs.org/doi/10.1021/acs.jpca.4c01287)). As will be presented in the corresponding section, other charge system can be used, however it this so-called *bcc-ml* is highly recommended. In order to be able to use *bcc-ml* fitted charges, an additional library should be installed as follows in the `TidyScreen` environment:

```python
# Install EspalomaCharge library
$ conda install -c conda-forge espaloma_charge openff-toolkit
```

In addition, in order to prepare molecules for docking studies `TidyScreen` workflows uses library called [Meeko](https://github.com/forlilab/Meeko), which was developed by the [ForliLab](https://forlilab.org/) at Scripps Research Institute.Currently, `TidyScreen` requires a specific development version of [Meeko](https://github.com/forlilab/Meeko) capable of reading .mol2 atomic charges. This development version should be installed in the corresponding `TidyScreen` enviroment as follows:


```python
# Install Meeko library
$ pip install git+https://github.com/forlilab/Meeko@develop 
