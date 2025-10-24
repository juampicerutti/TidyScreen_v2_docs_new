---
title: Molecules preparation
---

The preparation (i.e. parameterization) of molecules to be studied as candidates constitutes a critical step in the design and execution of of screening campaign. The selection of the parameterization procedures represents an important decision within the design and development stage of the screenign workflow.

After being read from a .csv file (or prepared within TidyScreen through virtual synthesis), and prior to performing molecular docking, molecules requires a preparation step within TidyScreen using the following function:

```python
> # Import required modules
>>> from tidyscreen import tidyscreen as ts
>>> from tidyscreen.chemspace import chemspace as chemspace

# A project named 'sample_project' is assumed to contain all the information
#### Activate the project
>>> sample_project = ts.ActivateProject("sample_project")

#### Activate when using chemspace functions!
>>> sample_project_cs = cs.ChemSpace(sample_project)

# Here we assume the a table name 'test_table' is available
# in the 'chemspace.db' and that it contains the corresponding
# molecules in the SMILES notation.

>>> sample_project_cs.generate_mols_in_table("test_table")
```
Upon executing this function, now `test_table` will contain additional columns identified as:

- `pdb_file`: the *.pdb* file generated from the corresponding SMILES input.
- `mol2_file_sybyl`: a *.mol2* file generated from the *.pdb* counterpart with SYBYL atom typing
- `mol2_file_gaff`: a *.mol2* file generated from the *.pdb* counterpart with GAFF atom typing
- `frcmod_file`: a *.frcmod* used for Amber like computations
- `pdbqt_file`: the *.pdbqt* file that will be used for docking runs
- `charge_model`: a string describing the `charge` and `pdbqt` method used to parameterize the *.pdbqt* file. Read below for a detailed explanation on this. 


Several methods for the computation of molecular parameters have been implemented in TidyScreen to compute charges and prepare `.pdbqt` files intended for molecular docking. These methods as basically controlled by two keywords:

- `charge_method`: indication the type of charge method to be used to assign atomic charges to the molecules. Three methods are available:
    - `"gas"`: corresponds to *gasteiger* atomic charges
    - `"bcc"`: stands for *"Bond Charge Correction"*, which are semi-empirical atomic charges as calculated by *Antechamber*. The user needs to consider that computation of this type of charge is takes a long time, so it should not be used to prepare large libraries of molecules.
    - `"bcc-ml"`: this is a machine learning-accelerated version of BCC (Bond Charge Correction) charges.
- `pdbqt_method`: refers to the methodology to be used to compute the above mentioned charge. Three methods are available within TidyScreen:
  - `"prepare_ligand_script"`: this method will use a script shipped with the AutodockTools package to prepare the '.pdbqt' file. This method always uses the `'gas'` charge method.
  - `meeko`: Will use Meeko capabilities to assign charges. Both `'gas'` and `'bcc-ml'` are accepted. If the later is select, the [espaloma](https://docs.espaloma.org/en/latest/) package will be used.
  - `third_party`: will use Antechamber to compute charges taking into account the indicated configuration. Both `'gas'` and `'bcc-ml'` are accepted. If the later is select, the [espaloma_charge](https://pypi.org/project/espaloma-charge/) package will be used.

:::tip Default method
If only the table name is passed to `generate_mols_in_table()`, the default configuration will be applied which is `charge_method = 'gas'` and `pdbqt_method = 'prepare_ligand_script'`
:::

Below the user can find some custom examples of the usage of the `generate_mols_in_table("test_table")` function:


```python
# Will use Meeko to compute gasteiger charges
sample_project_cs.generate_mols_in_table("test_table", charge_method='gas', pdbqt_method="meeko")

# Will use Meeko to compute espaloma type charges
sample_project_cs.generate_mols_in_table("test_table", charge_method='bcc-ml', pdbqt_method="meeko")

# Will use the espaloma_charge package to compute bcc-ml charges
sample_project_cs.generate_mols_in_table("test_table", charge_method='bcc-ml', pdbqt_method="third-party")

# Will use the antechamber package to compute gas charges
sample_project_cs.generate_mols_in_table("test_table", charge_method='gas', pdbqt_method="third-party")

```