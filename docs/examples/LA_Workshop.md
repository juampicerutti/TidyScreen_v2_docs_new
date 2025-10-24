---
title: 1st Latin American hands-on workshop on open source AI/ML guided drug discovery
---

**Authors:** ... *Juan Pablo Cerutti (jpcerutti@unc.edu.ar) – October 2025*


## 1. Replication of the crystallographic binding mode of K777

### 1.1. Receptor preparation

#### 1.1.1. Download, clean and parametrize the crystal

The first step consists of preparing the targets from the crystallographic structures so that they can be used by AutoDock.  
This involves downloading the structure, removing non-essential molecules (ligand, solvent, crystallization agents), adding hydrogens at physiological pH, and finally generating a parameterized `.pdbqt` file with the correct atom types and charges.  




```python

# Download PDB of hCatL in complex with K777
wget https://files.rcsb.org/download/8hfv.pdb 

# Download PDB of CZP in complex with K777
wget https://files.rcsb.org/download/2OZ2.pdb 

# Align
lovoalign -p1 2OZ2.pdb -p2 8hfv.pdb -o CZP_2oz2.pdb

# Debería "limpiarse", quedándose solo con la cadena A. Estaría bueno aislar y guardar K777 en ambos casos como "pose cristalográfica de referencia".

```

[Download 8hfv](/files/8hfv.pdb)

[Download 2OZ2](/files/2OZ2.pdb)

[Download hCatL (pdb file we've been using PDB: 3OF8)](/files/hCatL.pdb)

[Download CZP1 (pdb file we've been using PDB: 3IUT)](/files/CZP1.pdb)

Once downloaded, the raw crystallographic structures can be processed by TidyScreen using the `process_raw_pdb()` function from the `MolDock` class as follows:

```python
>>> from tidyscreen import tidyscreen as ts
>>> tidyscreen.moldock import moldock as md

#### Activate the project. We will assume that a project 'sample_project' has been created
sample_project = ts.ActivateProject("sample_project")

## Activate the MolDock dimension to process the receptor models
sample_project_moldock = md.MolDock(sample_project)

# Process the raw pdb file to generate the receptor mol2 and pdqbt files
czp_rx = "/PATH/TO/PDB/FILE/2OZ2.pdb"
catl_rx= "/PATH/TO/PDB/FILE/8hfv_aligned.pdb"

sample_project_moldock.process_raw_pdb(czp_rx, x_coord=6.7,  y_coord=-3.7, z_coord=-22.0, x_points=40, y_points=40, z_points=40)

# Upon execution this method, a couple of checkings and questions will be presented:
> # Outputs:
>>> MolDock dimension activated.
>>> Please provide a brief description to store with the receptor model: # enter the info that will be stored to a 'description.txt' file in the receptor folder
>>> Chains detected: {'C', 'A'} # The presence of multi chains is performed and informed
>>> Provide the chain identifier to mantain: ("all" to keep all chains) # Enter the chain to extract, in this example we will use chain "A" 
>>> Non-standard residues detected:
>>> The following non-standard residues were found in the pdb file: ['D1R', 'SO4']. Do you want to mantain ONE of them in the processed receptor pdb file? (y/n): 
# We don't want to include non standard residues as part of the receptor, so just indicate 'n'
>>> The following non-standard residues were found in the pdb file: ['D1R', 'SO4']. Do you want to mantain ONE of them as a REFERENCE pdb file? (y/n):
# We would like to save a reference .pdb file containing the D1R to evaluate docking assays, so indicate 'y'
>>> Type the 3-letter code of the residue you want to save as REFERENCE FILE:
# Type D1R to select the residue

## The system should inform the successful processing of the receptor.

```

After the above process has been completed, a folder within `$PATH/TO/PROJECT/docking/raw_data` should have been created containing the prepared receptor folder (2OZ2)