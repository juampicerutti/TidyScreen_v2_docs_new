---
title: Receptors modeling and validation
---

**Step 1**: Retrieve receptors crystallographic structures

Using bash after project creation

```python
# Retrieve a structure corresponding to CZP
$ mkdir -p /PATH/TO/PROJECT/docking/raw_data/2OZ2
$ wget -O /PATH/TO/PROJECT/docking/raw_data/2OZ2/2OZ2.pdb https://files.rcsb.org/download/2OZ2.pdb

```

**Step 2**: Processing of raw crystallographic structures

Processing of CZP crystal structure pdb code: 2OZ2

```python 
>>> from tidyscreen import tidyscreen as ts
>>> from tidyscreen.moldock import moldock as md

### Activate the project
>>> la_workshop = ts.ActivateProject("la_workshop_2025")

### Instantiate a MolDock object to process the receptor models
>>> la_workshop_moldock = md.MolDock(la_workshop)

# Process CZP crystal structure
structure_file = "/PATH/TO/PROJECT/docking/raw_data/2OZ2/2OZ2.pdb"
>>> la_workshop_moldock.process_raw_pdb(structure_file, x_coord=6.7,  y_coord=-3.7, z_coord=-22.0, x_points=40, y_points=40, z_points=40) # Note that coordinates and size of the grid box are provided


# Outputs
Please provide a brief description to store with the receptor model: # Provide the requested information

# Outputs
Chains found in the pdb file: {'C', 'A'} # Two chains were detected in the crystal
Chains detected: {'C', 'A'} # ID of the detected chains
Provide the chain identifier to mantain: ('all' to keep all chains) # Indicate the chain to retain. In this tutorial we will use chain 'A'

# Outputs
Non-standard residues detected:
The following non-standard residues were found in the pdb file: ['D1R', 'SO4']. 
Do you want to mantain ONE of them as a REFERENCE pdb file? (y/n): # answer 'y' if a reference ligand is to be kept in a separate file. 

# Outputs

Type the 3-letter code of the residue you want to save as REFERENCE FILE:  # provide 'D1R' to retain K777 as reference ligands

# Outputs
The following non-standard residues were found in the pdb file: ['D1R', 'SO4']. 
Do you want to mantain ONE of them in the processed receptor pdb file? (y/n): # Answer 'n'(as many times as non-standard residues exists)

# Outputs
# Informs the creation of:
# /PATH/TO/PROJECT/docking/raw_data/2OZ2/receptor.mol2
# /PATH/TO/PROJECT/docking/raw_data/2OZ2/receptor.pdbqt
```

**Step 3**: Compute Autodock4 grid files

```bash
$ cd /PATH/TO/PROJECT/docking/raw_data/2OZ2/
$ autogrid4 -p receptor.gpf -l receptor.glg

# The corresponding grid maps will be created in the receptor folder
```

**Step 4**: Import the prepared receptor into TidyScreen

```python 
>>> la_workshop_moldock.input_receptor("/PATH/TO/PROJECT/docking/raw_data/2OZ2/")

# Outputs
Provide a brief description of the receptor model: # Input the required info

# Outputs
Succesfully stored receptor model
Successfully stored the receptor model in the database.
```

