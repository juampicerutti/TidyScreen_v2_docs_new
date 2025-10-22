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