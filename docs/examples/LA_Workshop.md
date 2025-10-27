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

# Process the raw pdb file to generate the receptor mol2 and pdqbt files  ##JUAMPI##: pdbqt files with prepare_receptor4.py?
czp_rx = "/PATH/TO/PDB/FILE/2OZ2.pdb"
catl_rx= "/PATH/TO/PDB/FILE/8hfv_aligned.pdb"

sample_project_moldock.process_raw_pdb(czp_rx, x_coord=6.7, y_coord=-3.7, z_coord=-22.0, x_points=40, y_points=40, z_points=40)

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
##JUAMPI##: Which files are generated? Is the .gpf included? It'd be nice to give a short explanation of the output!

#### 1.1.2. AutoDock4 grids generation

Define a grid box centered on the crystallographic ligand (A-site) and generate AutoGrid maps for AutoDock4.

With the grid parameter file prepared, the next step is to run **AutoGrid**.  
This program uses the receptor and the defined grid box to generate the energy maps required by AutoDock4.  
These maps describe the interaction potentials of the ligand atom types within the receptor binding site and are essential for the subsequent docking step.

```bash
autogrid4 -p 2OZ2.gpf
```

After completion, AutoGrid will produce:

- One .map file per ligand atom type (e.g.: `2OZ2.C.map`, `2OZ2.HD.map`, `2OZ2.N.map`, `2OZ2.NA.map`, `2OZ2.OA.map`, etc.)

- One .map file for electrostatic and desolvation energies (i.e.: `2OZ2.d.map`, `2OZ2.e.map`).

- A .fld file (`2OZ2.maps.fld`) that collects all maps.


Up to this point, we should possess all the necessary information about the receptors to conduct the docking studies.


### 1.2. Ligand preparation: K777 

In order to replicate the crystal structure, we will run docking studies with K777 and the receptors we just parametrized.

To do that, in the `$PATH/chemspace/raw_data/` create a file called `K777.csv` containing the SMILES of the ligand.

```jsx title="K777.csv"
O=S(C1=CC=CC=C1)(/C=C/[C@H](CCC2=CC=CC=C2)NC([C@H](CC3=CC=CC=C3)NC(N4CCN(C)CC4)=O)=O)=O,K777,5
```

Then, activate the chemspace module within `Workshop.py`.

```python title="Workshop.py"
# Activate the ChemSpace section of the project
Workshop_chemspace = chemspace.ChemSpace(Workshop)

# Import the .csv file
Workshop_chemspace.input_csv("$PATH/chemspace/raw_data/K777.csv")
```
After preparing the receptor, the next step is to generate the **ligand files** required for docking.  
Specifically, we need to obtain the `.pdbqt` formats for AutoDock.  
Since these molecules may also be subjected to **molecular dynamics simulations** later on, we will additionally generate `.mol2` and `.frcmod` files.  

Atomic charges are assigned using `prepare_ligand4.py` by default, but different charges can be set.

```python title="Workshop.py"
# Parametrization
Workshop_chemspace.generate_mols_in_table("K777", timeout=1000000)
```

At this stage, all the necessary input files have been prepared, and we are ready to proceed with the docking study.


### 1.3. First docking validation

In order to validate our docking protocol, the first step is to compare the predicted binding pose with the crystallographic reference to evaluate whether the chosen docking parameters (grid size, spacing, scoring function, and search algorithm) are able to reproduce the native interaction mode. A successful docking provides confidence that the protocol is suitable to be extended to other congeners and to virtual screening campaigns.

To proceed, we first activate the `moldock` module from *TidyScreen* and register the prepared receptor by creating a database entry that points to the target files (`.map` grids, `.pdbqt`, `.fld`, and reference `.pdb`).

```python title="Workshop.py"
# Initialize MolDock
Workshop_docking = moldock.MolDock(Workshop) 

# Register the receptor in the database, indicating the path to the folder containing the receptor/grid files.
Workshop_docking.input_receptor("$PATH/TO/PROJECT/docking/raw_data/2OZ2")
```

:::note[📝 Receptor description]

When prompted, you will need to provide a short **description of the receptor**.  
Try to include as many relevant details as possible, such as:

- PDB ID
- Preprocessing (waters/ions removed, protonation state)  
- Grid setup (center, npts, spacing)  
- Any special considerations (e.g. retained Mg²⁺)  

This information will be very useful in the future, especially when your database grows and you manage multiple receptor models.
:::


Before running the docking experiment, it is necessary to define the docking parameters.  
The function `create_docking_params_set()` generates a database (`docking/params/docking_params.db`) containing the **default AutoDock conditions**.  


```python title="Workshop.py"
Workshop_docking.create_docking_params_set()
```

Once the default set has been created, you may duplicate and customize the conditions to better suit your study.
For example, you can increase the number of docking runs from 20 to 100 (*--nrun=100*) or enable AutoDock to identify and report ligand–receptor interactions for each pose (*--contact_analysis=1*).

Now we can launch the docking experiment with the `dock_table()` function. 

:::note[📝 AutoDock]
TidyScreen prepares all required input files and executes the docking through AutoDock-GPU, the GPU-accelerated implementation of AutoDock4. This allows a substantial speed-up compared to the CPU version, while maintaining compatibility with the same scoring function and parameter set.
:::

You must define 3 variables:
- *table_name*: name of the table containing the ligand(s) you want to use for the docking study, as stated in the `chemspace/processed_data/chemspace.db` database.
- *id_receptor_model*: ID of the target model you want to use for the docking study, as stated in the `docking/receptors/receptor.db` database.
- *id_docking_params*: ID of the docking conditions you want to apply for the study, as stated in the `docking/params/docking_params.db` database.

Once again you will be asked to provide a brief **description** of the docking assay. Try to include as many relevant details as possible.

:::note[📝 IDs]
Make sure that the IDs (*id_receptor_model* and *id_docking_params*) correspond to the entries you created earlier.
If in doubt, you can query the databases to list available IDs before running the docking.
:::

```python title="Workshop.py"
Workshop_docking.dock_table(table_name="K777", id_receptor_model=1, id_docking_params=2)
```

After running `dock_table()`, *TidyScreen* will:

- Create a summary entry of the assay (variables + description) in `docking/docking_registers/docking_registers.db`.
- Create a new assay folder under `docking/docking_assays/assay_<ASSAY_ID>` containing:
  - *receptor/* and *ligands/* subfolders with all input files,
  - a ready-to-run execution script: `docking_execution.sh`.

You can launch the assay in your terminal with:

```jsx
./docking_execution.sh
```
### 1.4. Inspecting AutoDock outputs

After execution, AutoDock produces a *docking log* (`.dlg`) for each ligand under: `docking/docking_assays/assay_<ASSAY_ID>/dlgs/`.

Each `.dlg` file **summarizes**:
- the predicted **poses** (coordinates),
- the **interaction report** per pose (if `--contact_analysis=1`),
- the **docking score** (binding energy) of each conformation,
- and a **cluster size summary** at the end of the file.

:::note[📝 Cluster size]
The cluster size indicates how many independent runs converged to (approximately) the **same pose**.  
For example, if `--nrun=100` and a cluster reports `size = 80`, it means **80/100** runs found that pose as the best (lowest energy) within that cluster. Larger clusters often suggest **higher pose stability/reproducibility** under the search settings.
:::

### 1.5. Docking analyses

#### 1.5.1. Autodock output

In this first example, analysing the `.dlg` file and the docked poses is relatively straightforward, since the study involved only a single ligand. However, when working with larger ligand sets, manual inspection quickly becomes impractical.  For this reason, *TidyScreen* provides a `DockingAnalysis` module and a `process_docking_assay()` function to automate the post-processing of docking results.

:::warning[Important]
When activating the `DockingAnalysis` module, you can specify the path to **AMBER** (via `amberhome`) and **VMD**.  
This way, the program will automatically use this path during the analyses, instead of prompting you for it every time you run a single analysis.
:::

```python title="Workshop.py"
#Activate the module
Workshop_docking_analysis = docking_analysis.DockingAnalysis(Workshop, amberhome="/PATH/TO/AMBER/amber<VERSION>")
vmd_path="PATH/TO/VMD"
```

The `process_docking_assay()` function parses the AutoDock outputs of a given assay, leveraging the [Ringtail](https://ringtail.readthedocs.io/en/latest/) library to efficiently handle the results. It ranks the poses by docking score (most negative = best) and stores the selected results for downstream analysis.

You should define the following variables:

- *assay_id* (int, required): ID of the docking assay to analyze (as registered by dock_table()).
- *max_poses* (int, default = 10): Maximum number of poses per ligand to store in the results database, ordered from best (lowest binding energy) to worse.
- *vmd_path* (str | None, default = None): Path to VMD (Visual Molecular Dynamics). If provided, helper files/commands for quick visualization are prepared.
- *extract_poses* (int [0,1], default = 0): If 1, writes PDB files for each selected pose. If 0, keeps only the database records and existing PDBQT/DLG outputs.


```python title="Workshop.py"
#Set and run the docking analysis
Workshop_docking_analysis.process_docking_assay(assay_id=1, max_poses=5, vmd_path=vmd_path, extract_poses=1)
```

As a result, a new database `assay_<ID>.db`  will be generated inside the corresponding assay docking folder.

This file contains several tables that organize the results of the analysis:

* *interaction_indices*: dictionary of interaction types (e.g., hydrogen bonds H and van der Waals V) that can be identified between ligands and receptor.

* *interactions*: records which interactions were detected for each docked pose of every ligand.

* *ligands*: stores ligand-specific data and a summary of their detected interactions.

* *receptor*: summary information of the target used in the assay.

* *Results*: pose-level data extracted from the .dlg file, including binding energies, cluster ranks, and other docking details.

If *extract_poses=1*, a new subfolder *docked_1_per_cluster/* will be created inside the assay folder, containing the PDB files of the top docked poses (one representative per cluster) for each ligand, making it easier to visualize and compare binding modes.

##JUAMPI## Analyze scores and cluster sizes? Calculate rmsd? Visualize poses vs. crystallographic structure? 


#### 1.5.2. Computing interaction fingerprints for a docked pose

The function `compute_fingerprints_for_docked_pose()` automates the post-processing of **one specific docked pose**, generating:
- per-residue **ProLIF** interaction fingerprints (binary interaction features),
- per-residue **MMGBSA** interaction fingerprints (AMBER MMPBSA decomposition).


```python title="Workshop.py"
Workshop_docking_analysis.compute_fingerprints_for_docked_pose(
    assay_id,
    results_pose_id,
    mmgbsa=1,
    prolif=1,
    clean_files=1,
    clean_folder=1,
    solvent="implicit",
    min_steps=5000,
    store_docked_poses=1,
    prolif_parameters_set=1,
    iteration=1,
    ligresname="UNL",
)
```

**Parameters**

- *assay_id* (int, required): ID of the docking assay containing the pose.

- *results_pose_id* (int, required): ID of the pose to analyze (from the assay database results, e.g.: assay_X.pdb).

- *mmgbsa* (0/1, default=1): compute (1) or not (0) MMGBSA per-residue decomposition for the pose.

- *prolif* (0/1, default=1): compute (1) or not (0) ProLIF interaction fingerprint for the pose.

- *clean_files* (0/1, default=1): remove (1) or not (0) temporary files generated by MMPBSA/LEaP after finishing.

- *clean_folder* (0/1, default=1): remove (1) or not (0) the whole fingerprints_analyses/working folder after storing outputs.

- *solvent* ({"implicit","explicit"}, default="implicit"): solvent model used during minimization/prep for fingerprints.

- *min_steps* (int, default=5000): minimization steps for preparing the complex prior to fingerprints.

- *store_docked_poses* (0/1, default=1): store (1) or not (0) the docked pose in the results DB.

- *prolif_parameters_set* (int, default=1): ID of the parameter set in docking/params/prolif_parameters.db.

- *iteration* (int, default=1): if 1, builds a residue-mapping dictionary (crystal vs tleap) and saves it; otherwise reuses it.

- *ligresname* (str, default="UNL"): residue name used for the ligand in the complex (needed by LEaP/MMPBSA).

🤔 **What does the function do?**

**i)** Sets up a working directory under *docking/docking_assays/assay_`<ID>`/fingerprints_analyses/`<pose_subfolder>`/*.

**ii)** If *mmgbsa == 1*, prepares AMBER input (LEaP), minimizes the complex (min_steps, solvent), runs MMPBSA decomposition and writes a per-residue CSV, and stores the results in the DB (*store_mmbgsa_fingerprints_results_in_db*).

**iii)** If *prolif == 1*, builds a reference ProLIF dataframe for the entire receptor, mapped to the crystallographic numbering. It reuses the minimized prmtop/crd files if *mmgbsa == 1*; otherwise, these files are generated specifically for the ProLIF computation. Then, it retrieves the parameter set from *docking/params/prolif_parameters.db* according to the prolif_parameters_set ID, computes the ProLIF fingerprint for the docked pose and maps it to the crystallographic residue numbering, and saves the results as *prolif_fingerprints_renum.csv* and registers them in the database.

**iv)** If *store_docked_poses == 1*, the pose itself is stored in the results database.

**v)** Depending on the values of clean_files and clean_folder, the function will also remove intermediate files and/or the entire workspace once the computations are complete.

**vi)** Finally, the total execution time is logged and printed to the console.

As a practical example, we will calculate the interaction fingerprint of *Gentamicin C1a docked pose 1*, in order to evaluate whether the interactions identified experimentally can also be reproduced *in silico*. For simplicity, we will use the *default conditions* of the function.

```python title="Workshop.py"
Workshop_docking_analysis.compute_fingerprints_for_docked_pose(assay_id=1, results_pose_id=1)
```

#### 1.5.2.1. Interaction fingerprints with ProLIF

[ProLIF](https://prolif.readthedocs.io/en/latest/) is a Python library that encodes **target–ligand interactions** into machine-readable fingerprints.  
For each residue, it records whether specific interaction types were detected with the ligand in a given pose.

**What information does it provide?**  
Each column of the output corresponds to a residue–interaction pair, for example:  
- `GLN19_A_HBDonor` → hydrogen bond donor detected between residue Glutamine 19 of CZP and the ligand.  
- `TRP184_A_PiStacking` → π–π stacking between residue Tryptophan 184 of CZP and the ligand.  
- `MET145_A_VdWContact` → van der Waals contact between Methionine 145 of CZP and the ligand.  

Values are binary (0 or 1), indicating the **presence (1)** or **absence (0)** of the interaction in a given pose.  
This way, each pose can be represented as a vector of interactions, allowing systematic comparisons between poses or ligands.

**How to analyze it?**  
- At the **single-pose level**, you can identify which key interactions (e.g., H-bonds, hydrophobic contacts) stabilize the binding mode.  
- At the **comparative level**, you can check whether the best-scoring poses reproduce the same interactions observed in the crystallographic structure (e.g., conserved H-bonds in the A-site).  
- At the **dataset level**, these fingerprints can be clustered or used as features in **QSAR/ML models** to link interaction patterns with biological activity.

#### 1.5.2.2. Interaction fingerprints with MMGBSA

In addition to ProLIF, *TidyScreen* provides the option to minimize the docked pose of interest and evaluate the energetic contribution of each residue to the complex stabilization through **MMGBSA per-residue decomposition**. This approach not only reveals which contacts are present, but also quantifies how favorable or unfavorable each interaction is. Furthermore, it allows the calculation of the overall **Effective Binding Energy (EBE)** of the complex. While docking scores are useful for ranking poses, MMGBSA provides a more physically grounded estimation of binding affinity, often better correlated with experimental data. 

By setting `mmgbsa = 1` in the `compute_fingerprints_for_docked_pose()` function, a file named `mmgbsa_DECOMP_renum.csv` is generated, and a new `mmgbsa_fingerprints` table is added to the `assay_<ID>.db` database.

##JUAMPI## Is there any function to plot and/or analyze the results easily? Or should we provide scripts for that? (e.g.: heatmaps)

### 1.6. Molecular Dynamics?

## 2. Training set

