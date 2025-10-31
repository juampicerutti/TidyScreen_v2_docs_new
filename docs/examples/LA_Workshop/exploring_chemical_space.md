---
title: Chemical space design and vHTS workflow
---

# Chemical Space Design and vHTS Workflow

In this final section, we will extend the previously validated docking methodology to generate and evaluate **new triazole-based Targeted Covalent Inhibitors (TCIs)** of Cruzipain (CZP).

Following the successful validation of the docking workflow using a training set of triazole-containing derivatives featuring a methyleneamine linker as bioisosteric replacements of peptide bonds and the well-known efficacy of the vinylsulfone-containing inhibitor **K777**, our goal is to design a new series of molecules that **combine these three structural elements**.

Given the absence of reported analogues combining these motifs, we will perform a massive virtual library generation and screening campaign aimed at identifying candidates with **high affinity for CZP** and **minimal affinity for hCatL**, while ensuring **synthetic feasibility** for subsequent *wet-lab* validation.

---

## Overview of the Protocol

The overall workflow can be divided into five major steps:

1. **Extraction of building blocks**  
   Retrieve potential reactants from the *eMolecules* chemical database.

2. **Chemical space reduction**  
   Filter the combinatorial chemical space based on:
   - *Drug-likeness* criteria 
   - *Synthetic feasibility*
   - *Affordability*, using models provided through **Ersilia Hub**.

3. **In silico synthesis of derivatives**  
   Combine the selected reactants through feasible *in silico* synthetic routes to generate the final compounds, replicating realistic reactions that can be performed in the wet lab.

4. **Toxicity filtering**  
   Further reduce the library by eliminating potentially cytotoxic molecules, applying an **Ersilia Hub cytotoxicity prediction model**.

5. **Docking and analysis**  
   Perform docking-based *vHTS* on the safest and most promising analogues (approximately **5000 compounds**) and analyze their predicted affinity and selectivity profiles toward CZP and hCatL.

---

🎯 This will allow prioritization of candidates for *in vitro* and *in vivo* validation, completing the *dry*-to-*wet lab* pipeline proposed in this tutorial.

---

## 1. **Extraction of building blocks**  


#### 🧭 Synthetic plan

We will combine three independently tractable substituents (**R<sup>1</sup>**, **R<sup>2</sup>**, **R<sup>3</sup>**) on a common scaffold through a concise, wet-lab-ready sequence:

*i.)* **R<sup>1</sup>** will be introduced from commercial **amino acids**, which will be derivatized to the corresponding **azides** by means of *diazotransfer* reactions.

*ii.)* **R<sup>2</sup>** and **R<sup>3</sup>** will be derived from **aldehyde** and **methyleneamine** derivatives, respectively, which will be combined via *A3 coupling* to obtain **propargylamine** intermediates.

*iii.)* The azides and alkynes obtained from (i) and (ii) will be joined through *CuAAC* to yield **1,4-disubstituted 1,2,3-triazole** derivatives.

*iv.)* Finally, to introduce the **phenyl-vinylsulfone warhead (WH)**, a *Horner-Wadsworth-Emmons (HWE) olefination* will be replicated *in silico*.

---

Thus, we will generate lists of potential **reactants** - namely **amino acids**, **aldehydes**, and **methyleneamines** - that can be purchased and used as starting points for the *in silico* synthesis.  
To this purpose, we will filter compounds from the **eMolecules** database. This will ensure that all designed structures are **synthetically feasible and accessible** for *wet-lab* validation.


:::note[📝 About eMolecules]

**[eMolecules](https://www.emolecules.com/)** is a large supplier-aggregated catalog of **purchasable reagents** and **make-on-demand** building blocks. It was founded in 2005 with a vision to reduce drug discovery timelines through improved efficiencies in the compound search and acquisition process.
:::

:::warning[Important]

Please check if the *chemspace module* is active within `workshop.py`.

```python title="workshop.py"
la_workshop_cs = chemspace.ChemSpace(la_workshop)
```
:::

### 1.1. 🗂️ Import eMolecules database

Import the chemical database file from **eMolecules** (usually stored in `chemspace/raw_data/`) into the local **TidyScreen** project database (`chemspace/processed_data/chemspace.db`).

This step converts the original `.csv` file into a structured SQLite database that can be efficiently queried by TidyScreen for compound selection and filtering.

```python title="workshop.py"
emolecules_database_file = "/PATH/TO/FILE/emolecules.csv" #Set your own path!
la_workshop_cs.input_csv(emolecules_database_file)
```

### 1.2. 🧩 Chemical scaffolds filtering

To subset the types of compounds we need, **TidyScreen** uses chemical filters defined by **SMARTS patterns**.  
By default, a list of these filters is already included and can be inspected using the following function:


```python title="workshop.py"
# List available SMARTS filters to obtain reactants
la_workshop_cs.list_available_smarts_filters()
```
You should get the following output: 


```Available SMARTS filters:```  
```Filter_id: 1, Filter_Name: Aminoacids, SMARTS: [NX3H2,NX4+H3][CX4H]([*])[CX3H0](=[OX1])[OX2H,OX1-]```  
```...```  
```Filter_id: 10, Filter_Name: Aldehydes, SMARTS: [CX3H1](=O)[#6]```  
```Filter_id: 11, Filter_Name: PrimAmines, SMARTS: [NX3;H2;!$(NC=O)]```  
```Filter_id: 12, Filter_Name: SecAmines, SMARTS: [NX3;H1;!$(NC=O)]```  
```...```  
```Filter_id: 53, Filter_Name: sulfonamide, SMARTS: [SX4](=O)(=O)N```  

One of the classes of reactants we need corresponds to α-amino acids (`Filter_id = 1`).
While it is possible to directly select all eMolecules entries containing at least one such group, this approach would be inefficient and chemically irrational. The output would likely include compounds with
multiple amino acid groups (unsuitable for diazotransfer), additional reactive functional groups interfering with subsequent steps, rare isotopes or heteroatoms uncommon in drug-like molecules, chemically unstable or undesirable moieties, etc.

To refine the chemical space, TidyScreen allows creating sequential filtering workflows combining multiple SMARTS rules.

The function `create_smarts_filters_workflow()` accepts a dictionary in the form:

`{ filter_id : maximum_allowed_occurrences }`

Each key-value pair defines how many times a particular substructure is permitted in a molecule (e.g., 1:1 means exactly one α-amino acid group; 2:0 means exclude all boron-containing compounds).

In this example, we will filter molecules that:

- Contain exactly one α-amino acid group (1:1)
- Contain exactly one primary amine (11:1)
- Contain exactly one carboxylic acid group (26:1)
- Exclude molecules containing any of the following: Boron (2:0), Silicon (3:0), Azides (4:0), Terminal Alkynes (5:0), Deuterium (6:0), Tritium (7:0), <sup>13</sup>C (8:0), <sup>15</sup>N (9:0), Selenium (17:0), Amides (21:0), Thiols (25:0), Esters (35:0), Sulfonamides (53:0).


```python title="workshop.py"
# Create a custom filtering workflow for α-amino acids
la_workshop_cs.create_smarts_filters_workflow({1:1,11:1,26:1,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,17:0,21:0,25:0,35:0,53:0})
```

When prompted, you will need to provide a short **description of the filter** (i.e.: *"Filtering of amino acids compatible with CuAAC-centered reaction scheme"*).

`Provide a description for the SMARTS filters workflow: `

Once confirmed, the workflow is stored in the internal database.  
You can verify that it was correctly created using the `list_available_smarts_filters_workflows()` function:

```python title="workshop.py"
# List all available SMARTS filtering workflows
la_workshop_cs.list_available_smarts_filters_workflows()
```

*Output*: `Workflow_id: 1, Filter_Specs: {"1": 1, "11": 1, "26": 1, "2": 0, "3": 0, "4": 0, "5": 0, "6": 0, "7": 0, "8": 0, "9": 0, "17": 0, "21": 0, "25": 0, "35": 0, "53": 0}, Description: Filtering of amino acids compatible with CuAAC-centered reaction scheme`


Now we can apply the corresponding workflow to filter the desired building blocks using the function `subset_table_by_smarts_workflow()`. 
You need to specify two parameters: *The table name* (as stored in the chemspace.db) and the *Workflow ID* to be applied.

```python title="workshop.py"
# Filter α-amino acids
la_workshop_cs.subset_table_by_smarts_workflow("emolecules",1)
```

When prompted, provide a short **description of the subset** (i.e.: *"Amino acids for diazotransfer reaction"*).

`Enter a description for the subset:`

A new table containing only the filtered compounds will be created in the `chemspace.db`. Each table follows the naming convention `emolecules_subset_X`, where "X" corresponds to the subset number (incremental per run), not to the Workflow ID. This ensures that your filtered set is stored as an independent subset, ready to be used for subsequent in silico synthesis steps.

`Table 'emolecules_subset_1' created in: '/PATH/TO/PROJECT/chemspace/processed_data/chemspace.db'`  
`Succesfully subseted table: 'emolecules' by SMARTS filters workflow with ID: '1'`


Next, we will subset the two remaining classes of starting materials - **aldehydes** and **methyleneamines** - following the same procedure used for α-amino acids.

The following workflow will retain only molecules containing **one aldehyde group** (`Filter_id = 10`) while excluding other reactive or undesired functional groups.

```python title="workshop.py"
# Create a custom filtering workflow for aldehydes
la_workshop_cs.create_smarts_filters_workflow({10:1,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,11:0,12:0,17:0,21:0,25:0,26:0,35:0,53:0})

# Filter aldehydes
la_workshop_cs.subset_table_by_smarts_workflow("emolecules",2)
```

`Table 'emolecules_subset_2' created in: '/PATH/TO/PROJECT/chemspace/processed_data/chemspace.db'`  
`Succesfully subseted table: 'emolecules' by SMARTS filters workflow with ID: '2'`

The `emolecules_subset_2` table will now contain only aldehyde derivatives suitable for A3 coupling with methyleneamines.


For **methyleneamines**, we do not have a pre-defined SMARTS filter.  
However, TidyScreen provides a specific function for adding custom SMARTS-based filters: `add_smarts_filter()`.  
You simply need to specify the SMARTS pattern and a descriptive filter name.


```python title="workshop.py"
# Add a custom SMARTS filter for primary methylene amines
la_workshop_cs.add_smarts_filter("[NX3;H2][CX4;H2]","Primary_Amines_custom")
```

*Output*: `Successfully added SMARTS filter: '[NX3;H2][CX4;H2]'`

To confirm that the new SMARTS filter was added successfully, you can list all available filters again:

```python title="workshop.py"
# List available SMARTS filters to obtain reactants
la_workshop_cs.list_available_smarts_filters()
```

```Available SMARTS filters:```  
```Filter_id: 1, Filter_Name: Aminoacids, SMARTS: [NX3H2,NX4+H3][CX4H]([*])[CX3H0](=[OX1])[OX2H,OX1-]```  
```...```  
```Filter_id: 10, Filter_Name: Aldehydes, SMARTS: [CX3H1](=O)[#6]```  
```Filter_id: 11, Filter_Name: PrimAmines, SMARTS: [NX3;H2;!$(NC=O)]```  
```Filter_id: 12, Filter_Name: SecAmines, SMARTS: [NX3;H1;!$(NC=O)]```  
```...```  
```Filter_id: 53, Filter_Name: sulfonamide, SMARTS: [SX4](=O)(=O)N```  
```Filter_id: 54, Filter_Name: Primary_Amines_custom, SMARTS: [NX3;H2][CX4;H2]```  

Now that the new filter is available, we can create and apply a corresponding workflow to subset the desired primary methyleneamine derivatives:

```python title="workshop.py"
# Create a custom filtering workflow for primary methylene amines
la_workshop_cs.create_smarts_filters_workflow({54:1,10:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,12:0,17:0,21:0,25:0,35:0,53:0})

# Filter primary methylene amines
la_workshop_cs.subset_table_by_smarts_workflow("emolecules",3)
```

`Table 'emolecules_subset_3' created in: '/PATH/TO/PROJECT/chemspace/processed_data/chemspace.db'`  
`Succesfully subseted table: 'emolecules' by SMARTS filters workflow with ID: '3'`

At this point, your chemspace database contains three curated subsets ready for in silico synthesis.

### 1.3. 📸 Inspection of a sample of filtered building blocks

Before moving forward, it is useful to visually inspect a random sample of compounds from each curated subset to ensure that the filtering process performed as expected.

TidyScreen provides the function `depict_ligand_table()`, which automatically generates molecular depictions for a given table in the **chemspace database** and saves them in the `processed_data/misc` directory.


```python title="workshop.py"
# Depict a random sample of amino acids
la_workshop_cs.depict_ligand_table("emolecules_subset_1", limit=25, random=True) # Outputs to '/PATH/TO/PROJECT/chemspace/processed_data/misc'
# Depict a random sample of aldehydes
la_workshop_cs.depict_ligand_table("emolecules_subset_2", limit=25, random=True) # Outputs to '/PATH/TO/PROJECT/chemspace/processed_data/misc'
# Depict a random sample of primary methyleneamines
la_workshop_cs.depict_ligand_table("emolecules_subset_3", limit=25, random=True) # Outputs to '/PATH/TO/PROJECT/chemspace/processed_data/misc'
```

## 2. Building blocks prioritization

Compute default properties for Aldehydes: ["MolWt","MolLogP","NumHDonors","NumHAcceptors","NumRotatableBonds","TPSA"]

```python title="workshop.py"
la_workshop_cs.compute_properties("emolecules_subset_1")

# Outputs
The table named 'emolecules_subset_1' already exists in the database. I will replace it, are you ok with that? (y/n):
# Upon prompted, select 'y' to store the table with the new computed properties
```

Compute default properties for Amines: ["MolWt","MolLogP","NumHDonors","NumHAcceptors","NumRotatableBonds","TPSA"]

```python title="workshop.py"
la_workshop_cs.compute_properties("emolecules_subset_2")

# Outputs
The table named 'emolecules_subset_2' already exists in the database. I will replace it, are you ok with that? (y/n):
# Upon prompted, select 'y' to store the table with the new computed properties
```

Compute default properties for Amino acids: ["MolWt","MolLogP","NumHDonors","NumHAcceptors","NumRotatableBonds","TPSA"]

```python title="workshop.py"
la_workshop_cs.compute_properties("emolecules_subset_3")

# Outputs
The table named 'emolecules_subset_3' already exists in the database. I will replace it, are you ok with that? (y/n):
# Upon prompted, select 'y' to store the table with the new computed properties
```

Subset by drug-like properties the Aldehydes table ('emolecules_subset_1')

```python title="workshop.py"
la_workshop_cs.subset_table_by_properties("emolecules_subset_1",["MolWt>=200","MolWt<=500","MolLogP>=1.5","MolLogP<=3","NumRotatableBonds<=2"])

# Outputs
Enter a description for the subset: # Input for example 'Priotirized Aldehydes'
Succesfully subseted table: 'emolecules_subset_1' by properties: '['MolWt>=200', 'MolWt<=500', 'MolLogP>=1.5', 'MolLogP<=3', 'NumRotatableBonds<=2']'
```

Subset by drug-like properties the Amines table ('emolecules_subset_2')

```python title="workshop.py"
la_workshop_cs.subset_table_by_properties("emolecules_subset_2",["MolWt>=200","MolWt<=500","MolLogP>=1.5","MolLogP<=3","NumRotatableBonds<=2"])

# Outputs
Enter a description for the subset: # Input for example 'Priotirized Amines'
Succesfully subseted table: 'emolecules_subset_2' by properties: '['MolWt>=200', 'MolWt<=500', 'MolLogP>=1.5', 'MolLogP<=3', 'NumRotatableBonds<=2']'
```

Subset by drug-like properties the Amino Acids table ('emolecules_subset_3')

```python title="workshop.py"
la_workshop_cs.subset_table_by_properties("emolecules_subset_3",["MolWt>=200","MolWt<=500","MolLogP>=1.5","MolLogP<=3","NumRotatableBonds<=2"])

# Outputs
Enter a description for the subset: # Input for example 'Priotirized Amino Acids'
Succesfully subseted table: 'emolecules_subset_3' by properties: '['MolWt>=200', 'MolWt<=500', 'MolLogP>=1.5', 'MolLogP<=3', 'NumRotatableBonds<=2']'
```

Compute the prices using EOS eos7a45 model for Aldehydes

```python title="workshop.py"
la_workshop_cs.apply_ersilia_model_on_table("emolecules_subset_1_subset_4","eos7a45")
```

Compute the prices using EOS eos7a45 model for Amines

```python title="workshop.py"
la_workshop_cs.apply_ersilia_model_on_table("emolecules_subset_2_subset_5","eos7a45")
```

Compute the prices using EOS eos7a45 model for Amino Acids
```python title="workshop.py"
la_workshop_cs.apply_ersilia_model_on_table("emolecules_subset_3_subset_6","eos7a45")
```

Get the price ranges for Aldehydes
```bash
sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "SELECT MIN(eos7a45_coprinet), MAX(eos7a45_coprinet) FROM emolecules_subset_1_subset_4 WHERE eos7a45_coprinet IS NOT NULL;"

# Outputs
1.24|7.24
```

Get the price ranges for Amines
```bash
sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "SELECT MIN(eos7a45_coprinet), MAX(eos7a45_coprinet) FROM emolecules_subset_2_subset_5 WHERE eos7a45_coprinet IS NOT NULL;"

# Outputs
2.19|7.09
```

Get the price ranges for Aminoacids
```bash
sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "SELECT MIN(eos7a45_coprinet), MAX(eos7a45_coprinet) FROM emolecules_subset_3_subset_6 WHERE eos7a45_coprinet IS NOT NULL;"

# Outputs
$ 2.66|5.6
```

Create tables of cheap drug-like Aldehydes
```bash
sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "CREATE TABLE aldehydes_druglike_cheap AS SELECT * FROM emolecules_subset_1_subset_4 WHERE eos7a45_coprinet BETWEEN 1 AND 2.3 AND eos7a45_coprinet IS NOT NULL;" # Returns 135 building blocks
```

Create tables of cheap drug-like Amines
```bash
sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "CREATE TABLE aldehydes_druglike_cheap AS SELECT * FROM emolecules_subset_2_subset_5 WHERE eos7a45_coprinet BETWEEN 1 AND 2.6 AND eos7a45_coprinet IS NOT NULL;" # Returns 138 building blocks
```

Create tables of cheap drug-like Amino Acids
```bash
sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "CREATE TABLE aldehydes_druglike_cheap AS SELECT * FROM emolecules_subset_3_subset_6 WHERE eos7a45_coprinet BETWEEN 1 AND 3.5 AND eos7a45_coprinet IS NOT NULL;" # Returns 83 building blocks
```

**Step 9:** Combinatorial virtual synthesis

List available SMARTS reactions (will show empty)

```python title="workshop.py"
la_workshop_cs.list_available_smarts_reactions()

# Outputs
SMARTS reactions table does not exist yet. Add reactions to the database first.
```

Add the SMARTS reaction for diazotransfer
```python title="workshop.py"
la_workshop_cs.add_smarts_reaction("[NX3;H2:1][CX4:2][CX3,H0:3]>>[N-]=[N+]=[NX2;H0:1][CX4:2][CX3,H0:3]", "Diazotransfer")
```

Add the SMARTS reaction for an acylation intermediate step
```python title="workshop.py"
la_workshop_cs.add_smarts_reaction("[CX3:1](=[O:2])[OX2H,OX1-:3]>>[C:1](=[O:2])[O:3][C]", "Acylation")
```

Add the SMARTS reaction for A3 A3 coupling ion chiral - Note a specific stereochemistry is defined for reaction products
```python title="workshop.py"
la_workshop_cs.add_smarts_reaction("[N:1].[CX3H1:2](=[O:3])>>[NX4+:1][C@H:2][C:4]#[C:5]", "A3 coupling ion chiral")
```

Add the SMARTS reaction for CuAAC
```python title="workshop.py"
la_workshop_cs.add_smarts_reaction("[NX1-:1]=[NX2+:2]=[NX2:3].[CX2H1:4]#[CX2H0:5]>>[NX2+0:1]1=[NX2+0:2][N:3]-[C:4]=[C:5]1", "CuAAC")
```

Add the SMARTS reaction for the DIBAL reduction
```python title="workshop.py"
la_workshop_cs.add_smarts_reaction("[CX3:1](=[O:2])[OX2,OX1-:3]>>[CX3H1:1](=[O:2])", "DIBAL reduction")
```

Add the SMARTS reaction for the Horner_Wadsworth_Emmons olefination reaction
```python title="workshop.py"
la_workshop_cs.add_smarts_reaction("[CX4:1][CX3H1:2](=[O:3])>>[CX4:1]\[CX3H1:2]=[CX3H1]\[S](=O)(=O)c1ccccc1", "Horner_Wadsworth_Emmons Olefination")
```

List available SMARTS reactions
```python title="workshop.py"
la_workshop_cs.list_available_smarts_reactions()
```

Create a reaction workflow to obtain the target triazoles library: Diazotransfer -> Acylation -> A3 coupling ion chiral -> CuAAC -> DIBAL reduction -> Horner_Wadsworth_Emmons Olefination

```python title="workshop.py"
la_workshop_cs.add_smarts_reaction_workflow([1,2,3,4,5,6])
```

List available reaction workflow
```python title="workshop.py"
la_workshop_cs.list_available_reactions_workflows()
```

Execute the reaction workflow
```python title="workshop.py"
la_workshop_cs.apply_reaction_workflow(1,[["aminoacids_druglike_cheap"],["->:-1"],["amines_druglike_cheap","aldehydes_druglike_cheap"],["->:-2","->:-1"],["->:-1"],["->:-1"]])
```

Depict a sample of the generated virtual library
```python title="workshop.py"
la_workshop_cs.depict_ligand_table("reaction_set_1", limit=25, random=True)
```


**Step 10:** Candidates prioritization

Predict Human cytotoxicity endpoints based on eos21q7
```python title="workshop.py"
la_workshop_cs.apply_ersilia_model_on_table("reaction_set_1","eos21q7")
```

Check the cytotocicity ranges for the compounds

```bash
sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "SELECT MIN(eos21q7_dili_probability), MAX(eos21q7_dili_probability) FROM reaction_set_1 WHERE eos21q7_dili_probability IS NOT NULL;

# Outputs
0.28|0.73
``` 

Select the 5K safest compounds (lowest cytotoxicity probability) for molecular docking

```bash
sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "CREATE TABLE candidates_for_docking AS SELECT * FROM reaction_set_1 WHERE eos21q7_dili_probability IS NOT NULL ORDER BY eos21q7_dili_probability ASC LIMIT 5000;"
```

**Step 11:** Molecular docking of candidates

Prepare molecules in the candidates_for_docking table for docking
`python
la_workshop_cs.generate_mols_in_table("candidates_for_docking")
`

Instantiate a MolDock object to activate docking methods
```
la_workshop_moldock = md.MolDock(la_workshop)
```

Create custom docking conditions

```python title="workshop.py"
la_workshop_moldock.create_docking_params_set()
```

Here we need to open sqlite_web and duplciate the default parameters to update n_runs=100

```bash
sqlite_web docking/params/docking_params.db

# Will open a local window to edit the database
```

Docking of selected candidates
```python title="workshop.py"
la_workshop_moldock.dock_table("candidates_for_docking",id_receptor_model=3,id_docking_params=2)

# Remenber: receptor model 3 corresponds to the refined CZP receptor
```

Upon creation of the docking assay, 5 executabl files containing 1000 docking runs are created, the can be launched from the terminal:
```bash
# Note: the execution of the script below requires the availability of a local GPU board.
for i in $(seq 1 11); do ./docking_execution_${i}.sh ; done
```

**Step 12:** Docking results analysis


**Step 13:** Selection of candidates for chemical synthesis


