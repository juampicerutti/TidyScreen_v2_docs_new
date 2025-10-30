---
title: Exploring a custom chemical space
---

**Step 1:** Import emolecules database

```python
# Import the emolecules chemical database from a CSV file
>>> emolecules_database_file = "/PATH/TO/FILE/emolecules.csv"
>>> la_workshop_cs.input_csv(emolecules_database_file)
```

**Step 2:** Inspection of built-in chemical scaffolds filters

```python
# List available SMARTS filters to obtain reactants
>>> la_workshop_cs.list_available_smarts_filters()

> # Outputs:
>>> Available SMARTS filters:
>>> Filter_id: 1, Filter_Name: Aminoacids, SMARTS: [NX3H2,NX4+H3][CX4H]([*])[CX3H0](=[OX1])[OX2H,OX1-]
>>> ...
>>> Filter_id: 10, Filter_Name: Aldehydes, SMARTS: [CX3H1](=O)[#6]
>>> Filter_id: 11, Filter_Name: PrimAmines, SMARTS: [NX3;H2;!$(NC=O)]
>>> Filter_id: 12, Filter_Name: SecAmines, SMARTS: [NX3;H1;!$(NC=O)]
...
```

**Step 3:** Add a custom SMARTS filter to extract primary and secondary amines
```python
# Add a custom SMARTS filter
>>> la_workshop_cs.add_smarts_filter("[NX3;H1,H2;!$(NC=O)]","Primary_and_Secondary_Amines")

> # Outputs:
>>> Successfully added SMARTS filter: '[NX3;H1,H2;!$(NC=O)]'
```

**Step 4:** Check that the new SMARTS filter has been successfully added

```python
# List available SMARTS filters to obtain reactants
>>> la_workshop_cs.list_available_smarts_filters()

> # Outputs:
>>> Available SMARTS filters:
>>> Filter_id: 1, Filter_Name: Aminoacids, SMARTS: [NX3H2,NX4+H3][CX4H]([*])[CX3H0](=[OX1])[OX2H,OX1-]
>>> ...
>>> Filter_id: 10, Filter_Name: Aldehydes, SMARTS: [CX3H1](=O)[#6]
>>> Filter_id: 11, Filter_Name: PrimAmines, SMARTS: [NX3;H2;!$(NC=O)]
>>> Filter_id: 12, Filter_Name: SecAmines, SMARTS: [NX3;H1;!$(NC=O)]
...
>>> Filter_id: 54, Filter_Name: Primary_and_Secondary_Amines, SMARTS: [NX3;H1,H2;!$(NC=O)]
```

**Step 5:** Create custom filtering workflows


```python
# Create a custom filtering workflow for aldehydes
>>> la_workshop_cs.create_smarts_filters_workflow({10:1,11:0,26:0,12:0,4:0,5:0,25:0,35:0,21:0,53:0,2:0,17:0,6:0,7:0,8:0,9:0})

# Outputs
>>> SMARTS workflows table does not exist yet. Creating it...
>>> Provide a description for the SMARTS filters workflow: 

# The user is prompted for a description of the objective of this filtering workflow.
# i.e.: "Filtering of aldehydes compatible with A3 coupling reactions"
```

```python
# Create a custom filtering workflow for primary and secondary amines
>>> la_workshop_cs.create_smarts_filters_workflow({54:1,10:0,4:0,5:0,25:0,35:0,21:0,53:0,2:0,17:0,6:0,7:0,8:0,9:0})

# Outputs
>>> Provide a description for the SMARTS filters workflow: 

# The user is prompted for a description of the objective of this filtering workflow.
# i.e.: "Filtering of amines compatible with A3 coupling reactions"
```

```python
# Create a custom filtering workflow for primary and secondary amines
>>> la_workshop_cs.create_smarts_filters_workflow({1:1,11:1,26:1,4:0,5:0,25:0,35:0,21:0,53:0,2:0,17:0,6:0,7:0,8:0,9:0})

# Outputs
>>> Provide a description for the SMARTS filters workflow: 

# The user is prompted for a description of the objective of this filtering workflow.
# i.e.: "Filtering of amino acids compatible with CuAAC coupling reactions"
```

**Step 6:** Inspect available filtering wokflows

```python
# Create a custom filtering workflow for primary and secondary amines
>>> la_workshop_cs.list_available_smarts_filters_workflows()

# Outputs
>>> Workflow_id: 1, Filter_Specs: {"10": 1, "11": 0, "26": 0, "12": 0, "4": 0, "5": 0, "25": 0, "35": 0, "21": 0, "53": 0, "2": 0, "17": 0, "6": 0, "7": 0, "8": 0, "9": 0}, Description: Filtering of aldehydes compatible with A3 coupling reactions 

>>> Workflow_id: 2, Filter_Specs: {"54": 1, "10": 0, "4": 0, "5": 0, "25": 0, "35": 0, "21": 0, "53": 0, "2": 0, "17": 0, "6": 0, "7": 0, "8": 0, "9": 0}, Description: Filtering of amines compatible with A3 coupling reactions 

>>> Workflow_id: 3, Filter_Specs: {"1": 1, "11": 1, "26": 1, "4": 0, "5": 0, "25": 0, "35": 0, "21": 0, "53": 0, "2": 0, "17": 0, "6": 0, "7": 0, "8": 0, "9": 0}, Description: Filtering of amino acids compatible with CuAAC coupling reactions 
```


**Step 7:** Apply the corresponding workflows to obtained filtered building blocks

Filter aldehydes

```python
# Filter aldehydes
la_workshop_cs.subset_table_by_smarts_workflow("emolecules",1)

# Outputs
>>> Enter a description for the subset:
# The user is prompted for a description of the resulting subset.
# i.e.: "Aldehydes for A3 coupling reactions"

# Outputs
>>> Table 'emolecules_subset_1' created in: '/PATH/TO/PROJECT/chemspace/processed_data/chemspace.db'
>>> Succesfully subseted table: 'emolecules' by SMARTS filters workflow with ID: '1'
```


Filter primary and secondary amines

```python
# Filter amines
la_workshop_cs.subset_table_by_smarts_workflow("emolecules",2)

# Outputs
>>> Enter a description for the subset:
# The user is prompted for a description of the resulting subset.
# i.e.: "Amines for A3 coupling reactions"

# Outputs
>>> Table 'emolecules_subset_2' created in: '/PATH/TO/PROJECT/chemspace/processed_data/chemspace.db'
>>> Succesfully subseted table: 'emolecules' by SMARTS filters workflow with ID: '2'
```

Filter amino acids

```python
# Filter amines
la_workshop_cs.subset_table_by_smarts_workflow("emolecules",3)

# Outputs
>>> Enter a description for the subset:
# The user is prompted for a description of the resulting subset.
# i.e.: "Amino acids for CuAAC reactions"

# Outputs
>>> Table 'emolecules_subset_3' created in: '/PATH/TO/PROJECT/chemspace/processed_data/chemspace.db'
>>> Succesfully subseted table: 'emolecules' by SMARTS filters workflow with ID: '3'
```

**Step 8:** Building blocks prioritization

Compute default properties for Aldehydes: ["MolWt","MolLogP","NumHDonors","NumHAcceptors","NumRotatableBonds","TPSA"]

```python
>>> la_workshop_cs.compute_properties("emolecules_subset_1")

# Outputs
>>> The table named 'emolecules_subset_1' already exists in the database. I will replace it, are you ok with that? (y/n):

# Upon prompted, select 'y' to store the table with the new computed properties
```
Compute default properties for Amines: ["MolWt","MolLogP","NumHDonors","NumHAcceptors","NumRotatableBonds","TPSA"]

```python
>>> la_workshop_cs.compute_properties("emolecules_subset_2")

# Outputs
>>> The table named 'emolecules_subset_2' already exists in the database. I will replace it, are you ok with that? (y/n):

# Upon prompted, select 'y' to store the table with the new computed properties
```

Compute default properties for Amino acids: ["MolWt","MolLogP","NumHDonors","NumHAcceptors","NumRotatableBonds","TPSA"]

```python
>>> la_workshop_cs.compute_properties("emolecules_subset_3")

# Outputs
>>> The table named 'emolecules_subset_3' already exists in the database. I will replace it, are you ok with that? (y/n):

# Upon prompted, select 'y' to store the table with the new computed properties
```

Subset by drug-like properties the Aldehydes table ('emolecules_subset_1')

```python
la_workshop_cs.subset_table_by_properties("emolecules_subset_1",["MolWt>=200","MolWt<=500","MolLogP>=1.5","MolLogP<=3","NumRotatableBonds<=2"])

# Outputs
>>> Enter a description for the subset: # Input for example 'Priotirized Aldehydes'
>>> Succesfully subseted table: 'emolecules_subset_1' by properties: '['MolWt>=200', 'MolWt<=500', 'MolLogP>=1.5', 'MolLogP<=3', 'NumRotatableBonds<=2']'
```

Subset by drug-like properties the Amines table ('emolecules_subset_2')

```python
la_workshop_cs.subset_table_by_properties("emolecules_subset_2",["MolWt>=200","MolWt<=500","MolLogP>=1.5","MolLogP<=3","NumRotatableBonds<=2"])

# Outputs
>>> Enter a description for the subset: # Input for example 'Priotirized Amines'
>>> Succesfully subseted table: 'emolecules_subset_2' by properties: '['MolWt>=200', 'MolWt<=500', 'MolLogP>=1.5', 'MolLogP<=3', 'NumRotatableBonds<=2']'
```

Subset by drug-like properties the Amino Acids table ('emolecules_subset_3')

```python
la_workshop_cs.subset_table_by_properties("emolecules_subset_3",["MolWt>=200","MolWt<=500","MolLogP>=1.5","MolLogP<=3","NumRotatableBonds<=2"])

# Outputs
>>> Enter a description for the subset: # Input for example 'Priotirized Amino Acids'
>>> Succesfully subseted table: 'emolecules_subset_3' by properties: '['MolWt>=200', 'MolWt<=500', 'MolLogP>=1.5', 'MolLogP<=3', 'NumRotatableBonds<=2']'
```

Compute the prices using EOS eos7a45 model for Aldehydes

```python
la_workshop_cs.apply_ersilia_model_on_table("emolecules_subset_1_subset_4","eos7a45")
```

Compute the prices using EOS eos7a45 model for Amines
```python
la_workshop_cs.apply_ersilia_model_on_table("emolecules_subset_2_subset_5","eos7a45")
```

Compute the prices using EOS eos7a45 model for Amino Acids
```python
la_workshop_cs.apply_ersilia_model_on_table("emolecules_subset_3_subset_6","eos7a45")
```

Get the price ranges for Aldehydes
```bash
$ sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "SELECT MIN(eos7a45_coprinet), MAX(eos7a45_coprinet) FROM emolecules_subset_1_subset_4 WHERE eos7a45_coprinet IS NOT NULL;"

# Outputs
$ 1.24|7.24
```

Get the price ranges for Amines
```bash
>>> sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "SELECT MIN(eos7a45_coprinet), MAX(eos7a45_coprinet) FROM emolecules_subset_2_subset_5 WHERE eos7a45_coprinet IS NOT NULL;"

# Outputs
>>> 1.53|7.47
```

Get the price ranges for Aminoacids
```bash
$ sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "SELECT MIN(eos7a45_coprinet), MAX(eos7a45_coprinet) FROM emolecules_subset_3_subset_6 WHERE eos7a45_coprinet IS NOT NULL;"

# Outputs
$ 2.66|5.6
```

Create tables of cheap drug-like Aldehydes
```bash
$ sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "CREATE TABLE aldehydes_druglike_cheap AS SELECT * FROM emolecules_subset_1_subset_4 WHERE eos7a45_coprinet BETWEEN 1 AND 2.3 AND eos7a45_coprinet IS NOT NULL;"
```

Create tables of cheap drug-like Amines
```bash
$ sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "CREATE TABLE aldehydes_druglike_cheap AS SELECT * FROM emolecules_subset_2_subset_5 WHERE eos7a45_coprinet BETWEEN 1 AND 2.3 AND eos7a45_coprinet IS NOT NULL;"
```

Create tables of cheap drug-like Amino Acids
```bash
$ sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "CREATE TABLE aldehydes_druglike_cheap AS SELECT * FROM emolecules_subset_3_subset_6 WHERE eos7a45_coprinet BETWEEN 1 AND 3.5 AND eos7a45_coprinet IS NOT NULL;"
```

**Step 9:** Combinatorial virtual synthesis

List available SMARTS reactions (will show empty)

```python
>>> la_workshop_cs.list_available_smarts_reactions()

# Outputs
>>> SMARTS reactions table does not exist yet. Add reactions to the database first.
```

Add the SMARTS reaction for diazotransfer
```python
la_workshop_cs.add_smarts_reaction("[NX3;H2:1][CX4:2][CX3,H0:3]>>[N-]=[N+]=[NX2;H0:1][CX4:2][CX3,H0:3]", "Diazotransfer")
```

Add the SMARTS reaction for A3 coupling
```python
la_workshop_cs.add_smarts_reaction("[N:1].[CX3H1:2](=[O:3])>>[N:1][C:2][C:4]#[C:5]", "A3_Coupling")
```

Add the SMARTS reaction for CuAAC
```python
la_workshop_cs.add_smarts_reaction("[NX1-:1]=[NX2+:2]=[NX2:3].[CX2H1:4]#[CX2H0:5]>>[NX2+0:1]1=[NX2+0:2][N:3]-[C:4]=[C:5]1", "CuAAC")
```

List available SMARTS reactions
```python
la_workshop_cs.list_available_smarts_reactions()
```

Create a reaction workflow A3: coupling -> diazotransfer -> CuAAC
```python
la_workshop_cs.add_smarts_reaction_workflow([2,1,3])
```

List available reaction workflow
```python
la_workshop_cs.list_available_reactions_workflows()
```

Execute the reaction workflow
```python
la_workshop_cs.apply_reaction_workflow(1,[["amines_druglike_cheap","aldehydes_druglike_cheap"],["aminoacids_druglike_cheap"],["->:-1","->:-2"]])
```

Predict Human cytotoxicity endpoints based on eos21q7
```python
la_workshop_cs.apply_ersilia_model_on_table("reaction_set_1","eos21q7")
```

**Step 10:** Candidates prioritization

Check the cytotocicity ranges for the compounds

```bash
$ sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "SELECT MIN(eos21q7_dili_probability), MAX(eos21q7_dili_probability) FROM reaction_set_1 WHERE eos21q7_dili_probability IS NOT NULL;

# Outputs
$ 0.3|0.81
``` 

Select the 20K safest compounds (lowest cytotoxicity probability) for molecular docking
```bash
$ sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "CREATE TABLE candidates_for_docking AS SELECT * FROM reaction_set_1 WHERE eos21q7_dili_probability IS NOT NULL ORDER BY eos21q7_dili_probability ASC LIMIT 20000;"


**Step 11:** Molecular docking of candidates





**Step 12:** Docking results analysis


**Step 13:** Selection of candidates for chemical synthesis


