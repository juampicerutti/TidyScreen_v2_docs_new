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

Subset by drug-like properties





**Step 9:** Combinatorial virtual synthesis


**Step 10:** Candidates prioritization


**Step 11:** Molecular docking of candidates


**Step 12:** Docking results analysis


**Step 13:** Selection of candidates for chemical synthesis


