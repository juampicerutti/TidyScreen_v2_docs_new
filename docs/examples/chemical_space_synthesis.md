---
title: Chemical space synthesis
---


### Creating a custom chemical space using in silico synthesis

In this worked example we are going to follow the steps required to reproduce the construction of a virtual chemical space as reported by [Cerutti et al](https://www.mdpi.com/1420-3049/29/17/4224).

The execution of a virtual drug screening campaing encompasess diverse medicinal chemisty objectives, with each project constituting a unique scientific scenario. In this example, the aim of the work was to design potential triazole-based targeted covalent inhibitors of the antichagasic target Cruzipain.


---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/CZP.png" alt="Description of image" width="1000"/>
  <figcaption>**Figure 1:** Three-dimensional structure of the CZP, the target enzime for which ligands/inhibitors are to be constructed.</figcaption>
  </p>
</figure>
---

In this case, the researchers plans to perform a combinatorial exloration of R1, R2 and R3 on a central scaffold of a *1,2,3*-triazole ring. The core synthetic strategy is based on the use of Cu(I)-catalyzed 1,3-dipolar cycloaddition (CuAAC) reaction as is shown in Figure 2.

---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/triazole_synthesis.png" alt="Description of image" width="1000"/>
  <figcaption>**Figure 2:** Synthetic workflow designed for the obtention of 1,4 disubstituted 1,2,3-triazoles.</figcaption>
  </p>
</figure>
---

As can be seen, three starting building blocks are required to construct the intended chemical space:

- Aminoacids
- Secondary amines
- Aldehides

Also a trimethylsilyl (TMS) is involved early in the reaction workflow in order protect aminoacids prior to the azide formation stage.

### Obtention of commercially available reactants

Obtaining the building blocks commercially available can be accomplished from diverse public databases. In this example, we searched the [eMolecules](https://www.emolecules.com/) database, which at the moment of writting this tutorial comprised more than 18 million available reactants. 

As a first step, lets create and activate a dedicated project to prepare the intended compounds

```python
# Import required modules
from tidyscreen import tidyscreen as ts
from tidyscreen.chemspace import chemspace as chemspace

# Create a dedicated project for the synthesis
ts.create_project("$HOME/Desktop/example", "synthesis_example")

# Activate the general project created
synthesis_example = ts.ActivateProject("synthesis_example")

# Activate the ChemSpace section of the project
synthesis_example_chemspace = chemspace.ChemSpace(synthesis_example)
```

Next, import the full [eMolecules](https://www.emolecules.com/) database into the ChemSpace section of the project

```python
# Input the eMolecules database into ChemSpace from the corresponding 'emolecules.csv' file"
synthesis_example_chemspace.input_csv("$PATH/TO/emolecules.csv") 

# Note_1: we are not providing the full archive corresponding to the emolecules database. The user can obtain it upon request to the authors.
# Note_2: the input of the whole eMolecules database may take a while to load, since it contains millions of molecules.
```
Upon reading the reactants, a table named `emolecules` in the `chemspace.db` database located within the `$PROJECT_PATH/chemspace/processed_data/` folder. A quick inspection of the generated table is shown un Figure 3:

---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/emolecules_table_screenshot.png" alt="Description of image" width="1000"/>
  <figcaption>**Figure 3:** Screenshot of the table corresponding to the read eMolecules database.</figcaption>
  </p>
</figure>
---


Once the reactants `emolecules_base` table has been created, SMARTS filtering procedures can be applied to select the required building blocks. As a first step, lest check the available SMARTS filters. In this case, filters provided upong TidyScreen installation are shown:

```python
# List available SMARTS filters to obtain reactants
synthesis_example_chemspace.list_available_smarts_filters()
```

```python
# TidyScreen output
Available SMARTS filters:
Filter_id: 1, Filter_Name: Aminoacids, SMARTS: [NX3H2,NX4+H3][CX4H]([*])[CX3H0](=[OX1])[OX2H,OX1-]
...
Filter_id: 10, Filter_Name: Aldehydes, SMARTS: [CX3H1](=O)[#6]
Filter_id: 12, Filter_Name: SecAmines, SMARTS: [NX3;H1;!$(NC=O)]
...
```

As can be seen in the output, we are interested in using three of the filter provided with TidyScreen installation: Filters `1` (aminoacids), `10` (aldehydes) and `12` (secondary amines).


From a synthetic point of view, apart from the requirement that the building block contains an aminoacid residue, there are some additional requisites that are imposed by underlying synthetic mechanism, such as:

- only one aminoacid scaffold is present. Otherwise the reaction will lead to several by products;
- no primary amines in R. Otherwise the azidation will lead to the mixture of two azides;
- no azides originally contained in the aminoacid. Same as above;
- no tefminal alkynes.
- no thiols;
- no esters;
- no exotic atoms such as: boron, selenium, 13C, 2H, 3H, 15N.

As can be seen trough the corresponding listing, all these filters are already available by default when installing TidyScreen. In case the used would like to add custom filter (i.e. a fluorine atom), it can be done by using:

```python
>>> synthesis_example_chemspace.add_smarts_filter("[F]","Fluorine atom") # Indicate the SMARTS, followed by a description of the filter

```
### Construction and execution of a reactants filtering workflow

It is possible to concatenate multiple filters in a single workflow to apply all the filtering criteria at once. For example:

```python
# Create a workflow to filter aldehydes (10:1) for A3 coupling reactions (no interfering groups as indicated: 12:0, 11:0, etc)
>>> synthesis_example_chemspace.create_smarts_filters_workflow({10:1,12:0,11:0,1:0,4:0,5:0,25:0,35:0,26:0,2:0,3:0,17:0,8:0,6:0,7:0,9:0}) 

# Create a workflow to filter primary amines (11:1) for A3 coupling reactions reactions
>>> synthesis_example_chemspace.create_smarts_filters_workflow({11:1,12:0,1:0,10:0,4:0,5:0,25:0,35:0,26:0,2:0,3:0,17:0,8:0,6:0,7:0,9:0}) 

# Create a workflow to filter aminoacids (1:1) for click reactions
>>> synthesis_example_chemspace.create_smarts_filters_workflow({1:1,11:1,26:1,4:0,5:0,25:0,35:0,2:0,3:0,17:0,8:0,6:0,7:0,9:0}) 
```

Once create, available reactants filtering workflows can be listed:

```python
>>> synthesis_example_chemspace.list_available_smarts_filters_workflows()

### Outputs to terminal
Available SMARTS filters workflows:

Workflow_id: 1, Filter_Specs: {"10": 1, "12": 0, "11": 0, "1": 0, "4": 0, "5": 0, "25": 0, "35": 0, "26": 0, "2": 0, "3": 0, "17": 0, "8": 0, "6": 0, "7": 0, "9": 0}, Description: Filter Aldehydes to perform A3 coupling reactions 

Workflow_id: 2, Filter_Specs: {"11": 1, "12": 0, "1": 0, "10": 0, "4": 0, "5": 0, "25": 0, "35": 0, "26": 0, "2": 0, "3": 0, "17": 0, "8": 0, "6": 0, "7": 0, "9": 0}, Description: Filter primary amines for A3 coupling reactions 

Workflow_id: 3, Filter_Specs: {"1": 1, "11": 1, "26": 1, "4": 0, "5": 0, "25": 0, "35": 0, "2": 0, "3": 0, "17": 0, "8": 0, "6": 0, "7": 0, "9": 0}, Description: Filter Aminoacids for Click reactions
```

A filtering workflow can be applied on a given table (i.e. `emolecules`) available in the `chemspace.db` using:

```python
# Apply filter_id: 1 to the 'emolecules' table
>>> synthesis_example_chemspace.subset_table_by_smarts_workflow("emolecules",1)
```

Upon executing the filter, a table names `tables_subsets` will created in `chemspace.db`:


---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/tables_subsets_screenshot.png" alt="Description of image" width="1000"/>
  <figcaption>**Figure 4:** Screenshot of the subsets table containing the filtering actions.</figcaption>
  </p>
</figure>
---

The info included in the `tables_subsets` is:
- `table_name`: the source table on which the filtering workflow was applied.
- `subset_name`: the destination table in which the filtered compounds were written.
- `filtering_type`: the kind of filtering that originated the destination table. In this case using SMARTS notation
- `prop_filter`: explicit indication ot the filtering workflow applied.
- `description`: this is requested as information upon executing the filtering workflow.


After filtering the required aldehides, primary amines and aminoacids, the `tables_subsets` looks like this:

---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/tables_subsets_ready_screenshot.png" alt="Description of image" width="1000"/>
  <figcaption>**Figure 5:** Screenshot of the subsets table after filtering aldehides, primary amines and aminoacids.</figcaption>
  </p>
</figure>
---

In this specific example, the filtered tables contains:
- Aldehides: 316750 compounds
- Primary amines: 800298
- Aminoacids: 5713

A random set of each subseted table can be depicted:

```python 
# Randomly depict 25 aldehydes
>>> synthesis_example_chemspace.depict_ligand_table("emolecules_subset_1", limit=25,random=True)
# Randomly depict 25 primary amines
>>> synthesis_example_chemspace.depict_ligand_table("emolecules_subset_2", limit=25,random=True)
# Randomly depict 25 aminoacids
>>> synthesis_example_chemspace.depict_ligand_table("emolecules_subset_3", limit=25,random=True)
```

The resulting depictions are:

- **Aldehydes**
---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/emolecules_subset_1_0.png" alt="Description of image" width="1000"/>
  <figcaption>**Figure 6:** Set of 25 randomly picked aldehydes.</figcaption>
  </p>
</figure>
---


- **Primary amines**

---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/emolecules_subset_2_0.png" alt="Description of image" width="1000"/>
  <figcaption>**Figure 6:** Set of 25 randomly picked primary amines.</figcaption>
  </p>
</figure>
---


- **Aminoacids**

---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/emolecules_subset_3_0.png" alt="Description of image" width="1000"/>
  <figcaption>**Figure 6:** Set of 25 randomly picked aminoacids.</figcaption>
  </p>
</figure>
---

It is clear that exploring the while set of combinatorial possibilities towards 1,4 disubstituted 1,2,3-triazoles is not possible in terms of computational costs. Thus, some filtering (such as by drug-like properties) on the resulting subsets can be performed. In this way, the corresponding properties can be computed for each table subset: 

```python
## By default the properties are computed: ["MolWt","MolLogP","NumHDonors","NumHAcceptors","NumRotatableBonds","TPSA"]
### You can specify whatever property computed by RDKit you wan using the 'properties_list' keyword
>>> synthesis_example_chemspace.compute_properties("emolecules_subset_1")
```

Now the `emolecules_subset_1` contains additional columns included de computed properties:

---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/table_with_computed_properties.png" alt="Description of image" width="1000"/>
  <figcaption>**Figure X:** --- .</figcaption>
  </p>
</figure>
---

Now imagine you want to subset a given table containing computed properties using the following criteria:
- MolWt \>= 200 and MolWt \<= 500
- MolLogP \>=1.5 and MolLogP \<= 3 
- NumRotatableBonds \<= 2

you should construct the list: `["MolWt>=200","MolWt<=500","MolLogP>=1.5","MolLogP<=3","NumRotatableBonds<=2"]` and pass it to the corresponding subsetting method:

```python 
>>> synthesis_example_chemspace.subset_table_by_properties("emolecules_subset_1",["MolWt>=200","MolWt<=500","MolLogP>=1.5","MolLogP<=3","NumRotatableBonds<=2"])
```

Upon execution, a record will be stored in the `tables_subsets` table:

---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/subsetting_by_props.png" alt="Description of image" width="1000"/>
  <figcaption>**Figure X:** --- .</figcaption>
  </p>
</figure>
---





### Definition of single SMARTS reactions

It is now the moment to work with single reaction to be applied to the workflow. It should be noted that SMARTS based reactions belongs to a given project. To list available reactions for the current project, used the corresponding method as follows:

```python
synthesis_example_chemspace.list_available_smarts_reactions() 
```

Note that the first time you run this listing, since no reactions are available for the project, the message `"SMARTS reactions table does not exist yet. Add reactions to the database first."` will be printed to the terminal. To add user defined reactions in the SMARTS formats use: 

```python
synthesis_example_chemspace.add_smarts_reaction("SMARTS_REACTION", "Description") 
```

In the context of the chemical synthesis shown in Figure 2, we are interested in performing the following reactions as part of our synthetic scheme:

1. **Diazotransfer on aminocids:** `[NX3;H2:1][CX4:2][CX3,H0:3]>>[N-]=[N+]=[NX2;H0:1][CX4:2][CX3,H0:3]`
2. **Acylation of carboxylic acids:** `[CX3:1](=[O:2])[OX2H,OX1-:3].[C:4][O]>>[C:1](=[O:2])[O:3][C:4]`
3. **A3 Coupling reaction:** `[N:1].[CX3H1:2](=[O:3])>>[N:1][C:2][C:4]#[C:5]`
4. **CuAAC reaction:** `[NX1-:1]=[NX2+:2]=[NX2:3].[CX2H1:4]#[CX2H0:5]>>[NX2+0:1]1=[NX2+0:2][N:3]-[C:4]=[C:5]1`
5. **DIBAL reduction:** `[CX3:1](=[O:2])[OX2H,OX1-:3]>>[CX3H1:1](=[O:2])`

Lets add the corresponding reactions:

```python
# Add diazotransfer
synthesis_example_chemspace.add_smarts_reaction("[NX3;H2:1][CX4:2][CX3,H0:3]>>[N-]=[N+]=[NX2;H0:1][CX4:2][CX3,H0:3]", "Diazotransfer") 
# Add acylation
synthesis_example_chemspace.add_smarts_reaction("[CX3:1](=[O:2])[OX2H,OX1-:3].[C:4][O]>>[C:1](=[O:2])[O:3][C:4]", "Acylation") 
# Add A3 coupling
synthesis_example_chemspace.add_smarts_reaction("[N:1].[CX3H1:2](=[O:3])>>[N:1][C:2][C:4]#[C:5]", "A3 coupling") 
# Add CuAAC
synthesis_example_chemspace.add_smarts_reaction("[NX1-:1]=[NX2+:2]=[NX2:3].[CX2H1:4]#[CX2H0:5]>>[NX2+0:1]1=[NX2+0:2][N:3]-[C:4]=[C:5]1", "CuAAC") 
#Add DIBAL reduction
synthesis_example_chemspace.add_smarts_reaction("[CX3:1](=[O:2])[OX2H,OX1-:3]>>[CX3H1:1](=[O:2])", "DIBAL reduction") 
```

Upon adding reactions, a table named `smarts_reactions` will be created in the `chemspace.db` database:

---
<figure>
  <p align="left">
  <img src="/TidyScreen_v2_docs_new/img/smarts_reactions_table.png" alt="Description of image" width="500"/>
  <figcaption>**Figure X:** Screenshot of the table containing the SMARTS reactions added to the project.</figcaption>
  </p>
</figure>
---



### Creation of a reaction workflow


### Execution of the reaction workflow