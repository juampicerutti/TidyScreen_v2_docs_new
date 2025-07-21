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
# Input the eMolecules database into ChemSpace from the corresponding 'emolecules_base.csv' file"
synthesis_example_chemspace.input_csv("$PATH/TO/emolecules_base.csv") 

# Note_1: we are not providing the full archive corresponding to the emolecules database. The user can obtain it upon request to the authors.
# Note_2: the input of the whole eMolecules database may take a while to load, since it contains millions of molecules.
```
Upon reading the reactants, a table named `emolecules_base` in the `chemspace.db` database located within the `$PROJECT_PATH/chemspace/processed_data/` folder. A quick inspection of the generated table is shown un Figure 3:

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


Lets filter out aminoacids. From a synthetic point of view, apart from the requirement that the building block contains an aminoacid residue, there are some additional requisites that are imposed by underlying synthetic mechanism, such as:

- only one aminoacid scaffold is present. Otherwise the reaction will lead to several by products;
- no primary amines in R. Otherwise the azidation will lead to the mixture of two azides;
- no azides originally contained in the aminoacid. Same as above;
- no tefminal alkynes.
- no thiols;
- no esters;
- no exotic atoms such as: boron, selenium, 13C, 2H, 3H, 15N.

As can be seen trough the corresponding listing, all these filters are already available by default when installing TidyScreen. In case the used would like to add custom filter (i.e. a fluorine atom), it can be done by using:

```python
synthesis_example_chemspace.add_smarts_filter("[F]","Fluorine atom") # Indicate the SMARTS, followed by a description of the filter

```


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
  <figcaption>**Figure 4:** Screenshot of the table containing the SMARTS reactions added to the project.</figcaption>
  </p>
</figure>
---



### Creation of a reaction workflow


### Execution of the reaction workflow