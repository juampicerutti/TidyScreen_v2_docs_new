---
title: Chemical space generation
sidebar_position: 2
---
#
## Construction of molecules tables using *in silico* synthesis

A frequent scenario is that a user may want to generate a table of molecules by enumerating a synthetic route. In this case, the flow of information required is shown in the figure below:

---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/synthesis_pipeline.png" alt="Description of image" width="1000"/>
  <figcaption>**Figure 2:** Diagram showing the steps involved in a virtual synthetic workflow.</figcaption>
  </p>
</figure>
---

#### Step 1: Reactants input

Reactants may be read by following standard `.csv` files reading as was explained before. In this case, since only reactants compatible with the reaction procedure may be required, the filtering of specific functional groups may be required. This can be accomplished by using [SMARTS](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html). 

Upon installation of TidyScreen, a set of custom SMARTS filters are included, which can be retrieved using:

```python
# Step #1: It is assumed that TidyScreen packages have been activated
# Step #2: A project named `tutorial` has been generated
# Step #3: A chemspace object named `tutorial_chemspace` has been created
$ tutorial_chemspace.list_available_smarts_filters()
```

Upon execution of this function, the following output will be provided:

```python
Available SMARTS filters:
Filter_id: 1, Filter_Name: Aminoacids, SMARTS: [NX3H2,NX4+H3][CX4H]([*])[CX3H0](=[OX1])[OX2H,OX1-]
Filter_id: 2, Filter_Name: Boron, SMARTS: B
Filter_id: 3, Filter_Name: Silicon, SMARTS: [Si]
...
```

As can be seen, each SMART filter is associated to an `id` which is the key to be used when filtering compounds. 

New SMARTS filters can be added to the system by using:

```python
# Step #1: It is assumed that TidyScreen packages have been activated
# Step #2: A project named `tutorial` has been generated
# Step #3: A chemspace object named `tutorial_chemspace` has been created
$ tutorial_chemspace.add_smarts_filter(SMARTS_PATTERN,DESCRIPTION)
```
:::info
It is important to note that once a given SMARTS filter is added, it will be available to **all projects** in the local computer, since they are stored in the global TidyScreen database.
:::

It is commmon use that a serie of multiple chemical filters are applied to table, consequently filtering of reactants tables within TidyScreen is implemented as **workflow**. In this way, the first set to apply SMARTS filtering on reactants is **to create the given workkflow** as follows:

```python
# Step #1: It is assumed that TidyScreen packages have been activated
# Step #2: A project named `tutorial` has been generated
# Step #3: A chemspace object named `tutorial_chemspace` has been created
$ tutorial_chemspace.create_smarts_filters_workflow(DICT_OF_FILTERS_DEFINITION,DESCRIPTION)
```

:::info
The `DICT_OF_FILTERS_DEFINITION` required to construct a given filtering workflows is composed by:
- filter_id
- instances

So for example, if you would like to match molecules containing **only once** instance of filter ***1*** and excluding all ocurrences of filter ***2***, the following `dict` will be required:

```python
$DICT_OF_FILTERS = {1:1,2:0}
```
:::

It should be considered that filtering wrokflows will only be available to the current project, since they are written to the local project database. To **list** available filtering workkflows, the following function can be used:

```python
# Step #1: It is assumed that TidyScreen packages have been activated
# Step #2: A project named `tutorial` has been generated
# Step #3: A chemspace object named `tutorial_chemspace` has been created
$ tutorial_chemspace.list_available_smarts_filters()
```

:::info
As can be seen, each filtering workflow is assigned with an ID, which is required for it usage.
:::

Once the filter workflows is created, it can be applied as follows:

```python
# Step #1: It is assumed that TidyScreen packages have been activated
# Step #2: A project named `tutorial` has been generated
# Step #3: A chemspace object named `tutorial_chemspace` has been created
$ tutorial_chemspace.subset_table_by_smarts_workflow(TABLE_NAME,WORKFLOW_ID)
```

The table containing molecules matching the filtering criteria will be written to a new table named: `TABLE_NAME_subset_#` and a filtering registry will be created into a table called `tables_subsets` which contains the following information:
- Source table name
- Destination table name
- The type of filtering applyied (in this example is: `by_smarts`)
- The applied property filter: (in this example the corresponding SMARTS filters)
- A description, which which provided by the used upon creating the filtering workflow.

:::tip
After working for some time in a same project, a high number of filtering procedures may be applied, so the `tables_subsets` is quite useful to remember the applied procedure and the destination table. In case the user would like to make a copy of the filtered table with a more *descriptive* name, the following function can be used:

```python
# Step #1: It is assumed that TidyScreen packages have been activated
# Step #2: A project named `tutorial` has been generated
# Step #3: A chemspace object named `tutorial_chemspace` has been created
$ tutorial_chemspace.copy_table_to_new_name(SOURCE_TABLE_NAME,TARGET_TABLE_NAME)
```
:::

#### Step 2: Create reactions

In order to be able to execute synthetic workflows, the corresponding individual reactions **must** be defined within the working **TidyScreen project**. As oposed to the case of SMARTS filters, default reactions are not created by default upon installing TidyScreen, since it is considered that they are very specific to a given screening project. 

As was the case for reactants, [**reactions are also defined by using SMARTS/SMIRKS**](https://www-daylight-com.translate.goog/dayhtml/doc/theory/theory.smirks.html?_x_tr_sl=en&_x_tr_tl=es&_x_tr_hl=es&_x_tr_pto=tc) chemical transformation language.

To include a new reaction in the chemical transformation toolbox of the active project, the following function is used:

```python
# Step #1: It is assumed that TidyScreen packages have been activated
# Step #2: A project named `tutorial` has been generated
# Step #3: A chemspace object named `tutorial_chemspace` has been created
$ tutorial_chemspace.add_smarts_reaction(SMARTS_PATTERN,DESCRIPTION)
```

:::info
It is important to note that the added reactions will only be available for use within the given active project (i.e. they are stored in the specific project database).
:::

Available reactions can be listed using:

```python
# Step #1: It is assumed that TidyScreen packages have been activated
# Step #2: A project named `tutorial` has been generated
# Step #3: A chemspace object named `tutorial_chemspace` has been created
$ tutorial_chemspace.list_available_smarts_reactions()
```

:::info
Listing the available chemical reactions will provide the used with an ID, which be be usefull to the user to construct one-step/multi-step synthetic routes as will be presented in the following section
:::

#### Step 3: Workflow definition

In TidyScreen, a multi-step synthetic plan can be constructed by coupling single reactions. Ech synthetic plan is identified by a unique ID, and will be indicated by the user to construct the required products.

To construct a synthesis workflow, the following function is used:

```python
# Step #1: It is assumed that TidyScreen packages have been activated
# Step #2: A project named `tutorial` has been generated
# Step #3: A chemspace object named `tutorial_chemspace` has been created
$ tutorial_chemspace.add_smarts_reaction_workflow(LIST_IF_REACTIONS)
```

:::info
The `LIST_OF_REACTIONS` corresponds to a Python list including the IDs of the single reactions that are to be chained. A single reaction can also be added, but it should be provided in the form of a `list`
:::

Available reaction workflows can be listed as follows:

```python
# Step #1: It is assumed that TidyScreen packages have been activated
# Step #2: A project named `tutorial` has been generated
# Step #3: A chemspace object named `tutorial_chemspace` has been created
$ tutorial_chemspace.list_available_reactions_workflows()
```

:::info
Upon listing the workflows, the user will get an ID that is required to execute it.
:::

#### Step 4: workflow execution

Once reactants have been defined, individual reactions made available and sythetic protocols designed as a workflow, it it time for setting up and execute the synthetic scheme. Using the workflow ID, it can be executed as follows:

```python
# Step #1: It is assumed that TidyScreen packages have been activated
# Step #2: A project named `tutorial` has been generated
# Step #3: A chemspace object named `tutorial_chemspace` has been created
$ tutorial_chemspace.apply_reaction_workflow(WORKFLOW_ID,REACTANTS_LIST)
```

There is nothing special about the `WORKFLOW_ID`, it is just the ID the user can find by listing the available workflows as explained above.

:::info
The definition of the `REACTANTS_LIST` need a carefull setup in order to comply with the chemical logic as required with the synthetic workflow.
This `REACTANTS_LIST` is a **list of list**, in which each item (i.e. list) corresponds to the reactants required for each step defined in the synthetic protocol. 

In order to better understand the way this list needs to be constructed, it is better to analyze the examples below.
:::

:::tip ***Example_1***:
`REACTANTS_LIST = [["reactants_1_table"]]` : This is a one-step unimolecular reaction, in which the reaction SMARTS defined in the correspondin workflow ID matched molecules stored in table `reactants_1`. The simplest case posible.
:::

:::tip ***Example_2***:
`REACTANTS_LIST = [["reactants_1_table"],["reactants_2_table"]]` : This is reaction scheme composed of two sequential and **independent** unimolecular steps. That is, once the reaction workflow is read, two reactions are expected to happen, the first one on `"reactants_1_table"` and the second one on `"reactants_2_table"`. **It should be considered that TidyScreen will only store the results of the last reaction accomplished as part of a workflow, since it is assumed that those are the products aimed for screening. In this example, only products obtained on `"reactants_2_table"` will be stored.  
:::

:::tip ***Example_3***:
`REACTANTS_LIST = [["reactants_1_table"],["->:-1"]]` : This is reaction scheme composed of two sequential and **dependent** unimolecular steps. The first reactant table has nothing new, however the second one contains: **`"->:-1"`**. This string is saying: *"use the products* (**`->`**) *of the immediate previous step (**`-1`**)"*.   
:::

:::tip ***Example_4***:
`REACTANTS_LIST = [["reactants_1_table","reactants_2_table"],["->:-1"]]` : This is reaction scheme composed of two sequential and **dependent** steps, **the first one being bimolecular and the second one unimolecular**.   
:::

:::tip ***Example_5***:
`REACTANTS_LIST = [["reactants_1_table","reactants_2_table"],["->:-1","reactants_3_table"],["->:-1"]]` : This is reaction scheme composed of two sequential and **dependent** steps, **the first one being bimolecular, the second again bimolecular and involving as the first set of reagent the products obtained in the previous step**, ending with a unimolecular reaction starting the the products of the second step.   
:::

:::tip ***Example_6***:
`REACTANTS_LIST = [["reactants_1_table"],["reactants_2_table"],["->:-2""->:-1"]]` : This example is more elaborated, involving 3 reaction stages: a) the first and second ones being independent unimolecular reactions, while the third is a dependent bimolecular reaction, using as first reactants the products obtained from step 1 (*note the -2 indexing*) combined the products from step 2 (*note the -1 indexing*).   
:::

:::warning
***Consider that the order in which the tables of reactants are defined in the each item of the `REACTANTS_LIST` should match that of the definition in the corresponding SMARTS reaction.
:::


