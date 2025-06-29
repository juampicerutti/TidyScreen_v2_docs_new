---
title: Molecules input
sidebar_position: 1
---

# Molecules input

Every single virtual screening campaign initiates with the identification of a chemical space that is intended to be explored. This chemical space may vary significantly in size, ranging from a single molecule to a large set of hundred millions compounds.

In this respect, the most simple procedure to read a set of compound into TidyScreen is by using a `.csv` file structured in the following way:

***field_#1:*** [*SMILES*](https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html) notation indicating the structure of the molecule/s to be read.

***field_#2:*** name of the molecule/s

***field_#3:*** Flag used for general purposes (such as filtering, etc)

In order to read a `.csv` file, the corresponding function should be invoked using the generated `chemspace_object` as follows:

```python
# Step #1: It is assumed that TidyScreen packages have been activated
# Step #2: A project named `tutorial` has been generated
# Step #3: A chemspace object named `tutorial_chemspace` has been created
$ tutorial_chemspace.input_csv($PATH_TO_FILE)
```

When the above function is executed, the first line of the `.csv` file is written to the console, asking the user to indicate if the files has a header or not (y/n).

Once the input files containing valid SMILES notations, a table named after the read `.csv` file will be created in the `chemspace.db`.

---

Molecules tables contained within the `chemspace.db` can be listed as follows:

```python
# Step #1: It is assumed that TidyScreen packages have been activated
# Step #2: A project named `tutorial` has been generated
# Step #3: A chemspace object named `tutorial_chemspace` has been created
$ tutorial_chemspace.list_ligand_tables()
```

Also a table can be deleted from `chemspace.db` using:

```python
# Step #1: It is assumed that TidyScreen packages have been activated
# Step #2: A project named `tutorial` has been generated
# Step #3: A chemspace object named `tutorial_chemspace` has been created
$ tutorial_chemspace.delete_table(TABLE_NAME)
```

In case a 2D depiction of compounds stored in a given table is required, the following function can be invoked:

```python
# Step #1: It is assumed that TidyScreen packages have been activated
# Step #2: A project named `tutorial` has been generated
# Step #3: A chemspace object named `tutorial_chemspace` has been created
$ tutorial_chemspace.depict_ligand_table(TABLE_NAME)
```

The depiction of a table of molecules will generate a set of `.png` files that will be stored within the `$PROJECT_PATH/chemspace/misc/$TABLE_NAME` folder. 

:::note
Frequently, during virtual screening campaigns a high number of molecules are stored within a table. Consequently, depicting all of them may lead to a large number or `.png` files. In those cases, the following are usefull

- ***max_mols_ppage***: indicates the number of depictions included per `.png` file. Default: `25`
- ***limit***: set the *maximum* number of molecules to depict. Default: `0` ; depicts all molecules stored in the table.
- ***random***: If set to `True` will depict the required molecules indicated with *limit* picking a randomized selection. Default: `False`  
:::

