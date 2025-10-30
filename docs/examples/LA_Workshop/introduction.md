---
title: Introduction to TidyScreen projects
---

Creation of a new project

```python
>>> from tidyscreen import tidyscreen as ts

### Create a new TidyScreen project named "la_workshop_2025" in the working directory
>>> ts.create_project("/PATH/TO/PROJECT/DIR","la_workshop_2025")
```

Project available locally can be listed

```python
>>> ts.projects()

# Outputs
Input project short description: # Provide the requested information

# Outputs
Project 'la_workshop_2025' created at: '/PATH/TO/PROJECT/DIR/la_workshop_2025'.
```

List locally available projects

```python
>>> ts.projects()

# Outputs
TidyScreen main database found! Continuing...
Project: la_workshop_2025
         located at /PATH/TO/PROJECT/DIR/la_workshop_2025
```

Importing TidyScreen projects

```python
>>> ts.import_project()

# Outputs
Input the full path to the project to import: # The use should enter the /PATH/TO/PROJECT/DIR
```

Structure of TidyScreen projects

- chemspace
- moldock
- docking_analysis