---
title: Information flow
sidebar_position: 4
---

The basic information flow in TidyScreen is presented in **Figure 1**.

<figure>
  <p align="center">
  <img src="/img/information_flow.png" alt="Description of image" width="1000"/>
  <figcaption>**Figure 1:** Diagram showing the information flow in TidyScreen.</figcaption>
  </p>
</figure>


Required user inputs are shown in dashed boxes (i.e. SMILES and Receptor). Boxes consisting of solid lines represents processes elicited by TidyScreen to process input information.

### Components description

- **SMILES**: a `.csv` file containing the molecules to be included in the virtual screening workflow. Check this input file format in the corresponding section under [ChemSpace actions -> SMILES input file](/chemspace_docs/smiles_input.md)

- **ChemSpace**: the module capable of processing molecules as part of the screening campaign. This processing includes SMILES sanitization, molecule storage and organization. Also, this module is responsible of managing reactant molecules when products to be screened are originated by applying custom chemical reactions. Overall, the whole chemical space associated to the screening is managed by this module. For further detail check the dedicated [ChemSpace actions -> ChemSpace actions](/chemspace_docs/chemspace_actions.md)

- **ChemSpace Analysis**: this module is responsible of processing (and/or sampling) the chemical space associated to the project. This processing includes analyses such as stereoisomers enumeration, filtering by molecular descriptors, prioritization using ML modelos, etc. For implementation details check the dedicated  

- **Receptor**:

- **MolDock**:

- **MolDock Analysis**:

- **MolDyn**:


