# Bivariate-Bicycle code
This project implements the Bivariate-Bicycle (BB) code introduced by IBM Quantum in 2024 (see Bravyi et al. https://arxiv.org/pdf/2308.07915) in Python.

The implementation mainly uses functionality from the Python library `panqec` (pronounced "pancake"), with additional usage of `ldpc` for decoding and `bposd` to extract logical operators. 

The framework for running simulations on codes should be fairly general. To implement your own code, follow the recipe from `BBCode_classes.py`, where you make your own motherclass as a subclass of the `StabilizerCode` class from `panqec`. The only thing you should need to specify is how the check-matrices are built. 

<img width="475.2" height="213.6" alt="image" src="https://github.com/user-attachments/assets/058833f8-9132-45b1-8462-0bf013b55ce0" />

## **Code-structure**
To be able to run every single file, you need the following libraries:

`panqec`, `bposd`, `ldpc`, `tqdm`, `networkx`, `numpy`, `scipy`, `matplotlib`.

### BBcode
- `GUI.py`
    - This code visualizes the BBcodes implemented in the source folder using PanQEC's GUI. To generate the visuals, a .json file needs to be specified. For the `BB2DCode` class,  `BBcode.json` is utilized and can be tweaked to change e.g. shapes and colors. 

- `tanner_graph.py`
    - Visualizes tanner-graphs of specified BBcode. Imports `tanner_graph_qubit_space.py` and `tanner_graph_group_space.py` to visualize code both in physical and abstract space.

- `make_figs.py`
    - Imports `analysis.py` from source to run simulations and make figures. Mainly to separate source code from interface. 

### js
- `BBcode.json`
    - Specifies how the code should render in the GUI from `GUI.py` for BBcode. Specify shapes, colors etc.

- `TileCode.json`
    - Specifies how the code should render in the GUI from `GUI.py` for TileCode. Specify shapes, colors etc.

### source
- `BBcode_classes.py` 
    - In this file you can implement your own BBcode by making a subclass of the `BB2DCode` class, which is a subclass of PanQEC's `StabilizerCode` class. This subclass only needs to contain the method `get_AB()`, where monomials A(x,y) and B(x,y) are specified.

- `TileCode_classes.py` 
    - In this file you can implement your own TileCode. Original code belongs to https://github.com/Bragit123/QEC, read about the code here for more information. 

- `errormodel_classes.py`
    - In this file you can build your own error-model class as a subclass of PanQEC's `BaseErrorModel` class. Inspired by analog QEC in the GKP qubit, the `GaussianPauliErrorModel` class aims to generate Pauli-errors from displacements in position and momentum drawn from Gaussian distributions. 

- `decoder_classes.py`
    - In this file you can make your own decoder class as a subclass of PanQEC's `BaseDecoder` class. To make use of the analog information generated in the `GaussianPauliErrorModel` class for decoding, the decoder class `BeliefPropagationLSDDecoder` is implemented in this file. This decoder uses Beleif Propagation (BP) with the Localized Statistics Decoding (LSD) algorithm using the `ldpc` library. 

- `analysis.py`
    - Simulates and plots logical vs. physical error-rates for specified codes using tools from `PanQEC`. Resulting figures from this file end up in \figures.

- `tanner_graph_group_space.py`
    - Represents the qubit grid in the abstract group space of the BBcode. Walks in the grid are done by applying group elements given by the A and B monomials. Shows both checks and logical operators for X and Z. 

- `tanner_graph_qubit_space.py`
    - Represents the qubit grid with physical distances. Visualizes which qubits the X-check and Z-check uses for stabilizers. Works nicely with `GUI.py`, where you can highlight corresponding qubits to see that they are indeed stabilizers. 