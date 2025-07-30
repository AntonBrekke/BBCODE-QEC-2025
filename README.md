# Bivariate-Bicycle code
This project implements the Bivariate-Bicycle (BB) code introduced by IBM Quantum in 2024 (see Bravyi et al. https://arxiv.org/pdf/2308.07915).

The implementation mainly uses functionality from the Python library `panqec` (pronounced "pancake"), with additional usage of `ldpc` for decoding and `bposd` to extract logical operators. 

<img width="475.2" height="213.6" alt="image" src="https://github.com/user-attachments/assets/058833f8-9132-45b1-8462-0bf013b55ce0" />

## **Code-structure**
To be able to run every single file, you need the following libraries:

`panqec`, `bposd`, `ldpc`, `tqdm`, `networkx`, `numpy`, `scipy`, `matplotlib`.

### BBcode
- `GUI.py`
    - This code visualizes the BBcodes implemented in the source folder using PanQEC's GUI. To generate the visuals, a .json file needs to be specified. For the `BB2DCode` class,  `_BBcode.json` is utilized and can be tweaked to change e.g. shapes and colors. 

- `tanner_graph.py`
    - Visualizes tanner-graphs of specified BBcode. Imports `tanner_graph_qubit_space.py` and `tanner_graph_group_space.py` to visualize code both in physical and abstract space.

- `analysis.py`
    - Simulates and plots logical vs. physical error-rates for specified codes using tools from `PanQEC`. Resulting figures from this file end up in \figures.

### source
- `BBcode_classes.py` 
    - In this file you can implement your own BBcode by making a subclass of the `BB2DCode` class, which is a subclass of PanQEC's `StabilizerCode` class. This subclass only needs to contain the method `get_AB()`, where monomials A(x,y) and B(x,y) are specified.

- `errormodel_classes.py`
    - In this file you can build your own error-model class as a subclass of PanQEC's `BaseErrorModel` class. Inspired by analog QEC in the GKP qubit, the `GaussianPauliErrorModel` class aims to generate Pauli-errors from displacements in position and momentum drawn from Gaussian distributions. 

- `decoder_classes.py`
    - In this file you can make your own decoder class as a subclass of PanQEC's `BaseDecoder` class. To make use of the analog information generated in the `GaussianPauliErrorModel` class for decoding, the decoder class `BeliefPropagationLSDDecoder` is implemented in this file. This decoder uses Beleif Propagation (BP) with the Localized Statistics Decoding (LSD) algorithm using the `ldpc` library. 

- `tanner_graph_group_space.py`
    - Represents the qubit grid in the abstract group space of the BBcode. Walks in the grid are done by applying group elements given by the A and B monomials. Shows both checks and logical operators for X and Z. 

- `tanner_graph_qubit_space.py`
    - Represents the qubit grid with physical distances. Visualizes which qubits the X-check and Z-check uses for stabilizers. Works nicely with `GUI.py`, where you can highlight corresponding qubits to see that they are indeed stabilizers. 