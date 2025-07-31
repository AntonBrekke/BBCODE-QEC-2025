import source.tanner_graph_qubit_space as tgq
import source.tanner_graph_group_space as tgg
import source.BBcode_classes as BBcode
import matplotlib.pyplot as plt 
import numpy as np

ell, m = 12, 6

# Make test on Toric surface code 
"""
code_toric = BBcode.BBcode_Toric(ell, m)
A_toric, B_toric = code_toric.get_AB()

# Visualize qubit-space tanner graph
tgq.make_tanner_graph(A_toric, B_toric, ell, m)
plt.show()
# Visualize group-space tanner graph
tgg.make_tanner_graph(A_toric, B_toric, ell, m)
plt.show()
"""

# If you want more direct control
I_ell = np.identity(ell, dtype=int)
I_m = np.identity(m, dtype=int)
I = np.identity(ell*m, dtype=int) 
# Neat way to make polynomial powers of cycle matrix
x = {}
y = {}
# x =  S_ell \otimes I_m
for i in range(ell*m):
    x[i] = np.kron(np.roll(I_ell, i, axis=1), I_m)
# y =  I_ell \otimes S_m
for i in range(ell*m):
    y[i] = np.kron(I_ell, np.roll(I_m, i, axis=1))

A = (I + y[1])%2
B = (I + x[1])%2

# Or if you want to get a specific code 
code_A312_B312 = BBcode.BBcode_A312_B312(ell, m)
A312, B312 = code_A312_B312.get_AB()

# Visualize qubit-space tanner graph
# tgq.make_tanner_graph(A312, B312, ell, m)
# plt.show()
# # Visualize group-space tanner graph
tgg.make_tanner_graph(A312, B312, ell, m)
plt.show()