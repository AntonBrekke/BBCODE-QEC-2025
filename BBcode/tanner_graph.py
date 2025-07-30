import source.tanner_graph_qubit_space as tgq
import source.tanner_graph_group_space as tgg
import source.BBcode_classes as BBcode
import matplotlib.pyplot as plt 

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

code_A312_B312 = BBcode.BBcode_A312_B312(ell, m)
A312, B312 = code_A312_B312.get_AB()

# Visualize qubit-space tanner graph
tgq.make_tanner_graph(A312, B312, ell, m)
# plt.show()
# # Visualize group-space tanner graph
tgg.make_tanner_graph(A312, B312, ell, m)
plt.show()