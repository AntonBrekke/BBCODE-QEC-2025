import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

"""
Index space refers to the space when each qubit is given an index, 
counting from left to right and top to bottom. This is the convention used in the PanQEC library.
"""

def make_tanner_graph(A, B, ell, m):
    # Create an empty graph
    G = nx.Graph()

    # code length
    n = 2*m*ell
    n2 = m*ell

    AT = np.transpose(A)
    BT = np.transpose(B)

    # Each row of HX defines an X-type check operator X(v) 
    HX = np.hstack((A,B))
    # Each row of HZ defines a Z-type check operator Z(v)
    HZ = np.hstack((BT,AT))
    zero_matrix = np.zeros_like(HX)
    H0 = np.vstack((HX, zero_matrix))
    H1 = np.vstack((zero_matrix, HZ))
    H = np.hstack((H0, H1))

    q_nodes = [f'q{i+1}' for i in range(len(HX[0,:]))]
    x_nodes = [f'X{i+1}' for i in range(len(HX[:,0]))]
    z_nodes = [f'Z{i+1}' for i in range(len(HZ[:,0]))]

    # Add nodes with a 'bipartite' attribute to distinguish sets
    G.add_nodes_from(q_nodes, bipartite=0) # Set 0 for variable nodes
    G.add_nodes_from(x_nodes, bipartite=1) # Set 1 for check nodes
    G.add_nodes_from(z_nodes, bipartite=1) # Set 1 for check nodes

    # Add edges connecting variable nodes to check nodes
    # These represent the connections in the parity-check matrix
    draw_edgesX = [(f'q{i+1}', f'X{j+1}', {'color': "#9500FF" if i<ell*m else "#FFA200"}) for j in range(len(HX[:,0])) for i in HX[j,:].nonzero()[0]]
    draw_edgesZ = [(f'q{i+1}', f'Z{j+1}', {'color': "#9500FF" if i<ell*m else "#FFA200"}) for j in range(len(HZ[:,0])) for i in HZ[j,:].nonzero()[0]]


    vertex_num_X = ell + 2 
    # 6 is the weight of the checks
    weight = np.sum(HX[0,:])
    G.add_edges_from(draw_edgesX[weight*(vertex_num_X-1):weight*vertex_num_X])
    vertex_num_Z = n2+1 - (ell + 2)
    G.add_edges_from(draw_edgesZ[weight*(vertex_num_Z-1):weight*vertex_num_Z])

    # Define positions for the nodes for better visualization
    # You can adjust these manually or use a layout algorithm
    pos = {}
    # Position variable nodes in a line

    """
    PanQEC convention for labelling and positioning 
    """
    for i, node in enumerate(q_nodes):
        if i < ell*m:
            pos[node] = (i//m+0.5, i%m)
        else:
            pos[node] = (i//m-ell, i%m+0.5)
    # Position check nodes in another line
    vertex = []
    for i, node in enumerate(z_nodes):
        pos[node] = (i//m, i%m) # Offset for better alignment
        vertex.append(pos[node])
    for i, node in enumerate(x_nodes):
        pos[node] = (i//m+0.5, i%m+0.5) # Offset for better alignment
        # vertex.append(pos[node])
    """
    Article convention 
    Is the "transposed" of the PanQEC convention
    A and B switch roles <--> interchange A and B
    E.g. article convention with 
    A = x[4] + y[1] + y[2]
    B = y[3] + x[1] + x[2]
    is the same as PanQEC convention with
    A = y[3] + x[1] + x[2]
    B = x[4] + y[1] + y[2]
    Or maybe even just (ell, m, x, y) <--> (m, ell, y, x) wihtout changing A, B.
    """
    # for i, node in enumerate(q_nodes):
    #     if i < ell*m:
    #         pos[node] = (i//m, i%m+0.5)
    #     else:
    #         pos[node] = (i//m-ell+0.5, i%m)
    # # Position check nodes in another line
    # vertex = []
    # for i, node in enumerate(z_nodes):
    #     pos[node] = (i//m, i%m) # Offset for better alignment
    #     vertex.append(pos[node])
    # for i, node in enumerate(x_nodes):
    #     pos[node] = (i//m+0.5, i%m+0.5) # Offset for better alignment
    #     # vertex.append(pos[node])

    edges = G.edges()
    edges_colors = [G[u][v]['color'] for u,v in edges]

    q_nodes_colors = ["#2C7FE4" if i<n2 else "#F4C12A" for i in range(n)]
    # Draw the graph

    fig = plt.figure(figsize=(ell, m))
    for i in range(len(vertex)):
        plt.axvline(vertex[i][0], color='grey', alpha=0.1, zorder=0)
        plt.axhline(vertex[i][1], color='grey', alpha=0.1, zorder=0)
    nx.draw_networkx_nodes(G, pos, nodelist=q_nodes, node_color=q_nodes_colors, alpha=0.85, node_size=100)
    nx.draw_networkx_nodes(G, pos, nodelist=x_nodes, node_color="#F95F5F", alpha=0.85, node_size=100, node_shape='s')
    nx.draw_networkx_nodes(G, pos, nodelist=z_nodes, node_color="#37C656", alpha=0.85, node_size=100, node_shape='s')
    nx.draw_networkx_edges(G, pos, edge_color=edges_colors, width=2.5, connectionstyle="arc3,rad=0.11", arrows=True, arrowstyle='-')
    # Uncomment to get labelling/indexing of qubit grid
    nx.draw_networkx_labels(G, pos, font_size=5, font_color='black')


    # plt.title("Tanner Graph Example")
    # plt.legend()
    plt.axis('off') # Hide axes
    plt.tight_layout()
    # plt.show()