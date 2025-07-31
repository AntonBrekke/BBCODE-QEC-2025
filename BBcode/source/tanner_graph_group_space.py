import numpy as np
import matplotlib.pyplot as plt
from bposd.css import css_code
from matplotlib.patches import Rectangle, Circle, Patch
import matplotlib.patches as mpatches
from matplotlib.widgets import Button
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def draw_bivariate_tanner_from_checks(fig, ax,
                                     hx: np.ndarray,
                                     hz: np.ndarray,
                                     ell: int,
                                     m: int,
                                     q: int,
                                     field: int,
                                     name: str=None):
    """
    Draw the Bivariate Bicycle Tanner graph using parity-check matrices
    hx = [A | B], hz = [B^T | A^T], replicating the original A/B polynomial-based plot.

    hx, hz : shape (ell*m, 2*ell*m)
    ell, m   : code parameters (so ell*m = n/2)
    q      : field characteristic multiplier for Z-edge weights
    field  : size of finite field for dash computation
    name   : optional title
    """

    short_offsets = {( 1,  0),  # one step +x
                     (0,  0),  # no steps (for both x and y)
                     ( 0,  1),  # one step +y
                    }

    # Dimensions
    lm, n = hx.shape
    assert n == 2*lm and lm == ell*m, "Invalid shapes: require hx shape (ell*m,2*ell*m)"
    n2 = lm  # number of columns in A or B block
    nx_cells = lm // m  # = ell
    ny_cells = lm // ell   # = m

    # 1) Extract A and B coefficient arrays from hx row 0
    """
    a_coeff = Coeff{[[x^0y^0, x^0y^1, x^0y^2, ... x^0y^(m-1)], 
               [x^1y^0, x^1y^1, x^1y^2, ..., x^1y^(m-1)], 
               ..., 
               [x^(ell-1)y^0, x^(ell-1)y^1, ..., x^(ell-1)y^(m-1)]]}
    e.g. A(x, y) = x + x^2 + y^3, B(x, y) = x^3 + y + y^2
    a_coeff = [[0,0,0,1, ...], 
               [1,0,0,0, ...], 
               [1,0,0,0, ...], ... ] 
    a_coeff = [[0,1,1, ...], 
               [0,0,0, ...], 
               [0,0,0, ...], 
               [1,0,0, ...], ... ] 
    where ... represents 0's filling in the remaining shape.
    """
    a_coefficients = np.zeros((ell, m), dtype=int)
    b_coefficients = np.zeros((ell, m), dtype=int)
    row0 = hx[0, :]
    for c in range(n2):
        w = int(row0[c])
        if w:
            i = c // m
            j = c %  m
            a_coefficients[i, j] = w
    for c in range(n2, 2*n2):
        w = int(row0[c])
        if w:
            idx = c - n2
            i = idx // m
            j = idx %  m
            b_coefficients[i, j] = w

    # 2) Compute min/max degrees for A and B
    def factor(coeff):
        "Find index of smallest and largest (k, l) for x^k*y^l with non-zero coeff in polynomial"
        nz = np.array(coeff.nonzero()).T  # list of (k,l)
        deg = nz.sum(axis=1)    # deg = k + l 
        idx_min = np.argmin(deg)        # where nz has smallest k + l
        idx_max = np.argmax(deg)        # where nz har largest k + l
        # Return (k,l) such that k+l is smallest, (k,l) such that k+l is largest
        return tuple(nz[idx_min, :]), tuple(nz[idx_max, :])
    
    (a_k_min, a_l_min), (a_k_max, a_l_max) = factor(a_coefficients)
    (b_k_min, b_l_min), (b_k_max, b_l_max) = factor(b_coefficients)
    x_max = max(a_k_max, b_k_max)   # Largest x-power in (A, B)
    y_max = max(a_l_max, b_l_max)   # largest y-power in (A, B)

    # 3) Set up plot
    # fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(-0.3, nx_cells-0.2)
    ax.set_ylim(-0.3, ny_cells-0.2)
    ax.set_aspect('equal', adjustable='box')
    ax.set_axis_off()
    if name:
        ax.set_title(f'Tanner Graph of {name}')

    # Node shapes
    def x_stabiliser(x, y): 
        return Rectangle((x, y), 0.1, 0.1, edgecolor='lightcoral', facecolor='lightcoral', zorder=3)
    def z_stabiliser(x, y): 
        return Rectangle((x, y), 0.1, 0.1, edgecolor="#47BE6B", facecolor='#47BE6B', zorder=3)
    def l_data(x, y): 
        # bLue
        return Circle((x, y), 0.06, edgecolor='royalblue', facecolor='royalblue', zorder=3)
    def r_data(x, y): 
        # oRange
        return Circle((x, y), 0.06, edgecolor="#F2C531", facecolor='#F2C531', zorder=3)

    # Draw nodes in grid
    for i in range(nx_cells):
        for j in range(ny_cells):
            # ax.add_patch(x_stabiliser(i+0.5-0.05, j-0.05))
            # ax.add_patch(z_stabiliser(i-0.05, j+0.5-0.05))
            # ax.add_patch(l_data(i+0.5, j+0.5))
            # ax.add_patch(r_data(i, j))
            ax.add_patch(x_stabiliser(i+0.5-0.05, j+0.5-0.05))
            ax.add_patch(z_stabiliser(i-0.05, j-0.05))
            ax.add_patch(l_data(i+0.5, j))
            ax.add_patch(r_data(i, j+0.5))
            ax.axvline(i, color='gray', alpha=0.05, zorder=-1)
            ax.axhline(j+0.5, color='gray', alpha=0.05, zorder=-1)
            ax.axvline(i+0.5, color='gray', alpha=0.05, zorder=-1)
            ax.axhline(j, color='gray', alpha=0.05, zorder=-1)

    # Style helper
    def style_x(w):
        if w == 1:
            return {'color':'k','dashes':None}
        dash = 16/(w**2)
        return {'color':'k','dashes':[dash, 2, dash, 2]}
    def style_z(w):
        div = (q * w) % field
        if div == 1:
            return {'color':'k','dashes':None}
        dash = 16/(div**2) if div else 4
        return {'color':'k','dashes':[dash, 2, dash, 2]}

    # 4) Draw X-edges from A and B supports
    # Add margins to x_max, y_max to make sure support is drawn
    for i in range(-x_max, x_max + nx_cells):
        for j in range(-y_max, y_max + ny_cells):
            # A-block X edges
            for k in range(ell):        # x^k, k in [0, ell-1]
                for l in range(m):        # y^l, l in [0, m-1]
                    w = a_coefficients[k, l]
                    if not w: continue
                    dx = k - a_k_min
                    dy = l - a_l_min
                    style = style_x(w)
                    if i == 5 and j == 3:
                        if (dx, dy) in short_offsets or True:
                            line, = ax.plot([0.5+i, 0.5+i+dx], [0.5+j, 1+j-dy], color=style['color'])
                            if style['dashes']: line.set_dashes(style['dashes'])
                        else: 
                            lr = mpatches.FancyArrowPatch([0.5+i, j], [0.5+i+dx, 0.5+j-dy], connectionstyle="arc3,rad=-0.1", color='k', lw=1.5)
                            plt.gca().add_patch(lr)
            # B-block X edges
            for k in range(ell):
                for l in range(m):
                    w = b_coefficients[k, l]
                    if not w: continue
                    dx = k - b_k_min
                    dy = l - b_l_min
                    style = style_x(w)
                    if i == 5 and j == 3:
                        if (dx, dy) in short_offsets or True:
                            line, = ax.plot([0.5+i, i+dx], [0.5+j, 0.5+j-dy], color=style['color'])
                            if style['dashes']: line.set_dashes(style['dashes'])
                        else:
                            lr = mpatches.FancyArrowPatch([0.5+i, 0.5+j], [i+dx, 0.5+j-dy], connectionstyle="arc3,rad=0.1", color='k', lw=1.5)
                            plt.gca().add_patch(lr)

    # 5) Draw Z-edges from A and B supports
    for i in range(-x_max, x_max + nx_cells):
        for j in range(-y_max, y_max + ny_cells):
            # A-block Z edges
            for k in range(ell):
                for l in range(m):
                    w = a_coefficients[k, l]
                    if not w: continue
                    dx = k - a_k_min
                    dy = l - a_l_min
                    style = style_z(w)
                    if i == 7 and j == 2:
                        if (dx, dy) in short_offsets or True:
                            line, = ax.plot([i, i-dx], [j, -0.5+j+dy], color=style['color'])
                            if style['dashes']: line.set_dashes(style['dashes'])
                        else:
                            lr = mpatches.FancyArrowPatch([i, j], [i-dx, -0.5+j+dy], connectionstyle="arc3,rad=-0.1", color='k', lw=1.5)
                            plt.gca().add_patch(lr)
            # B-block Z edges
            for k in range(ell):
                for l in range(m):
                    w = b_coefficients[k, l]
                    if not w: continue
                    dx = k - b_k_min
                    dy = l - b_l_min
                    style = style_z(w)
                    if i == 7 and j == 2:
                        if (dx, dy) in short_offsets or True:
                            line, = ax.plot([i, 0.5+i-dx], [j, j+dy], color=style['color'])
                            if style['dashes']: line.set_dashes(style['dashes'])
                        else:
                            lr = mpatches.FancyArrowPatch([i, 0.5+j], [0.5+i-dx, 0.5+j+dy], connectionstyle="arc3,rad=0.1", color='k', lw=1.5)
                            plt.gca().add_patch(lr)

    # 6) Draw boundary arrows
    ax.plot([-0.25, -0.25], [-0.25, ny_cells-0.25], color='black', linewidth=0.7)
    ax.arrow(-0.25, -0.25, 0, ny_cells/2, head_width=0.1, head_length=0.1,
             color='black', linewidth=0.05)
    ax.plot([-0.25, nx_cells-0.25], [-0.25, -0.25], color='black', linewidth=0.7)
    ax.arrow(-0.25, -0.25, (nx_cells)/2-0.05, 0, head_width=0.1, head_length=0.1,
             color='black', linewidth=0.05)
    ax.arrow(-0.25, -0.25, (nx_cells)/2+0.05, 0, head_width=0.1, head_length=0.1,
             color='black', linewidth=0.05)
    ax.plot([-0.25, nx_cells-0.25], [ny_cells-0.25, ny_cells-0.25], color='black', linewidth=0.7)
    ax.arrow(-0.25, ny_cells-0.25, (nx_cells)/2-0.05, 0, head_width=0.1, head_length=0.1,
             color='black', linewidth=0.05)
    ax.arrow(-0.25, ny_cells-0.25, (nx_cells)/2+0.05, 0, head_width=0.1, head_length=0.1,
             color='black', linewidth=0.05)
    ax.plot([nx_cells-0.25, nx_cells-0.25], [-0.25, ny_cells-0.25], color='black', linewidth=0.7)
    ax.arrow(nx_cells-0.25, -0.25, 0, ny_cells/2, head_width=0.1, head_length=0.1,
             color='black', linewidth=0.05)

    # 7) Legend
    handles = [Patch(color='lightcoral'), Patch(color='#47BE6B'),
               Patch(color='royalblue'), Patch(color="#F2C531")]
    labels = ['X stabiliser','Z stabiliser','Left data','Right data']
    # add weight dashes
    for i in range(1, field):
        style = style_x(i)
        # xline = mpatches.Patch(edgecolor=style['color'], linestyle='--' if style['dashes'] else '-', label=f'$X^{i}$')
        xline, = plt.plot([0], [0], color=style['color'])
        handles.append(xline)
        style = style_z(i)
        # zline = mpatches.Patch(edgecolor=style['color'], linestyle='--' if style['dashes'] else '-', label=f'$Z^{i}$')
        zline, = plt.plot([0], [0], color=style['color'])
        handles.append(zline)
        labels.extend([f'$X^{i}$', f'$Z^{i}$'])
    ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1,1), handlelength=2.4)

def plot_logicals(ax, hx, hz, ell, m, x_logical_idx=0, z_logical_idx=0, show_x=True, show_z=True):

    n2, n = hx.shape

    qcode = css_code(hx, hz)
    # Logical check matrices -- for logical operators later
    logicals_x = qcode.lx.toarray()
    logicals_z = qcode.lz.toarray()

    # 1) Extract A/B supports from H_X row 0
    a0 = hx[0, :n2]
    b0 = hx[0, n2:]
    """
    A_support = {(k,l) | x^k*y^l in A(x,y)} = [(k,l) for k in range(ell) for l in range(m) if a_coefficients[k,l]]
    Sneaky way to do so: A_support = [(c//m, c%m) for c,w in enumerate(a0) if w]
    """
    A_support = [(c//m, c%m) for c,w in enumerate(a0) if w]
    B_support = [(c//m, c%m) for c,w in enumerate(b0) if w]
    # include inverses so we catch loops in both directions
    # (-k) % N = (N-k) % N 
    invA = [((-dx)%ell,(-dy)%m) for dx,dy in A_support]
    invB = [((-dx)%ell,(-dy)%m) for dx,dy in B_support]
    full_X_support = A_support + B_support + invA + invB

    # 2) Extract B/A supports from H_Z row 0 for the X‐graph
    #    H_Z = [B^T | A^T], so left n2 of hz row0 is B^T support,
    #    right n2 is A^T support.
    b0T = hz[0, :n2]
    a0T = hz[0, n2:]
    BT_support = [(c//m, c%m) for c,w in enumerate(b0T) if w]
    AT_support = [(c//m, c%m) for c,w in enumerate(a0T) if w]
    invBT = [((-dx)%ell,(-dy)%m) for dx,dy in BT_support]
    invAT = [((-dx)%ell,(-dy)%m) for dx,dy in AT_support]
    full_Z_support = BT_support + AT_support + invBT + invAT

    # Define colors for logical x and logical z
    col_log_x = 'red'
    col_log_z = 'blue'
    # 4) Overlay *all* logical‐Z loops (in ker H_X) using full_X_support
    for vec in logicals_z:
        if not np.all(vec == logicals_z[z_logical_idx]) or not show_z:
            continue
        # support = set of qubit‐indices in first n2
        idxs = np.nonzero(vec[:n2])[0]
        Z_nodes = {(idx//m, idx%m) for idx in idxs}  # set of (i,j), ith column, jth row of grid
        for (i,j) in Z_nodes:
            for dx,dy in full_X_support:
                ni, nj = (i+dx)%ell, (j+dy)%m
                if (ni,nj) in Z_nodes:
                    xi, yj = i, j 
                    if (i+dx)>ell-0.5:
                        ax.arrow(i+0.5, j, dx, dy, head_width=0, head_length=0, fc=col_log_z, ec=col_log_z, linewidth=2, alpha=0.8, length_includes_head=True)
                        xi = (i+dx)%ell - dx
                    if (j+dy)>m-0.5:
                        ax.arrow(i+0.5, j, dx, dy, head_width=0, head_length=0, fc=col_log_z, ec=col_log_z, linewidth=2, alpha=0.8, length_includes_head=True)
                        yj = (j+dy)%m - dy
                    ax.arrow(xi+0.5, yj, dx, dy, head_width=0.1, head_length=0.1, fc=col_log_z, ec=col_log_z, linewidth=2, alpha=0.8, length_includes_head=True)

        idxs = np.nonzero(vec[n2:])[0]
        Z_nodes = {(idx//m, idx%m) for idx in idxs}  # set of (i,j)
        for (i,j) in Z_nodes:
            for dx,dy in full_X_support:
                ni, nj = (i+dx)%ell, (j+dy)%m
                if (ni,nj) in Z_nodes:
                    xi, yj = i, j
                    if (i+dx)>ell-0.5:
                        ax.arrow(i, j+0.5, dx, dy, head_width=0, head_length=0, fc=col_log_z, ec=col_log_z, linewidth=2, alpha=0.8, length_includes_head=True)
                        xi = (i+dx)%ell - dx
                    if (j+dy)>m-0.5:
                        ax.arrow(i, j+0.5, dx, dy, head_width=0, head_length=0, fc=col_log_z, ec=col_log_z, linewidth=2, alpha=0.8, length_includes_head=True)
                        yj = (j+dy)%m - dy
                    ax.arrow(xi, yj+0.5, dx, dy, head_width=0.1, head_length=0.1, fc=col_log_z, ec=col_log_z, linewidth=2, alpha=0.8, length_includes_head=True)

    # 4) Overlay *all* logical‐X loops (in ker H_X) using full_X_support
    for vec in logicals_x:
        if not np.all(vec == logicals_x[x_logical_idx]) or not show_x:
            continue
        # support = set of qubit‐indices in first n2
        idxs = np.nonzero(vec[:n2])[0]
        X_nodes = {(idx//m, idx%m) for idx in idxs}  # set of (i,j)
        for (i,j) in X_nodes:
            for dx,dy in full_Z_support:
                ni, nj = (i+dx)%ell, (j+dy)%m
                if (ni,nj) in X_nodes:
                    xi, yj = i, j
                    if (i+dx)>ell-0.5:
                        ax.arrow(i+0.5, j, dx, dy, head_width=0, head_length=0, fc=col_log_x, ec=col_log_x, linewidth=2, alpha=0.8, length_includes_head=True)
                        xi = (i+dx)%ell - dx
                    if (j+dy)>m-0.5:
                        ax.arrow(i+0.5, j, dx, dy, head_width=0, head_length=0, fc=col_log_x, ec=col_log_x, linewidth=2, alpha=0.8, length_includes_head=True)
                        yj = (j+dy)%m - dy
                    ax.arrow(xi+0.5, yj, dx, dy, head_width=0.1, head_length=0.1, fc=col_log_x, ec=col_log_x, linewidth=2, alpha=0.8, length_includes_head=True)

        idxs = np.nonzero(vec[n2:])[0]
        X_nodes = {(idx//m, idx%m) for idx in idxs}  # set of (i,j), ith column, jth row
        for (i,j) in X_nodes:
            for dx,dy in full_Z_support:
                ni, nj = (i+dx)%ell, (j+dy)%m
                if (ni,nj) in X_nodes:
                    xi, yj = i, j
                    if (i+dx)>ell-0.5:
                        ax.arrow(i, j+0.5, dx, dy, head_width=0, head_length=0, fc=col_log_x, ec=col_log_x, linewidth=2, alpha=0.8, length_includes_head=True)
                        xi = (i+dx)%ell - dx
                    if (j+dy)>m-0.5:
                        ax.arrow(i, j+0.5, dx, dy, head_width=0, head_length=0, fc=col_log_x, ec=col_log_x, linewidth=2, alpha=0.8, length_includes_head=True)
                        yj = (j+dy)%m - dy
                    ax.arrow(xi, yj+0.5, dx, dy, head_width=0.1, head_length=0.1, fc=col_log_x, ec=col_log_x, linewidth=2, alpha=0.8, length_includes_head=True)



def make_tanner_graph(A, B, ell, m):
    AT = np.transpose(A)
    BT = np.transpose(B)

    # Each row of HX defines an X-type check operator X(v) 
    HX = np.hstack((A,B))
    # Each row of HZ defines a Z-type check operator Z(v)
    HZ = np.hstack((BT,AT))

    from bposd.css import css_code
    qcode=css_code(HX, HZ)
    qcode = qcode
    # Logical check matrices -- for logical operators later
    logicals_x = qcode.lx
    logicals_z = qcode.lz

    # draw_bivariate_tanner_from_checks(HX, HZ, ell, m, 1, 2)

    num_logicals = logicals_x.shape[0]

    fig = plt.figure(figsize=(ell if ell>5 else 5, m if m>5 else 5), constrained_layout=True)
    ax = fig.add_subplot()

    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    # plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    # params = {'axes.labelsize': 14,
    #         'axes.titlesize': 14,
    #         'xtick.labelsize': 10,
    #         'ytick.labelsize': 10,
    #         'legend.title_fontsize': 10,
    #         'legend.fontsize': 10,
    #         'font.size': 10,
    #         'figure.titlesize': 16} # extend as needed
    # # print(plt.rcParams.keys())
    # plt.rcParams.update(params)

    ax_size = ax.get_window_extent()._bbox
    x0 = ax_size.x0
    x1 = ax_size.x1
    y0 = ax_size.y0
    y1 = ax_size.y1

    class Index:
        ind_x = 0
        ind_z = 0
        show_x = 1
        show_z = 1

        def next_x(self, event):
            ax.cla()
            if event.button == 1: self.ind_x += 1       # Left-click
            if event.button == 3: self.ind_x -= 1       # Right-click
            i_x = self.ind_x % num_logicals
            i_z = self.ind_z % num_logicals
            # This line makes bnext_x global, otherwise it is deleted and Button will not work when imported
            bnext_x.label.set_text(f'Next X: {i_x+1}')
            draw_bivariate_tanner_from_checks(fig, ax, HX, HZ, ell, m, 1, 2)
            plot_logicals(ax, HX, HZ, ell, m, i_x, i_z, self.show_x, self.show_z)
            plt.draw()

        def next_z(self, event):
            ax.cla()
            if event.button == 1: self.ind_z += 1       # Left-click
            if event.button == 3: self.ind_z -= 1       # Right-click
            i_x = self.ind_x % num_logicals
            i_z = self.ind_z % num_logicals
            # This line makes bnext_z global, otherwise it is deleted and Button will not work when imported
            bnext_z.label.set_text(f'Next Z: {i_z+1}')
            draw_bivariate_tanner_from_checks(fig, ax, HX, HZ, ell, m, 1, 2)
            plot_logicals(ax, HX, HZ, ell, m, i_x, i_z, self.show_x, self.show_z)
            # ax.text(ell, y1, f'Z: {self.ind_z}', fontsize=14)
            plt.draw()

        def toggle_x(self, event):
            ax.cla()
            self.show_x += 1
            self.show_x = self.show_x % 2
            if self.show_x: state = 'on'
            else: state = 'off'
            i_x = self.ind_x % num_logicals
            i_z = self.ind_z % num_logicals
            # This line makes btoggle_x global, otherwise it is deleted and Button will not work when imported
            btoggle_x.label.set_text(f'Toggle X: {state}')
            draw_bivariate_tanner_from_checks(fig, ax, HX, HZ, ell, m, 1, 2)
            plot_logicals(ax, HX, HZ, ell, m, i_x, i_z, self.show_x, self.show_z)
            plt.draw()

        def toggle_z(self, event):
            ax.cla()
            self.show_z += 1
            self.show_z = self.show_z % 2 
            if self.show_z: state = 'on'
            else: state = 'off'
            i_x = self.ind_x % num_logicals
            i_z = self.ind_z % num_logicals
            # This line makes btoggle_x global, otherwise it is deleted and Button will not work when imported
            btoggle_z.label.set_text(f'Toggle Z: {state}')
            draw_bivariate_tanner_from_checks(fig, ax, HX, HZ, ell, m, 1, 2)
            plot_logicals(ax, HX, HZ, ell, m, i_x, i_z, self.show_x, self.show_z)
            plt.draw()



    draw_bivariate_tanner_from_checks(fig, ax, hx=HX, hz=HZ, ell=ell,m= m, q=1, field=2)
    plot_logicals(ax, hx=HX, hz=HZ, ell=ell, m=m, x_logical_idx=0, z_logical_idx=0)

    leg = ax.get_legend().get_window_extent().transformed(fig.transFigure.inverted())
    leg_y0 = leg.y0
    leg_x0 = leg.x0

    callback = Index()
    b_width = 0.8
    b_height = 0.3

    ax_xnext = inset_axes(ax, width=b_width, height=b_height, bbox_to_anchor=(0.975, 0.5), bbox_transform=fig.transFigure)
    ax_znext = inset_axes(ax, width=b_width, height=b_height, bbox_to_anchor=(0.975, 0.4), bbox_transform=fig.transFigure)
    ax_xtoggle = inset_axes(ax, width=b_width, height=b_height, bbox_to_anchor=(0.975, 0.3), bbox_transform=fig.transFigure)
    ax_ztoggle = inset_axes(ax, width=b_width, height=b_height, bbox_to_anchor=(0.975, 0.2), bbox_transform=fig.transFigure)
    # ax_xnext = fig.add_axes([x1, leg_y0-0.25*b_height, b_width, b_height])
    # ax_znext = fig.add_axes([x1, leg_y0-1.25*b_height, b_width, b_height])
    # ax_xtoggle = fig.add_axes([x1, leg_y0-2.25*b_height, b_width, b_height])
    # ax_ztoggle = fig.add_axes([x1, leg_y0-3.25*b_height, b_width, b_height])
    bnext_x = Button(ax_xnext, f'Next X: {callback.ind_x+1}')
    bnext_z = Button(ax_znext, f'Next Z: {callback.ind_z+1}')
    btoggle_x = Button(ax_xtoggle, 'Toggle X: on')
    btoggle_z = Button(ax_ztoggle, 'Toggle Z: on')
    bnext_x.on_clicked(callback.next_x)
    bnext_z.on_clicked(callback.next_z)
    btoggle_x.on_clicked(callback.toggle_x)
    btoggle_z.on_clicked(callback.toggle_z)

    b_fontsize = 8
    bnext_x.label.set_fontsize(b_fontsize)
    bnext_z.label.set_fontsize(b_fontsize)
    btoggle_x.label.set_fontsize(b_fontsize)
    btoggle_z.label.set_fontsize(b_fontsize)
                    
    # fig.tight_layout()
    # plt.show()