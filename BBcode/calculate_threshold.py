import BBcode_classes as BBcode
from BP_LSD import BeliefPropagationLSDDecoder
from panqec.config import CODES, DECODERS
import numpy as np
import os 
from tqdm import tqdm
import matplotlib.pyplot as plt

"""
* Code written by Anton Brekke * 

This file calculates the error threshold given a BBCode class from 'BBcode_classes.py'. 
"""

### Calculate threshold and get error-rate plot 
from panqec.simulation import read_input_dict
from panqec.analysis import Analysis


def calculate_threshold(BBclass: BBcode.AntonBB2DCode=BBcode.BBcode_Toric,
                        decoder_dict: dict = {'name': 'BeliefPropagationLSDDecoder',  #  Class name of the decoder
                                              'parameters': [{'max_bp_iter': 1e3, 'lsd_order': 10, 
                                              'channel_update': False, 'bp_method': 'minimum_sum'}]}, 
                        n_trials: int=1e2, 
                        grids: list[dict]=[{'L_x':10,'L_y':10}],
                        p_range: tuple=(0.1, 0.25, 40),
                        r: tuple=(1/3, 1/3, 1/3)
                        ):

    n_trials = int(n_trials)  # Ensure n_trials is an integer
    p_min, p_max, n_points = p_range
    r_x, r_y, r_z = r
    p = np.linspace(p_min, p_max, n_points)

    # Define which code-class to use 
    # code_class = BBcode.BBcode_Toric
    # code_class = BBcode.BBcode_ArXiV_example
    code_class = BBclass
    code_name = code_class.__name__
    decoder_name = decoder_dict['name']

    # Check if parity checks we implement are the same as in PanQEC 
    test_code = code_class(4,4)
    print(test_code.HX)
    print(test_code.Hx.toarray())
    print(np.all(test_code.HX == test_code.Hx.toarray()))
    print(np.all(test_code.HZ == test_code.Hz.toarray()))

    # Must register the new code in panQEC 
    CODES[f'{code_name}'] = code_class
    DECODERS['BeliefPropagationLSDDecoder'] = BeliefPropagationLSDDecoder

    save_frequency = 10  # Frequency of saving to file
    n_trials_str = f'{n_trials:.0e}'.replace('+0', '')
    grids_str = f'{grids}'.replace(' ', '').replace(':',';').replace(':', ';').replace("'", "")
    p_range_str = f'{p_range}'.replace(' ', '')
    r_tup_str =  f'({r[0]:.1f},' + f'{r[1]:.1f},' + f'{r[2]:.1f})'
    decoder_dict_str = f'{decoder_dict}'.replace(' ', '').replace(':', ';').replace("'", "")
    filename = f"data\{code_name};n_trials;{n_trials_str};grids;{grids_str};p_range;{p_range_str};r;{r_tup_str};decoder;{decoder_dict_str}.json"

    # This magically fixes the fact that the filename is too long... 
    filename = u"\\\\?\\" + os.path.abspath(filename)

    rewrite_data = True
    if os.path.exists(filename):
        advance = False
        while not advance:
            answer = input(f'Filename {filename} already exists. Do you want to write over the existing one (y/n)? ')
            if answer.lower() == 'y':
                advance = True
            elif answer.lower() == 'n':
                rewrite_data = False
                advance = True 

    if rewrite_data:
        input_data = {
            'ranges': {
                'label': 'BB 2D Experiment',  # Can be any name you want
                'code': {
                    'name': f'{code_name}',  # Class name of the code
                    'parameters': grids  # List of dictionaries with code parameters
                },
                'error_model': {
                    'name': 'PauliErrorModel',  #  Class name of the error model
                    'parameters': [
                        {'r_x': r_x, 'r_y': r_y, 'r_z': r_z}  #  Ratios of X, Y and Z errors
                    ],
                },
                'decoder': decoder_dict,
                'error_rate': p.tolist()  #  List of physical error rates
            }
        }

        # If data-file already exists, rewrite the file by setting to 'True'
        if os.path.exists(filename):
            os.remove(filename)

        # We create a BatchSimulation by reading the input dictionary
        batch_sim = read_input_dict(
            input_data,
            output_file=filename  # Where to store the simulation results
            # save_frequency=save_frequency
        ) 

        batch_sim.run(n_trials, progress=tqdm)
    analysis = Analysis(filename)
    return analysis, filename 


def plot_thresholds(analysis, savefig: bool=False, filename: str=None):
    if filename is None:
        filename = 'figures/no_filename.pdf'
    ### Plot resulting data 
    plt.style.use('seaborn-v0_8')
    # Comment back in to get LaTeX font 
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    params = {'axes.labelsize': 14,
            'axes.titlesize': 14,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.title_fontsize': 10,
            'legend.fontsize': 10,
            'font.size': 10} # extend as needed
    # print(plt.rcParams.keys())
    plt.rcParams.update(params)

    custom_cycle = ["#6667AB", "#C74375", "#F0C05A", "#009473", 'red', 'pink']
    plt.rcParams["axes.prop_cycle"] = plt.cycler(color=custom_cycle)

    fig, ax = plt.subplots(ncols=3, figsize=(15, 5))

    plt.sca(ax[0])
    analysis.plot_thresholds()
    plt.sca(ax[1])
    analysis.plot_thresholds(sector='X')
    plt.sca(ax[2])
    analysis.plot_thresholds(sector='Z')

    fig.tight_layout()
    plt.show()

    # analysis.make_collapse_plots()

    fig, ax = plt.subplots(ncols=3, figsize=(15, 5))

    results = analysis.get_results()
    # Count number of combinations of L_x and L_y
    dict_arr = np.array([[*L_dict.values()][:-1] for L_dict in results['code_params']])
    n_Ls = 0 
    com_list = []
    for d in dict_arr:
        lx, ly = d
        if (lx, ly) not in com_list:
            com_list.append((lx, ly))
            n_Ls += 1

    # Divide by total number of choices of L_x x L_y in parameters from input_data to only get each code param. one time 
    n_trials_pr = int(len(results['n_trials'])/n_Ls)
    code_params = results['code_params'][::n_trials_pr]

    capsize = 5
    ms = 5
    legend_lines = []
    legend_labels = []
    for Ls in code_params:
        Lx, Ly, Lz = Ls.values()
        code = eval('BBcode.' + results['code'][0] + f'({Lx}, {Ly})')
        index = results['code_params'] == Ls
        k = code.num_logical_qubits 
        lx = code.lx.toarray()
        w_lx = np.sum(lx, axis=1)
        k = np.max(w_lx)
        p_phys = results[index]['error_rate']
        kp = 1 - (1 - p_phys)**k


        line = ax[0].errorbar(p_phys, results[index]['p_est'], results[index]['p_se'],
                    label=rf'$L_x\!: {Lx}$, $L_y\!: {Ly}$', capsize=capsize, marker='o', ms=ms)
        linecolor = line[0].get_color()
        ax[0].plot(p_phys, kp, color=linecolor, linestyle=(0,(3,6)))
        ax[0].plot(p_phys, kp, color='k', linestyle=(4.5,(3,6)))

        ax[1].errorbar(p_phys, results[index]['p_est_X'], results[index]['p_se_X'],
                    label=rf'$L_x\!: {Lx}$, $L_y\!: {Ly}$', capsize=capsize, marker='o', ms=ms)
        ax[2].errorbar(p_phys, results[index]['p_est_Z'], results[index]['p_se_Z'],
                    label=rf'$L_x\!: {Lx}$, $L_y\!: {Ly}$', capsize=capsize, marker='o', ms=ms)
        
        legend_lines.append(line)
        legend_labels.append(rf'$L_x\!: {Lx}$, $L_y\!: {Ly}$')

    from matplotlib import lines

    th_line1 = lines.Line2D([], [], color='gray', linestyle=(0,(3,6)))
    th_line2 = lines.Line2D([], [], color='k', linestyle=(4.5,(3,6)))

    legend_lines.append((th_line1, th_line2))
    legend_labels.append('pseudo-threshold')

    ax[0].axvline(analysis.thresholds.iloc[0]['p_th_fss'], color='red', linestyle='--')
    ax[0].axvspan(analysis.thresholds.iloc[0]['p_th_fss_left'], analysis.thresholds.iloc[0]['p_th_fss_right'],
                alpha=0.5, color='pink')
    result_X = analysis.sector_thresholds['X']
    result_Z = analysis.sector_thresholds['Z']
    ax[1].axvline(result_X['p_th_fss'][0], color='red', linestyle='--')
    ax[1].axvspan(result_X['p_th_fss_left'][0], result_X['p_th_fss_right'][0],
                alpha=0.5, color='pink')
    ax[2].axvline(result_Z['p_th_fss'][0], color='red', linestyle='--')
    ax[2].axvspan(result_Z['p_th_fss_left'][0], result_Z['p_th_fss_right'][0],
                alpha=0.5, color='pink')

    pth_str_1 = r'$p_{\rm th}' + f' = ({100*analysis.thresholds.iloc[0]["p_th_fss"]:.2f}' + '\pm' + f'{100*analysis.thresholds.iloc[0]["p_th_fss_se"]:.2f})\%$'
    pth_str_2 = r'$p_{\rm th}' + f' = ({100*result_X["p_th_fss"][0]:.2f}' + '\pm' + f'{100*result_X["p_th_fss_se"][0]:.2f})\%$'
    pth_str_3 = r'$p_{\rm th}' + f' = ({100*result_Z["p_th_fss"][0]:.2f}' + '\pm' + f'{100*result_Z["p_th_fss_se"][0]:.2f})\%$'

    # ax[0].set_xscale('log')
    # ax[0].set_yscale('log')
    dist = p_phys.max() - p_phys.min()
    ax[0].set_xlim(p_phys.min()-0.05*dist, p_phys.max()+0.05*dist)
    ax[0].set_ylim(ymax=1.1)

    ax[0].set_title('All errors')
    ax[1].set_title('X errors')
    ax[2].set_title('Z errors')

    ax[0].set_xlabel('Physical error rate')
    ax[0].set_ylabel('Logical error rate')
    ax[0].legend(legend_lines, legend_labels, title=pth_str_1)


    ax[1].set_xlabel('Physical error rate')
    ax[1].set_ylabel('Logical error rate')
    ax[1].legend(title=pth_str_2)

    ax[2].set_xlabel('Physical error rate')
    ax[2].set_ylabel('Logical error rate')
    ax[2].legend(title=pth_str_3)

    fig.tight_layout()
    if savefig: plt.savefig(filename)
    plt.show()

analysis, filename = calculate_threshold(BBclass=BBcode.BBcode_A312_B312,
                    decoder_dict={'name': 'BeliefPropagationOSDDecoder',
                                  'parameters': [{'max_bp_iter': int(1e3), 'osd_order': 10}]},
                    n_trials=5e2, 
                    grids=[{'L_x':6,'L_y':6},
                           {'L_x': 9,'L_y':6},
                           {'L_x': 12,'L_y':6},
                           {'L_x': 15,'L_y':6}],
                    p_range=(0, 0.25, 60),
                    r=(1/3, 1/3, 1/3))

# analysis, filename = calculate_threshold(BBclass=BBcode.BBcode_Toric,
#                     decoder_dict={'name': 'MatchingDecoder',
#                                   'parameters': [{}]},
#                     n_trials=1e2, 
#                     grids=[{'L_x':6,'L_y':6},
#                            {'L_x': 9,'L_y':6},
#                            {'L_x': 12,'L_y':6},
#                            {'L_x': 15,'L_y':6}],
#                     p_range=(0, 0.25, 60),
#                     r=(1/3, 1/3, 1/3))

plot_thresholds(analysis, savefig=True, filename=filename.replace('data', 'figures').replace('.json', '.pdf'))


# print(calculate_threshold.__annotations__)