import BBcode_classes as BBcode
from decoder_classes import BeliefPropagationLSDDecoder
from errormodel_classes import GaussianPauliErrorModel

from panqec.config import CODES, DECODERS, ERROR_MODELS
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


def simulate_code(BBclass: BBcode.AntonBB2DCode=BBcode.BBcode_Toric,
                  error_model_dict: dict = {'name': 'GaussianPauliErrorModel',  #  Class name of the error model
                                                  'parameters': [{'r_x': 1/3, 'r_y': 1/3, 'r_z': 1/3}]},
                  decoder_dict: dict = {'name': 'BeliefPropagationLSDDecoder',  #  Class name of the decoder
                                              'parameters': [{'max_bp_iter': 1e3, 'lsd_order': 10, 
                                              'channel_update': False, 'bp_method': 'minimum_sum'}]}, 
                  n_trials: int=1e2, 
                  grids: list[dict]=[{'L_x':10,'L_y':10}],
                  p_range: tuple=(0.1, 0.25, 40),
                  ask_overwrite: bool=True):

    n_trials = int(n_trials)  # Ensure n_trials is an integer
    p_min, p_max, n_points = p_range
    p = np.linspace(p_min, p_max, n_points)

    # Define which code-class to use 
    # code_class = BBcode.BBcode_Toric
    # code_class = BBcode.BBcode_ArXiV_example
    code_class = BBclass
    code_name = code_class.__name__
    decoder_name = decoder_dict['name']

    # Check if parity checks we implement are the same as in PanQEC 
    # test_code = code_class(4,4)
    # print(test_code.HX)
    # print(test_code.Hx.toarray())
    # print(np.all(test_code.HX == test_code.Hx.toarray()))
    # print(np.all(test_code.HZ == test_code.Hz.toarray()))

    # Must register the new code in panQEC 
    CODES[f'{code_name}'] = code_class
    DECODERS['BeliefPropagationLSDDecoder'] = BeliefPropagationLSDDecoder
    ERROR_MODELS['GaussianPauliErrorModel'] = GaussianPauliErrorModel

    save_frequency = 10  # Frequency of saving to file
    n_trials_str = f'{n_trials:.0e}'.replace('+0', '')
    grids_str = f'{grids}'.replace(' ', '').replace(':',';').replace(':', ';').replace("'", "").replace('_', '')
    p_range_str = f'{p_range}'.replace(' ', '')
    # Must copy and edit parameters from error_model_dict in a SAFE way. Trivial .copy() and edit does not work due to deepcopy
    parameters_copy = error_model_dict['parameters'][0].copy()
    rx = parameters_copy['r_x']
    ry = parameters_copy['r_y']
    rz = parameters_copy['r_z']
    parameters_copy['r_x'] = f"{rx:.2f}"
    parameters_copy['r_y'] = f"{ry:.2f}"
    parameters_copy['r_z'] = f"{rz:.2f}"
    parameters_str = f'[{parameters_copy}]'
    error_model_dict_str = f'{error_model_dict}'.replace(f"{error_model_dict['parameters']}", parameters_str).replace(' ', '').replace("'name':", '').replace(':', ';').replace("'", "").replace('_', '')
    decoder_dict_str = f'{decoder_dict}'.replace(' ', '').replace("'name':", '').replace(':', ';').replace("'", "")

    filename = f"data\{code_name};{grids_str};{error_model_dict_str};{decoder_dict_str}.json"

    # This magically fixes the fact that the filename is too long... 
    filename = u"\\\\?\\" + os.path.abspath(filename)

    rewrite_data = True
    if os.path.exists(filename):
        if ask_overwrite: advance = False
        else: advance = True
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
                'error_model': error_model_dict,
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


def plot_error_rates(analysis, 
                     savefig: bool=False, 
                     filename: str=None, 
                     include_threshold_estimate: bool=True):
    
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
            'font.size': 10,
            'figure.titlesize': 16} # extend as needed
    # print(plt.rcParams.keys())
    plt.rcParams.update(params)

    # Get colors from https://en.wikipedia.org/wiki/Pantone#Color_of_the_Year
    # Very Peri, Fuchsia Rose, Mimosa, Emerald, Classic Blue, Chili Pepper
    custom_cycle = ["#6667AB", "#C74375", "#F0C05A", "#009473", '#0F4C81', '#9B1B30']
    plt.rcParams["axes.prop_cycle"] = plt.cycler(color=custom_cycle)

    fig, ax = plt.subplots(ncols=3, figsize=(15, 5))

    plt.sca(ax[0])
    analysis.plot_thresholds(include_threshold_estimate=include_threshold_estimate)
    plt.sca(ax[1])
    analysis.plot_thresholds(sector='X', include_threshold_estimate=include_threshold_estimate)
    plt.sca(ax[2])
    analysis.plot_thresholds(sector='Z', include_threshold_estimate=include_threshold_estimate)

    fig.tight_layout()
    plt.show()

    # analysis.make_collapse_plots()

    fig, ax = plt.subplots(ncols=3, figsize=(15, 5))

    results = analysis.get_results()        # same as analysis.trunc_results['total']
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

    code_names =  results['code'][::n_trials_pr]
    code_params = results['code_params'][::n_trials_pr]
    error_models = results['error_model'][::n_trials_pr]
    decoders = results['decoder'][::n_trials_pr]
    biases = results['bias'][::n_trials_pr]
    num_logical_qubits = results['k'][::n_trials_pr].values

    capsize = 5
    ms = 5
    legend_lines = []
    legend_labels = []
    for i, Ls in enumerate(code_params):
        Lx, Ly, Lz = Ls.values()
        # code = eval('BBcode.' + results['code'][0] + f'({Lx}, {Ly})')
        index = results['code_params'] == Ls
        k = num_logical_qubits[i] 
        # lx = code.lx.toarray()
        # w_lx = np.sum(lx, axis=1)
        # k = np.max(w_lx)
        p_phys = results[index]['error_rate']
        kp = 1 - (1 - p_phys)**k        # 1 - (1-p)^k = k*p to first order

        line = ax[0].errorbar(p_phys, results[index]['p_est'], results[index]['p_se'],
                    label=rf'$L_x\!: {Lx}$, $L_y\!: {Ly}$, $k\!: {k}$', capsize=capsize, marker='o', ms=ms)
        linecolor = line[0].get_color()
        ax[0].plot(p_phys, kp, color=linecolor, linestyle=(0,(3,6)))
        ax[0].plot(p_phys, kp, color='k', linestyle=(4.5,(3,6)))

        ax[1].errorbar(p_phys, results[index]['p_est_X'], results[index]['p_se_X'],
                    label=rf'$L_x\!: {Lx}$, $L_y\!: {Ly}$', capsize=capsize, marker='o', ms=ms)
        ax[2].errorbar(p_phys, results[index]['p_est_Z'], results[index]['p_se_Z'],
                    label=rf'$L_x\!: {Lx}$, $L_y\!: {Ly}$', capsize=capsize, marker='o', ms=ms)
        
        legend_lines.append(line)
        legend_labels.append(rf'$L_x\!: {Lx}$, $L_y\!: {Ly}$, $k\!: {k}$')

    from matplotlib import lines

    th_line1 = lines.Line2D([], [], color='gray', linestyle=(0,(3,6)))
    th_line2 = lines.Line2D([], [], color='k', linestyle=(4.5,(3,6)))

    legend_lines.append((th_line1, th_line2))
    legend_labels.append('pseudo-threshold')

    result_X = analysis.sector_thresholds['X']
    result_Z = analysis.sector_thresholds['Z']
    if include_threshold_estimate:
        ax[0].axvline(analysis.thresholds.iloc[0]['p_th_fss'], color='red', linestyle='--')
        ax[0].axvspan(analysis.thresholds.iloc[0]['p_th_fss_left'], analysis.thresholds.iloc[0]['p_th_fss_right'],
                    alpha=0.5, color='pink')
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

    code_name = code_names[0]
    error_model = error_models[0]
    bias_label = str(biases[0]).replace('inf', '\\infty')
    decoder = decoders[0]
    fig.suptitle(f'{error_model} {code_name}\n'f'$\\eta_Z={bias_label}$, {decoder}\n')

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



def plot_compare_models(analysis1, analysis2, 
                        relevant_error_params: list=['r_x', 'r_y', 'r_z'], 
                        relevant_decoder_params: list=['max_bp_iter', 'lsd_order'], 
                        savefig=False):
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
            'font.size': 10,
            'figure.titlesize': 16} # extend as needed
    # print(plt.rcParams.keys())
    plt.rcParams.update(params)

    # Get colors from https://en.wikipedia.org/wiki/Pantone#Color_of_the_Year
    custom_cycle = ["#6667AB", "#C74375", "#F0C05A", "#009473", '#0F4C81', '#9B1B30']
    plt.rcParams["axes.prop_cycle"] = plt.cycler(color=custom_cycle)

    fig = plt.figure(figsize=(9, 5))
    gs = fig.add_gridspec(1, 2)
    ax = fig.add_subplot(gs[0, :])

    results1 = analysis1.get_results()        # same as analysis.trunc_results['total']
    # Count number of combinations of L_x and L_y
    dict_arr1 = np.array([[*L_dict.values()][:-1] for L_dict in results1['code_params']])
    n_Ls1 = 0 
    com_list = []
    for d in dict_arr1:
        lx, ly = d
        if (lx, ly) not in com_list:
            com_list.append((lx, ly))
            n_Ls1 += 1

    # Divide by total number of choices of L_x x L_y in parameters from input_data to only get each code param. one time 
    n_trials_pr1 = int(len(results1['n_trials'])/n_Ls1)

    code_names1 =  results1['code'][::n_trials_pr1]
    code_params1 = results1['code_params'][::n_trials_pr1]
    error_models1 = results1['error_model'][::n_trials_pr1]
    decoders1 = results1['decoder'][::n_trials_pr1]
    biases1 = results1['bias'][::n_trials_pr1]
    num_logical_qubits1 = results1['k'][::n_trials_pr1].values

    results2 = analysis2.get_results()        # same as analysis.trunc_results['total']
    # Count number of combinations of L_x and L_y
    dict_arr2 = np.array([[*L_dict.values()][:-1] for L_dict in results2['code_params']])
    n_Ls2 = 0 
    com_list = []
    for d in dict_arr2:
        lx, ly = d
        if (lx, ly) not in com_list:
            com_list.append((lx, ly))
            n_Ls2 += 1

    # Divide by total number of choices of L_x x L_y in parameters from input_data to only get each code param. one time 
    n_trials_pr2 = int(len(results2['n_trials'])/n_Ls2)

    code_names2 =  results2['code'][::n_trials_pr2]
    code_params2 = results2['code_params'][::n_trials_pr2]
    error_models2 = results2['error_model'][::n_trials_pr2]
    decoders2 = results2['decoder'][::n_trials_pr2]
    biases2 = results2['bias'][::n_trials_pr2]
    num_logical_qubits2 = results2['k'][::n_trials_pr2].values

    code_name1 = code_names1[0]
    error_model1 = error_models1[0]
    bias_label1 = str(biases1[0]).replace('inf', '\\infty')
    decoder1 = decoders1[0]
    code_name2 = code_names2[0]
    error_model2 = error_models2[0]
    bias_label2 = str(biases2[0]).replace('inf', '\\infty')
    decoder2 = decoders2[0]

    relevant_error_params = ['r_x', 'r_y', 'r_z']
    error_params_dict1 = results1['error_model_params'][0]
    relevant_error_param_dict1 = dict([(key, value if value%1==0 else f'{value:.2f}') for key, value in zip(error_params_dict1.keys(), error_params_dict1.values()) if key in relevant_error_params])
    if len(relevant_error_param_dict1) == 0: relevant_error_param_dict_str1 = ''
    else: relevant_error_param_dict_str1 = f"{relevant_error_param_dict1}"

    error_params_dict2 = results2['error_model_params'][0]
    relevant_error_param_dict2 = dict([(key, value if value%1==0 else f'{value:.2f}') for key, value in zip(error_params_dict2.keys(), error_params_dict2.values()) if key in relevant_error_params])
    if len(relevant_error_param_dict2) == 0: relevant_error_param_dict_str2 = ''
    else: relevant_error_param_dict_str2 = f"{relevant_error_param_dict2}"

    relevant_decoder_params = ['gaussian']
    decoder_params_dict1 = results1['decoder_params'][0]
    relevant_decoder_param_dict1 = dict([(key, value) for key, value in zip(decoder_params_dict1.keys(), decoder_params_dict1.values()) if key in relevant_decoder_params])
    if len(relevant_decoder_param_dict1) == 0: relevant_decoder_param_dict_str1 = ''
    else: relevant_decoder_param_dict_str1 = f"{relevant_decoder_param_dict1}"

    decoder_params_dict2 = results2['decoder_params'][0]
    relevant_decoder_param_dict2 = dict([(key, value) for key, value in zip(decoder_params_dict2.keys(), decoder_params_dict2.values()) if key in relevant_decoder_params])
    if len(relevant_decoder_param_dict2) == 0: relevant_decoder_param_dict_str2 = ''
    else: relevant_decoder_param_dict_str2 = f"{relevant_decoder_param_dict2}"

    error_params_str1 = f"{relevant_error_param_dict_str1}".replace("'", '').replace(':', ';').replace(' ', '')
    error_params_str2 = f"{relevant_error_param_dict_str2}".replace("'", '').replace(':', ';').replace(' ', '')
    decoder_params_str1 = f"{relevant_decoder_param_dict_str1}".replace("'", '').replace(':', ';').replace(' ', '')
    decoder_params_str2 = f"{relevant_decoder_param_dict_str2}".replace("'", '').replace(':', ';').replace(' ', '')
    filename = f"figures\compare_{code_name1}_{error_params_str1}_{decoder_params_str1}_{code_name2}_{error_params_str2}_{decoder_params_str2}.pdf"

    rewrite_plot = True
    if os.path.exists(filename) and savefig:
        advance = False
        while not advance:
            answer = input(f'Filename {filename} already exists. Do you want to write over the existing one (y/n)? ')
            if answer.lower() == 'y':
                advance = True
            elif answer.lower() == 'n':
                rewrite_plot = False
                advance = True 

    capsize = 5
    ms = 5
    from matplotlib import lines 
    analysis1_line = lines.Line2D([], [], color='gray', linestyle='solid')
    legend_lines1 = [analysis1_line]
    legend_labels1 = [f'Analysis 1']
    for i, Ls in enumerate(code_params1):
        Lx, Ly, Lz = Ls.values()
        # code = eval('BBcode.' + results1['code'][0] + f'({Lx}, {Ly})')
        index = results1['code_params'] == Ls
        k = num_logical_qubits1[i] 
        # lx = code.lx.toarray()
        # w_lx = np.sum(lx, axis=1)
        # k = np.max(w_lx)
        p_phys1 = results1[index]['error_rate']
        kp = 1 - (1 - p_phys1)**k        # 1 - (1-p)^k = k*p to first order

        line = ax.errorbar(p_phys1, results1[index]['p_est'], results1[index]['p_se'],
                    label=rf'$L_x\!: {Lx}$, $L_y\!: {Ly}$, $k\!: {k}$', capsize=capsize, marker='o', ms=ms)
        linecolor = line[0].get_color()
        ax.plot(p_phys1, kp, color=linecolor, linestyle=(0,(3,6)))
        ax.plot(p_phys1, kp, color='k', linestyle=(4.5,(3,6)))
        
        legend_lines1.append(line)
        legend_labels1.append(rf'$L_x\!: {Lx}$, $L_y\!: {Ly}$, $k\!: {k}$')
    # print(plt.rcParams.keys())

    # Reset matplotlib color cycle 
    plt.gca().set_prop_cycle(plt.cycler(color=custom_cycle))
    analysis2_line = lines.Line2D([], [], color='gray', linestyle=(0,(2,2)))
    legend_lines2 = [analysis2_line]
    legend_labels2 = [f'Analysis 2']
    for i, Ls in enumerate(code_params2):
        Lx, Ly, Lz = Ls.values()
        # code = eval('BBcode.' + results2['code'][0] + f'({Lx}, {Ly})')
        index = results2['code_params'] == Ls
        k = num_logical_qubits2[i] 
        # lx = code.lx.toarray()
        # w_lx = np.sum(lx, axis=1)
        # k = np.max(w_lx)
        p_phys2 = results2[index]['error_rate']
        kp = 1 - (1 - p_phys2)**k        # 1 - (1-p)^k = k*p to first order

        line = ax.errorbar(p_phys2, results2[index]['p_est'], results2[index]['p_se'],
                    label=rf'$L_x\!: {Lx}$, $L_y\!: {Ly}$, $k\!: {k}$', linestyle=(0,(2,2)), capsize=capsize, marker='o', ms=ms)
        linecolor = line[0].get_color()
        ax.plot(p_phys2, kp, color=linecolor, linestyle=(0,(3,6)))
        ax.plot(p_phys2, kp, color='k', linestyle=(4.5,(3,6)))
        
        legend_lines2.append(line)
        legend_labels2.append(rf'$L_x\!: {Lx}$, $L_y\!: {Ly}$, $k\!: {k}$')

    th_line1 = lines.Line2D([], [], color='gray', linestyle=(0,(3,6)))
    th_line2 = lines.Line2D([], [], color='k', linestyle=(4.5,(3,6)))

    legend_lines2.append((th_line1, th_line2))
    legend_labels2.append('pseudo-threshold')

    dist = max(p_phys1.max() - p_phys1.min(), p_phys2.max() - p_phys2.min())
    p_phys_min = min(p_phys1.min(), p_phys2.min())
    p_phys_max = max(p_phys1.max(), p_phys2.max())
    ax.set_xlim(p_phys_min-0.05*dist, p_phys_max+0.05*dist)
    ax.set_ylim(ymax=1.1)

    if 'gaussian' in results1['decoder_params'][::n_trials_pr1].values[0].keys():
        gaussian_decoding1 = results1['decoder_params'][::n_trials_pr1].values[0]['gaussian']
        gauss_decoder_str1 = f'(Gaussian={gaussian_decoding1}) '
    else: 
        gauss_decoder_str1 = ''
    if 'gaussian' in results2['decoder_params'][::n_trials_pr2].values[0].keys():
        gaussian_decoding2 = results2['decoder_params'][::n_trials_pr2].values[0]['gaussian']
        gauss_decoder_str2 = f'(Gaussian={gaussian_decoding2}) '
    else: 
        gauss_decoder_str2 = ''

    # ax.set_title(f'Analysis 1: {error_model1} {code_name1},{decoder1}' + gauss_decoder_str1 + f'$\\eta_Z={bias_label1}$', loc='left', wrap=True)
    # ax.set_title(f'Analysis 2: {error_model2} {code_name2},{decoder2}' + gauss_decoder_str2 + f'$\\eta_Z={bias_label2}$', loc='right', wrap=True)
    title = (
        r"$\begin{array}{l}"
        r"\text{Analysis 1: }" + r"\text{" + f"{error_model1} {code_name1}, $\\eta_Z={bias_label1}$" + r"}\\"
        r"\text{\phantom{Analysis 1: }}" + r"\text{" + f"{decoder1}" + gauss_decoder_str1 + r"}\\"
        r"\text{Analysis 2: }" + r"\text{" + f"{error_model2} {code_name2}, $\\eta_Z={bias_label2}$" + r"}\\"
        r"\text{\phantom{Analysis 2: }}" + r"\text{" + f"{decoder2}" + gauss_decoder_str2 + r"}\\"
        r"\end{array}$"
    )
    title1 = (
        r"$\begin{array}{l}"
        r"\text{Analysis 1: }" + r"\text{" + f"{code_name1}" + r"}\\"
        r"\text{" + f"{error_model1}, $\\eta_Z={bias_label1}$" + r"}\\"
        r"\text{" + f"{decoder1}" + gauss_decoder_str1 + r"}"
        r"\end{array}$"
    )
    title2 = (
        r"$\begin{array}{l}"
        r"\text{Analysis 2: }" + r"\text{" + f"{code_name2}" + r"}\\"
        r"\text{" + f"{error_model2}, $\\eta_Z={bias_label2}$" + r"}\\"
        r"\text{" + f"{decoder2}" + gauss_decoder_str2 + r"}"
        r"\end{array}$"
    )
    # ax.set_title(f'Analysis 1: {error_model1} {code_name1},\n {decoder1}' + gauss_decoder_str1 + f'$\\eta_Z={bias_label1}$\n' + f'Analysis 2: {error_model2} {code_name2},\n{decoder2}' + gauss_decoder_str2 + f'$\\eta_Z={bias_label2}$')
    ax.set_title(title1, loc='left')
    ax.set_title(title2, loc='right')

    ax.set_xlabel('Physical error rate')
    ax.set_ylabel('Logical error rate')
    ax.legend(legend_lines1 + legend_lines2, legend_labels1 + legend_labels2)

    fig.tight_layout()
    if savefig and rewrite_plot: plt.savefig(filename)
    plt.show()

"""
Make function calls for analysis and plotting 
"""


# analysis, filename = simulate_code(BBclass=BBcode.BBcode_A312_B312,
#                     decoder_dict={'name': 'BeliefPropagationLSDDecoder',
#                                   'parameters': [{'max_bp_iter': int(1e3), 'lsd_order': 10}]},
#                     error_model_dict={'name': 'GaussianPauliErrorModel', 
#                                       'parameters': [{'r_x': 1, 'r_y': 0, 'r_z': 0}]},
#                     n_trials=5e2, 
#                     grids=[{'L_x':6,'L_y':6},
#                            {'L_x': 12,'L_y':6},
#                            {'L_x': 18,'L_y':6},
#                            {'L_x': 24,'L_y':6}],
#                     p_range=(0, 0.25, 40))

# analysis, filename = simulate_code(BBclass=BBcode.BBcode_Ay3x1x2_Bx3y7y2,
#                     decoder_dict={'name': 'BeliefPropagationLSDDecoder',
#                                   'parameters': [{'max_bp_iter': int(1e3), 'lsd_order': 10}]},
#                     error_model_dict={'name': 'GaussianPauliErrorModel', 
#                                       'parameters': [{'r_x': 1, 'r_y': 0, 'r_z': 0}]},
#                     n_trials=5e2, 
#                     grids=[{'L_x':12,'L_y':12},
#                            {'L_x':18,'L_y': 18},
#                            {'L_x':24,'L_y':24},],
#                     p_range=(0, 0.25, 40))

# analysis, filename = simulate_code(BBclass=BBcode.BBcode_Toric,
#                     decoder_dict={'name': 'BeliefPropagationLSDDecoder',
#                                   'parameters': [{}]},
#                     error_model_dict={'name': 'GaussianPauliErrorModel', 
#                                       'parameters': [{'r_x': 1, 'r_y': 0, 'r_z': 0}]},
#                     n_trials=1e2, 
#                     grids=[{'L_x':10,'L_y':10},
#                            {'L_x': 20,'L_y':20},
#                            {'L_x': 30,'L_y':30}],
#                     p_range=(0, 0.25, 40))

# plot_error_rates(analysis, savefig=False, filename=filename.replace('data', 'figures').replace('.json', '.pdf'), include_threshold_estimate=False)


# analysis1, filename = simulate_code(BBclass=BBcode.BBcode_Toric,
#                     decoder_dict={'name': 'BeliefPropagationLSDDecoder',
#                                   'parameters': [{'gaussian': True}]},
#                     error_model_dict={'name': 'GaussianPauliErrorModel', 
#                                       'parameters': [{'r_x': 1, 'r_y': 0, 'r_z': 0}]},
#                     n_trials=1e2, 
#                     grids=[{'L_x':10,'L_y':10}, 
#                            {'L_x':20,'L_y':20}, 
#                            {'L_x':30,'L_y':30}],
#                     p_range=(0, 0.2, 40),
#                     ask_overwrite=False)

# analysis2, filename = simulate_code(BBclass=BBcode.BBcode_Toric,
#                     decoder_dict={'name': 'BeliefPropagationLSDDecoder',
#                                   'parameters': [{'gaussian': False}]},
#                     error_model_dict={'name': 'GaussianPauliErrorModel', 
#                                       'parameters': [{'r_x': 1, 'r_y': 0, 'r_z': 0}]},
#                     n_trials=1e2, 
#                     grids=[{'L_x':10,'L_y':10}, 
#                            {'L_x':20,'L_y':20}, 
#                            {'L_x':30,'L_y':30}],
#                     p_range=(0, 0.2, 40), 
#                     ask_overwrite=False)

analysis1, filename1 = simulate_code(BBclass=BBcode.BBcode_A312_B312,
                    decoder_dict={'name': 'BeliefPropagationLSDDecoder',
                                  'parameters': [{'max_bp_iter': int(1e3), 'lsd_order': 10, 'gaussian': True}]},
                    error_model_dict={'name': 'GaussianPauliErrorModel', 
                                      'parameters': [{'r_x': 1/3, 'r_y': 1/3, 'r_z': 1/3}]},
                    n_trials=1e2, 
                    grids=[{'L_x':6,'L_y':6},
                           {'L_x': 12,'L_y':6},
                           {'L_x': 18,'L_y':6},
                           {'L_x': 24,'L_y':6}],
                    p_range=(0, 0.3, 60),
                    ask_overwrite=True)

analysis2, filename2 = simulate_code(BBclass=BBcode.BBcode_A312_B312,
                    decoder_dict={'name': 'BeliefPropagationLSDDecoder',
                                  'parameters': [{'max_bp_iter': int(1e3), 'lsd_order': 10, 'gaussian': False}]},
                    error_model_dict={'name': 'GaussianPauliErrorModel', 
                                      'parameters': [{'r_x': 1/3, 'r_y': 1/3, 'r_z': 1/3}]},
                    n_trials=1e2, 
                    grids=[{'L_x':6,'L_y':6},
                           {'L_x': 12,'L_y':6},
                           {'L_x': 18,'L_y':6},
                           {'L_x': 24,'L_y':6}],
                    p_range=(0, 0.3, 60), 
                    ask_overwrite=True)

plot_compare_models(analysis1, analysis2, relevant_error_params=['r_x', 'r_y', 'r_z'], relevant_decoder_params=['gaussian'], savefig=True)
# plot_error_rates(analysis1, savefig=False, filename=filename1.replace('data', 'figures').replace('.json', '.pdf'), include_threshold_estimate=True)
