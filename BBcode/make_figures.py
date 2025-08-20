import src.BBcode_classes as BBcode
import src.TileCode_classes as TileCode
from panqec.codes import Toric2DCode

from src.analysis import simulate_code, plot_error_rates, plot_compare_models, plot_all

"""
* Code written by Anton Brekke * 

This file calculates the error threshold given a BBCode class from 'BBcode_classes.py'. 
"""

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
#                     n_trials=1e3, 
#                     grids=[{'L_x':10,'L_y':10}, 
#                            {'L_x':20,'L_y':20}, 
#                            {'L_x':30,'L_y':30}],
#                     p_range=(0, 0.2, 40), 
#                     ask_overwrite=False)

# analysis2, input_data2, filename2 = simulate_code(
#                     BBclass=BBcode.BBcode_A312_B312,
#                     decoder_dict={'name': 'BeliefPropagationLSDDecoder',
#                                   'parameters': [{'max_bp_iter': int(1e3), 'lsd_order': 10, 'gaussian': False}]},
#                     error_model_dict={'name': 'GaussianPauliErrorModel', 
#                                       'parameters': [{'r_x': 1/3, 'r_y': 1/3, 'r_z': 1/3}]},
#                     n_trials=1e3, 
#                     grids=[{'L_x': 6,'L_y':6},
#                            {'L_x': 12,'L_y':6},
#                            {'L_x': 18,'L_y':6},
#                            {'L_x': 24,'L_y':6}],
#                     p_range=(0, 0.3, 60), 
#                     ask_overwrite=True)

"""
Compare BBcode with Toric code
Same error model, same decoder, same number of physical qubits, same number of logical qubits
"""

# k = 12 logical qubits
# analysis1, input_data1, filename1 = simulate_code(
#                     BBclass=BBcode.BBcode_A312_B312,
#                     decoder_dict={'name': 'BeliefPropagationLSDDecoder',
#                                   'parameters': [{'max_bp_iter': int(1e3), 'lsd_order': 10, 'gaussian': True}]},
#                     error_model_dict={'name': 'GaussianPauliErrorModel', 
#                                       'parameters': [{'r_x': 1, 'r_y': 0, 'r_z': 0}]},
#                     n_trials=1e2, 
#                     grids=[{'L_x':18,'L_y':12}],
#                     p_range=(0, 0.3, 60),
#                     ask_overwrite=True)

# # Since k = 12, simulate 6 Toric codes (since they have k=2 each) and collapse results later 
# analysis2, input_data2, filename2 = simulate_code(
#                     BBclass=BBcode.BBcode_Toric,
#                     decoder_dict={'name': 'BeliefPropagationLSDDecoder',
#                                   'parameters': [{'max_bp_iter': int(1e3), 'lsd_order': 10, 'gaussian': True}]},
#                     error_model_dict={'name': 'GaussianPauliErrorModel', 
#                                       'parameters': [{'r_x': 1, 'r_y': 0, 'r_z': 0}]},
#                     n_trials=1e2, 
#                     grids=6*[{'L_x':6,'L_y':6}],
#                     p_range=(0, 0.3, 60), 
#                     ask_overwrite=True)

# k = 16 logical qubits 
# analysis1, input_data1, filename1 = simulate_code(
#                     BBclass=BBcode.BBcode_A312_B312,
#                     decoder_dict={'name': 'BeliefPropagationLSDDecoder',
#                                   'parameters': [{'max_bp_iter': int(1e3), 'lsd_order': 10, 'gaussian': True}]},
#                     error_model_dict={'name': 'GaussianPauliErrorModel', 
#                                       'parameters': [{'r_x': 1, 'r_y': 0, 'r_z': 0}]},
#                     n_trials=1e2, 
#                     grids=[{'L_x':24,'L_y':24}],
#                     p_range=(0, 0.3, 60),
#                     ask_overwrite=True)

# analysis2, input_data2, filename2 = simulate_code(
#                     BBclass=BBcode.BBcode_Toric,
#                     decoder_dict={'name': 'BeliefPropagationLSDDecoder',
#                                   'parameters': [{'max_bp_iter': int(1e3), 'lsd_order': 10, 'gaussian': True}]},
#                     error_model_dict={'name': 'GaussianPauliErrorModel', 
#                                       'parameters': [{'r_x': 1, 'r_y': 0, 'r_z': 0}]},
#                     n_trials=1e2, 
#                     grids=8*[{'L_x':9,'L_y':8}],
#                     p_range=(0, 0.3, 60), 
#                     ask_overwrite=True)

# analysis1, input_data1, filename1 = simulate_code(
#                     BBclass=BBcode.BBcode_A312_B312,
#                     decoder_dict={'name': 'BeliefPropagationLSDDecoder',
#                                   'parameters': [{'max_bp_iter': int(1e3), 'lsd_order': 10, 'gaussian': True}]},
#                     error_model_dict={'name': 'GaussianPauliErrorModel', 
#                                       'parameters': [{'r_x': 1, 'r_y': 0, 'r_z': 0}]},
#                     n_trials=1e2, 
#                     grids=[{'L_x':24,'L_y':12}],
#                     p_range=(0, 0.3, 60),
#                     ask_overwrite=True)

# analysis2, input_data2, filename2 = simulate_code(
#                     BBclass=BBcode.BBcode_Toric,
#                     decoder_dict={'name': 'BeliefPropagationLSDDecoder',
#                                   'parameters': [{'max_bp_iter': int(1e3), 'lsd_order': 10, 'gaussian': True}]},
#                     error_model_dict={'name': 'GaussianPauliErrorModel', 
#                                       'parameters': [{'r_x': 1, 'r_y': 0, 'r_z': 0}]},
#                     n_trials=1e2, 
#                     grids=8*[{'L_x':6,'L_y':6}],
#                     p_range=(0, 0.3, 60), 
#                     ask_overwrite=True)

# k = 8 logical qubits 
# analysis1, input_data1, filename1 = simulate_code(
#                     BBclass=BBcode.BBcode_A312_B312,
#                     decoder_dict={'name': 'BeliefPropagationLSDDecoder',
#                                   'parameters': [{'max_bp_iter': int(1e3), 'lsd_order': 10, 'gaussian': True}]},
#                     error_model_dict={'name': 'GaussianPauliErrorModel', 
#                                       'parameters': [{'r_x': 1, 'r_y': 0, 'r_z': 0}]},
#                     n_trials=1e2, 
#                     grids=[{'L_x':15,'L_y':15}],
#                     p_range=(0, 0.3, 60),
#                     ask_overwrite=True)

# analysis2, input_data2, filename2 = simulate_code(
#                     BBclass=BBcode.BBcode_Toric,
#                     decoder_dict={'name': 'BeliefPropagationLSDDecoder',
#                                   'parameters': [{'max_bp_iter': int(1e3), 'lsd_order': 10, 'gaussian': True}]},
#                     error_model_dict={'name': 'GaussianPauliErrorModel', 
#                                       'parameters': [{'r_x': 1, 'r_y': 0, 'r_z': 0}]},
#                     n_trials=1e2, 
#                     grids=4*[{'L_x':8,'L_y':7}],
#                     p_range=(0, 0.3, 60), 
#                     ask_overwrite=True)


"""Compare BBcode with Gaussian decoding vs. without Gaussian decoding"""

# analysis1, input_data1, filename1 = simulate_code(
#                     BBclass=BBcode.BBcode_A312_B312,
#                     decoder_dict={'name': 'BeliefPropagationLSDDecoder',
#                                   'parameters': [{'max_bp_iter': int(1e3), 'lsd_order': 10, 'gaussian': True}]},
#                     error_model_dict={'name': 'GaussianPauliErrorModel', 
#                                       'parameters': [{'r_x': 1, 'r_y': 0, 'r_z': 0}]},
#                     n_trials=1e3, 
#                     grids=[{'L_x': 6,'L_y':6},
#                            {'L_x': 12,'L_y':6},
#                            {'L_x': 18,'L_y':12},
#                            {'L_x': 18,'L_y':18}],
#                     p_range=(0, 0.3, 60), 
#                     ask_overwrite=True)

# analysis2, input_data2, filename2 = simulate_code(
#                     BBclass=BBcode.BBcode_A312_B312,
#                     decoder_dict={'name': 'BeliefPropagationLSDDecoder',
#                                   'parameters': [{'max_bp_iter': int(1e3), 'lsd_order': 10, 'gaussian': False}]},
#                     error_model_dict={'name': 'GaussianPauliErrorModel', 
#                                       'parameters': [{'r_x': 1, 'r_y': 0, 'r_z': 0}]},
#                     n_trials=1e3, 
#                     grids=[{'L_x': 6,'L_y':6},
#                            {'L_x': 12,'L_y':6},
#                            {'L_x': 18,'L_y':12},
#                            {'L_x': 18,'L_y':18}],
#                     p_range=(0, 0.3, 60), 
#                     ask_overwrite=True)

"""Compare BBcode with Tile code"""

analysis1, input_data1, filename1 = simulate_code(
                    BBclass=BBcode.BBcode_A312_B312,
                    decoder_dict={'name': 'BeliefPropagationLSDDecoder',
                                  'parameters': [{'max_bp_iter': int(1e3), 'lsd_order': 10, 'gaussian': True}]},
                    error_model_dict={'name': 'GaussianPauliErrorModel', 
                                      'parameters': [{'r_x': 1, 'r_y': 0, 'r_z': 0}]},
                    n_trials=5e2, 
                    grids=[{'L_x':15,'L_y':12}],
                    p_range=(0, 0.3, 60), 
                    ask_overwrite=True)

analysis2, input_data2, filename2 = simulate_code(
                    BBclass=TileCode.TileCode_B3_W6,
                    decoder_dict={'name': 'BeliefPropagationLSDDecoder',
                                  'parameters': [{'max_bp_iter': int(1e3), 'lsd_order': 10, 'gaussian': True}]},
                    error_model_dict={'name': 'GaussianPauliErrorModel', 
                                      'parameters': [{'r_x': 1, 'r_y': 0, 'r_z': 0}]},
                    n_trials=5e2, 
                    grids=[{'L_x':15,'L_y':12}],
                    p_range=(0, 0.3, 60), 
                    ask_overwrite=True)

analysis3, input_data3, filename3 = simulate_code(
                    BBclass=Toric2DCode,
                    decoder_dict={'name': 'BeliefPropagationLSDDecoder',
                                  'parameters': [{'max_bp_iter': int(1e3), 'lsd_order': 10, 'gaussian': True}]},
                    error_model_dict={'name': 'GaussianPauliErrorModel', 
                                      'parameters': [{'r_x': 1, 'r_y': 0, 'r_z': 0}]},
                    n_trials=5e2, 
                    grids=4*[{'L_x':9,'L_y':5}],
                    p_range=(0, 0.3, 60), 
                    ask_overwrite=True)

# plot_error_rates(analysis1, savefig=False, filename=filename1.replace('data', 'figures').replace('.json', '.pdf'), include_threshold_estimate=True)
# plot_compare_models((analysis1, input_data1), (analysis2, input_data2), relevant_error_params=['r_x', 'r_y', 'r_z'], relevant_decoder_params=['gaussian'], savefig=True, collapse_1=False, collapse_2=False)
plot_all([analysis1, analysis2, analysis3], [input_data1, input_data2, input_data3])
