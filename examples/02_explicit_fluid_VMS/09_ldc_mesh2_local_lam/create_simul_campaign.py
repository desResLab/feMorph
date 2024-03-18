import os
import itertools
from create_fem_file_re100 import create_simulation_input

all_vals = [[1,2,3], # dt
            [1.0,10.0,100.0], # alpha
            [0,1], # use_euler
            [0,1], # use_B_hat
            [1], # include_K1
            [1], # include_K2
            [0,1], # include_K3
            [0]] # include_K4

all_combos = list(itertools.product(*all_vals))

# for loopA in range(len(all_combos)):
for loopA in range(2):

    params_dict = {}
    params_dict['dt']         = all_combos[loopA][0]
    params_dict['alpha']      = all_combos[loopA][1]
    params_dict['use_euler']  = all_combos[loopA][2]
    params_dict['use_B_hat']  = all_combos[loopA][3]
    params_dict['include_K1'] = all_combos[loopA][4]
    params_dict['include_K2'] = all_combos[loopA][5]
    params_dict['include_K3'] = all_combos[loopA][6]
    params_dict['include_K4'] = all_combos[loopA][7]

    # Create run folder name
    out_fld = 'run_'+str(loopA).zfill(3)

    # Create folder if it does not exist
    if not os.path.exists(out_fld):
        os.makedirs(out_fld)

    # Create run file
    create_simulation_input(params_dict,out_fld)