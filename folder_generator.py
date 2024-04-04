import os
path = 'C:/Users/ubillusj/Desktop/CMG/Ubillus_simulations/hyskrg_midpoint/'
n_scenarios = 4
# Create results folder
if not os.path.exists(path + 'results'):
    os.makedirs(path + 'results_realistic')

    # create scenarions folder inside results folder from 1 to n_scenarios
    for i in range(1, n_scenarios + 1):
        os.makedirs(path + 'results_realistic/S' + str(i))
        '''
        # create subfolders inside each scenario folder
        for j in range(1, 4):
            os.makedirs(path + 'results/S' + str(i) + '/S' + str(i) + str(j))
        '''
    

