import numpy as np
import matplotlib.pyplot as plt
def estimate_avg_gas(gas_dict):
    ret_dict = {}
    for key in gas_dict.keys():
        temp = gas_dict[key].copy()
        temp[temp<0.0001] = np.nan
        ret_dict[key] = np.nanmean(temp,axis=(0,1))

    return ret_dict
    

# Compare gas saturations in different set of simulations
def compare_avg_gas(sim_dict):
    '''
    
    '''
    r_dict = {}

    # Estimate average gas saturation for each simulation
    for key in sim_dict.keys():
        r_dict[key] = estimate_avg_gas(sim_dict[key])

    # Plot average gas saturation for each simulation
    fig, ax = plt.subplots(1,1,figsize=(10,10),dpi=300)
    for key in sim_dict.keys():
        for key2 in r_dict.keys():
            ax.plot(r_dict[key2][key])
    plt.show()

        

        

    
