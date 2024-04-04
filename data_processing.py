import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from scipy import stats

# Insert files path
path = 'C:/Users/ubillusj/Desktop/CMG/Ubillus_simulations/Test_scenario/'

# Read GSlib file into a pandas dataframe
poro = pd.read_csv(path + 'test Porosity 2010-Jan-01.gslib',sep = ' ',skiprows= 9,index_col=False, header = None,names= ['i_index','j_index','k_index','x','y','z','porosity'])
volume = pd.read_csv(path + 'test Gross Block Volume 2010-Jan-01.gslib',sep = ' ',skiprows= 9,index_col=False, header = None,names= ['i_index','j_index','k_index','x','y','z','volume'])
pres_bup = pd.read_csv(path + 'test Pres_BUp 2030-Dec-31.gslib',sep = ' ',skiprows= 9,index_col=False, header = None,names= ['i_index','j_index','k_index','x','y','z','pres_bup'])
gas_saturation = pd.read_csv(path + 'test Gas Saturation.gslib',sep = ' ',skiprows= 9,nrows=33000,index_col=False, header = None,names= ['i_index','j_index','k_index','x','y','z','gas_saturation'])

# Estimate average gas saturation of the domain from the pandas dataframe
mean_copy = gas_saturation.copy()
mean_copy['gas_saturation'][mean_copy['gas_saturation']<0.003] = np.nan
gas_saturation_mean = np.nanmean(mean_copy['gas_saturation'])

# Create a grid using the i,j,k indexes and the gas saturation values
sat_array = np.reshape(gas_saturation['gas_saturation'].values,newshape = (np.max(gas_saturation['k_index']),np.max(gas_saturation['i_index'])),order='C')
sat_array = np.flip(sat_array,axis=0)
poro_array = np.reshape(poro['porosity'].values,newshape = (np.max(gas_saturation['k_index']),np.max(gas_saturation['i_index'])),order='C')
poro_array = np.flip(poro_array,axis=0)
vol_array = np.reshape(volume['volume'].values,newshape = (np.max(gas_saturation['k_index']),np.max(gas_saturation['i_index'])),order='C')
vol_array = np.flip(vol_array,axis=0)

# Create position arrays for the grid
position_x = np.reshape(poro['x'].values,newshape = (np.max(gas_saturation['k_index']),np.max(gas_saturation['i_index'])),order='C')
position_x = np.flip(position_x,axis=0)
position_z = np.reshape(poro['z'].values,newshape = (np.max(gas_saturation['k_index']),np.max(gas_saturation['i_index'])),order='C')
position_z = np.flip(position_z,axis=0)

# Moments
def moments(sat_array,poro_array,vol_array,position_x,position_z):
    sat_array[sat_array<0.003] = 0
    m00 = np.round(np.sum(vol_array*sat_array*poro_array),2)
    m10 = np.sum(position_x*vol_array*sat_array*poro_array)
    m01 = np.sum(position_z*vol_array*sat_array*poro_array)
    m20 = np.sum(position_x**2*vol_array*sat_array*poro_array)
    m02 = np.sum(position_z**2*vol_array*sat_array*poro_array)

    #centroid of mass
    cx = np.round(m10/m00,2)
    cz = np.round(m01/m00,2)
    #lateral spreading (mm2)
    lx = np.round(m20/m00 - cx**2,2)
    lz = np.round(m02/m00 - cz**2,2)

    return m00,cx,cz,lx,lz









    




