# Import libraries
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import colors 
import os
# Variables
geocmap = colors.ListedColormap(['gold','blue'])
geobounds = [1,1.5,2]
geonorm = colors.BoundaryNorm(geobounds, geocmap.N)

# Builders and table reading
def geo_builder(path,file_name):
    data = np.loadtxt(path + file_name,skiprows= 9)
    shape = (np.int16(np.max(data[:,2])),np.int16(np.max(data[:,0])))

    z = np.flip(np.reshape(data[:,5],newshape = shape,order='C'),axis=0)
    x = np.flip(np.reshape(data[:,3],newshape = shape,order='C'),axis=0)

    rel_perm = np.flip(np.reshape(data[:,6],newshape = shape,order='C'),axis=0)

    return rel_perm,x,z,shape 

def list_to_array(path,file_name):
    data = np.loadtxt(path + file_name,skiprows= 9)
    shape = (np.int16(np.max(data[:,2])),np.int16(np.max(data[:,0])))
    data = np.reshape(data[:,6],newshape = shape,order='C')
    data = np.flip(data,axis=0)

    return data,shape

def mul_list_to_array(path,file_name,shape,size=8):
    matrix = np.zeros((shape[0],shape[1],size))
    for i in range(size):
        skiprows = 9 + i*33009
        max_rows = 33000
        data = np.loadtxt(path + file_name,skiprows= skiprows,max_rows=max_rows)
        data = np.reshape(data[:,6],newshape = shape,order='C')
        data = np.flip(data,axis=0)
        matrix[:,:,i] = data

    return matrix

def special_table(path,file_name):
    special_df = pd.read_excel(path + file_name,sheet_name=0,skiprows= 5,names=['Time (days)','Date','Dissolved','Super-critical','Trapped'])

    return special_df

# Data processing
def indexes(diss,trapped,spcrit,total_mol):
    diss = diss.values
    trapped = trapped.values
    spcrit = spcrit.values

    total = total_mol

    rti = trapped/total #residual trapping index
    sti = diss/total #solubility trapping index
    tei = rti + sti #total trapping efficiency index

    return rti, sti, tei



# Build geological model using relative permeability set number for each simulation
def geo_properties(rel_perm,porosity,perm):

    titles= ['Geologic model','Porosity','Permeability']

    fig,ax = plt.subplots(1,3,figsize=(10,2),sharey=True,dpi =300)
    im0 = ax[0].imshow(rel_perm,cmap=geocmap,norm=geonorm,origin='lower',interpolation='nearest',aspect='auto')
    im1 = ax[1].imshow(porosity,cmap='jet',aspect='auto')
    im2 = ax[2].imshow(perm,cmap='jet',aspect='auto',vmin=0,vmax=1000)
    ax[0].spines['left'].set_position(('axes', -0.05))
    ax[0].set_yticks([0,219])

    for i in range(3):
        ax[i].set_title(titles[i])
        ax[i].spines['bottom'].set_position(('axes', -0.05))
        ax[i].set_xticks([0,149])
        ax[i].set_xticklabels(['0','32850 ft'])
        ax[i].set_yticklabels(['8000','8680 ft'])

    cbar = plt.colorbar(im0,cmap=geocmap,norm=geonorm,boundaries=geobounds,ticks=[1,2],ax=ax[0],fraction=0.046, pad=0.04,shrink=0.5)
    cbar.set_ticklabels(['Sand','Shale'])

    plt.colorbar(im1,ax=ax[1],fraction=0.046, pad=0.04,shrink=0.5)
    cbar1 = plt.colorbar(im2,ax=ax[2],fraction=0.046, pad=0.04,shrink=0.5)
    cbar1.set_label('k [mD]')

    plt.show()

# Pressure buildup plot
def pres_bup(p_init,pbup,z):
    
    mean_pinit = np.mean(p_init,axis=1)
    mean_pbup = np.mean(pbup,axis=1)

    fig, ax = plt.subplots(figsize=(5,5))
    ax.plot(mean_pbup,z[:,0],color='blue',linewidth=1,linestyle='--',label='Pressure')
    ax.plot(mean_pinit,z[:,0],color='skyblue',linewidth=1,label='Initial Pressure')
    ax.set_ylim(8000,8680)
    ax.set_xlabel('Pressure [psi]')
    ax.set_ylabel('Depth [ft]')
    ax.set_title('Pressure buildup')
    ax.invert_yaxis()
    ax.legend()
    
    plt.show()

# Saturation field along depth (free and trapped gas)
def saturation_profile(gas_sat,z,scrit,gas_trapped = {},trapped = False):
    sat_copy = gas_sat.copy()
    sat_copy[sat_copy<0.0003] = np.nan

    mean_sat = np.nanmean(sat_copy,axis=1)
    if trapped == True:
        trp_copy = gas_trapped.copy()
        trp_copy[trp_copy<0.0003] = np.nan

        mean_trp = np.nanmean(trp_copy,axis=1)

    fig, ax = plt.subplots(figsize=(5,5))
    s1 = ax.plot(mean_sat[:,1],z[:,0],color='blue',linewidth=1,linestyle='--',label='2030')
    ax.plot(mean_sat[:,4],z[:,0],color='green',linewidth=1,linestyle='--',label='2130')
    ax.plot(mean_sat[:,7],z[:,0],color='red',linewidth=1,linestyle='--',label='2330')
    ax.fill_betweenx(z[:,0],mean_sat[:,1],scrit,where=mean_sat[:,1]>scrit,color='blue',alpha=0.5)
    ax.fill_betweenx(z[:,0],mean_sat[:,4],scrit,where=mean_sat[:,4]>scrit,color='green',alpha=0.5)
    ax.fill_betweenx(z[:,0],mean_sat[:,7],scrit,where=mean_sat[:,7]>scrit,color='red',alpha=0.5,label='Free gas')
    if trapped == True:
        ax.plot(mean_trp[:,7],z[:,0],color='orange',linewidth=1)
        ax.fill_betweenx(z[:,0],mean_trp[:,7],scrit,where=mean_trp[:,7]>scrit,color='orange',label='Residual trapped gas (2330)')
    v = ax.vlines(scrit,8000,8680,linestyle='--',linewidth=1,color='black')
   
    ax.text(scrit+0.01,8120,'Critical saturation',rotation=90,fontsize=5)
    ax.set_ylim(8000,8680)
    ax.set_xlabel('Gas saturation')
    ax.set_ylabel('Depth [ft]')
    ax.invert_yaxis()

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=3, fancybox=True, shadow=True)
    plt.show()

    return fig

# Plot trapping indexes
def plot_indexes(date,rti,sti):
    fig, ax = plt.subplots(figsize=(5,5))
    ax.bar(date,width=100,height=sti,label = 'STI')
    ax.bar(date,width=100,height=rti,label = 'RTI',bottom=sti)
    ax.set_xlabel('Year')
    ax.set_ylabel('Total Trapping index = RTI + STI')
    ax.set_ylim(0,1)

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2, fancybox=True, shadow=True)

    plt.show()

    return fig

def time_series_plot(diss,trapped,spcrit,dates):
    diss = diss.values * (44*0.0005) #tons
    trapped = trapped.values * (44*0.0005) #tons
    spcrit = spcrit.values * (44*0.0005) #tons
    fig, ax = plt.subplots(figsize=(10,5))
    ax.plot(dates,trapped,linewidth=3,color='black',linestyle='dashed',label='CO2 Trapped (tons)')
    ax.plot(dates,spcrit,linewidth=3,color='red',label='Supercritical CO2 (tons)')
    ax.set_title('CO2 Trapped vs Supercritical CO2')
    ax.set_xlabel('Date')
    ax.set_ylabel('CO2 (tons)')
    ax.legend()
    plt.show()
    return 

# Plot plume centroid over time
def plot_plume_centroid(cx,cz,dates):
    fig,ax = plt.subplots(1,2,figsize=(10,5))
    ax[0].scatter(dates,cx,c='red',marker='o')
    ax[0].set_ylim(0,32850)
    ax[0].set_xlabel('Year')
    ax[0].set_title('Plume centroid x [ft] over time')

    ax[1].scatter(dates,cz,c='red',marker='o')
    ax[1].set_ylim(8000,8680)
    ax[1].invert_yaxis()
    ax[1].set_xlabel('Year')
    ax[1].set_title('Plume centroid z [ft] over time')

    plt.show()

# Batch plots pressure, moments and plume
def plot_batch_pressure(p_dict,z,size):
    colors = ['blue','green','red','orange','purple','brown','pink']
    mean_pinit = np.mean(p_dict['S1'][:,:,0],axis=1)
    fig, ax = plt.subplots(figsize=(5,5))
    for i in range(size):
        mean_pbup = np.mean(p_dict[f'S{i+1}'][:,:,1],axis=1)
        ax.plot(mean_pbup,z[:,0],color=colors[i],linewidth=1,linestyle='--',label=f'S{i+1}')
    ax.plot(mean_pinit,z[:,0],color='skyblue',linewidth=1,label='Pi')
    ax.set_ylim(8000,8680)
    ax.set_xlabel('Pressure [psi]')
    ax.set_ylabel('Depth [ft]')
    ax.set_title('Pressure buildup')
    ax.invert_yaxis()
    ax.legend()

    plt.show()

def batch_plot_moments(sat,poro,vol,x,z,dates,scrit,cl):
    por_vol = poro*vol
    c10 = por_vol*x
    c01 = por_vol*z
    c20 = por_vol*x**2
    c02 = por_vol*z**2

    mts_dict = {}
    fig, ax = plt.subplots(1,3,figsize=(15,5))
    for i in range(len(sat)):
        
        m00 = np.sum(sat[f'S{i+1}']*por_vol[:,:,None],axis=(0,1))
        m10, m01 = np.sum(sat[f'S{i+1}']*c10[:,:,None],axis=(0,1)), np.sum(sat[f'S{i+1}']*c01[:,:,None],axis=(0,1))
        m20, m02 = np.sum(sat[f'S{i+1}']*c20[:,:,None],axis=(0,1)), np.sum(sat[f'S{i+1}']*c02[:,:,None],axis=(0,1))

        #center of mass
        cx = m10/m00
        cz = m01/m00

        #spreading  
        sx = m20/m00 - cx**2
        sz = m02/m00 - cz**2

        ax[0].scatter(dates,m00,c=cl[i],s=5,label=f'S{i+1}-{scrit[i]}')
        ax[1].scatter(dates,cx,c=cl[i],s=5,label=f'S{i+1}-{scrit[i]}')
        ax[2].scatter(dates,sx,c=cl[i],s=5,label=f'S{i+1}-{scrit[i]}')
        mts_dict[f'S{i+1}'] = [m00,cx,cz,sx,sz]

    ax[0].vlines(2030,0,1e7,linestyle='--',linewidth=0.5,color='black')
    ax[0].text(2033,0,'End of injection',rotation=90,fontsize=7)
    ax[0].set_title('CO2 volume over time')
    ax[1].set_title('CO2 centroid (lateral extent) over time')
    ax[2].set_title('CO2 lateral spreading over time')

    for j in range(3):
        ax[j].set_xlabel('Year')

    ax[0].set_ylabel('Volume [ft3]')
    ax[1].set_ylabel('Lateral extent [ft]')
    ax[2].set_ylabel('Lateral spreading [ft2]')

    handles, labels = ax[0].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys(), loc='lower center',ncol=7,
            fontsize='large',bbox_to_anchor = (0.5,-0.07),fancybox=True, shadow=True)
    
    plt.show()

    return mts_dict

def plot_plume(rel_perm,gas_dict,size):
    fig,ax = plt.subplots(3,size,figsize=(20,5),dpi=300)
    for i in range(size):
        ax[0,i].imshow(rel_perm,cmap='binary',aspect='auto')
        im = ax[0,i].imshow(gas_dict[f'S{i+1}'][:,:,1],cmap='jet',vmin=0.0,vmax=1,aspect='auto',alpha=0.9)
        ax[0,i].set_title(f'Scrit: 0.{i+1}')
        ax[0,i].set_xticks([])
        ax[0,i].set_yticks([])

        ax[1,i].imshow(rel_perm,cmap='binary',aspect='auto')
        ax[1,i].imshow(gas_dict[f'S{i+1}'][:,:,4],cmap='jet',vmin=0.0,vmax=1,aspect='auto',alpha=0.9)
        ax[1,i].set_xticks([])
        ax[1,i].set_yticks([])

        ax[2,i].imshow(rel_perm,cmap='binary',aspect='auto')
        ax[2,i].imshow(gas_dict[f'S{i+1}'][:,:,6],cmap='jet',vmin=0.0,vmax=1,aspect='auto',alpha=0.9)
        ax[2,i].set_xticks([])
        ax[2,i].set_yticks([])

        ax[0,0].set_ylabel('2030')
        ax[1,0].set_ylabel('2130')
        ax[2,0].set_ylabel('2330')

    fig.colorbar(im, ax=ax.ravel().tolist(), orientation = 'vertical',shrink=0.5,pad=0.04,ticks=[0,0.2,0.4,0.6,0.8,1])

    return plt.show()
                                                      
def plot_contribution_residual_sat(scrit,trapped_dict,size,scenarios):
    mean_trp = []
    gas_sat = np.arange(0.1,1.1,0.1)
    for i in range(size):
        trp = trapped_dict[f'S{i+1}'][:,:,-1]
        trp_copy = trp.copy()
        trp_copy[trp_copy<=scrit[i]] = np.nan
        mean_trp.append(np.nanmean(trp_copy))

    mean_trp = np.array(mean_trp)
    cont = mean_trp - scrit

    fig, ax = plt.subplots(figsize=(5,5))
    ax.bar(scrit,height=cont,width=0.05,color='blue',label='Imbibition hysteresis ',bottom=scrit)
    ax.bar(scrit,height=scrit,width=0.05,color='red',label='Capillary Heterogeneity Trapping')
    ax.set_xlabel('Scenarios')
    ax.set_ylabel('Contribution to total trapping')
    ax.set_xticks(scrit)
    ax.set_xticklabels(scenarios)
    ax.set_ylim(0,1)
    ax.set_yticks(gas_sat)

    ax.legend(fontsize='small')
    fig.tight_layout()
    plt.show()

    return mean_trp

def count_grid_cells(gas_dict):
   
    xcount = {}
    zcount = {}
    satcount = {}
    for key in gas_dict.keys():
        gasc = gas_dict[key].copy()
        gasc[gasc<0.003] = 0
        xcount[key] = np.max(np.count_nonzero(gasc[:,:,-1],axis=1,keepdims=True))
        zcount[key] = np.max(np.count_nonzero(gasc[:,:,-1],axis=0,keepdims=True))
        satcount[key] = np.count_nonzero(gasc[:,:,-1])

    return xcount, zcount, satcount

        


