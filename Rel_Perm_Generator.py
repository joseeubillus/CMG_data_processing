# Batch generator
import numpy as np
import matplotlib.pyplot as plt

# Build relative permeability curves using Modified Brooks-Corey model 

def rel_perm_batch(path_to_save,init_crit_co2,final_crit_co2,
                   swi=0.2,pce_1=1.93,pce_2=2.14e3,kx_s=500,kx_sh=0.01,cap_p_model='MBC',print_bool=False):
    '''
        This function generates relative permeability curves and capillary pressure curves for a range of critical
        CO2 saturations (Sgc).

        Parameters:
        -----------
        path_to_save: str
            Path to save the generated files
        init_crit_co2: float
            Initial critical CO2 saturation
        final_crit_co2: float
            Final critical CO2 saturation
        cap_p_model: str
            Capillary pressure model to use. Options: 'MBC' (Modified Brooks Corey), BC (Brooks Corey) and VG (Van Genuchten)
        swi: float
            Irreducible water saturation
        pce_1: float
            Entry pressure for sand (psi)
        pce_2: float
            Entry pressure for shale (psi)
        kx_s: float
            Permeability for sand (mD)
        kx_sh: float
            Permeability for shale (mD)
        print: bool
            If True, the generated files are saved in the path_to_save folder

        Returns:
        --------
        fig: matplotlib figure
            Figure with the relative permeability and capillary pressure curves
    
    '''
    colors = ['b','r','g','c','m','y','k']
    labels = ['Sc = 0.1','Sc = 0.2','Sc = 0.3','Sc = 0.4','Sc = 0.5','Sc = 0.6','Sc = 0.7']

    crit_c02_list = np.arange(init_crit_co2,final_crit_co2 + 0.1,0.1)
    wat_sat = np.arange(swi,1.01,0.01)
    sg= 1-wat_sat[::-1]
        
    fig, ((ax0,ax1),(ax2,ax3)) = plt.subplots(2,2,figsize=(20,20),dpi=300)

    for snc in crit_c02_list:

        #Relative permeability curves
        S_ws = (wat_sat - swi)/(1-swi)
        S_wh = (wat_sat - swi)/(1-swi-snc)

        krl = (S_ws**4)
        krg = 0.4*(1-S_wh**2)*(1-S_wh)**2 
        krg[krg<0]=0
        wat_sat[0] = wat_sat[0] + 0.001

        if cap_p_model == 'MBC':
            idx = np.where(np.isclose(wat_sat,(1 - snc)))[0][0]
            SwP = wat_sat[:idx+1]
            SwPb = wat_sat[idx+1:]
        # Modified Brooks Corey model

            Pc_facies1a = pce_1 * ((SwP - swi)/(1-swi))**(-0.5) # Capillary pressure (psi)
            Pc_facies1b = ((pce_1/snc)*((1-snc-swi)/(1-swi))**(-0.5))*(1-SwPb) # Capillar pressure for the saturations lower than Snc
            pc1 = np.concatenate((Pc_facies1a,Pc_facies1b))

            Pc_facies2a = pce_2 * ((SwP - swi)/(1-swi))**(-0.5) # Capillary pressure (psi)
            Pc_facies2b = ((pce_2/snc)*((1-snc-swi)/(1-swi))**(-0.5))*(1-SwPb) # Capillar pressure for the saturations lower than Snc

            pc2 = np.concatenate((Pc_facies2a,Pc_facies2b))
        elif cap_p_model == 'BC':
            pc1 = pce_1 *((wat_sat - swi)/(1-swi))**(-0.5) # Sand Capillary pressure (psi)
            pc2 = pce_2 *((wat_sat - swi)/(1-swi))**(-0.5) # Shale Capillary pressure (psi)

        elif cap_p_model == 'VG':
            pc1 = pce_1*(((wat_sat - swi)/(1-swi))**(-1/0.464)-1)**(1/1.864) # Sand Capillary pressure (psi)
            pc2 = pce_2*(((wat_sat - swi)/(1-swi))**(-1/0.559)-1)**(1/2.266) # Shale Capillary pressure (psi)

        # J-functions for sand and shale
        Phi_s = (kx_s/7e7)**(1/9.6)

        J_s = pc1*np.sqrt(kx_s/Phi_s)

        Phi_sh = (kx_sh/7e7)**(1/9.6)

        J_sh = pc2*np.sqrt(kx_sh/Phi_sh)

        ## Imbibition
        sgrmax = snc + 0.5*(1-snc-swi)
        sgh = sg[sg>sgrmax]
        zgh = sg[sg<=sgrmax] # zeros for reference
        zgh[:] = 0

        #Land trapping model    
        c = 1/(sgrmax - snc) - 1/(sgh[-1]-snc)
        sgrh = snc +(sgh[-1]-snc)/(1+c*(sgh[-1]-snc))
        sgf = snc + 0.5*((sgh-sgrh)+np.sqrt((sgh-sgrh)**2 + (4/c)*(sgh-sgrh)))
        
        # Relative permeability
        swf = 1 - sgf
        
        swfh = (swf-swi)/(1-swi-snc)

        krgf = 0.4*(1-swfh**2)*(1-swfh)**2
        krgf[krgf<0] = 0

        # Capillary pressure
 
        if (cap_p_model == 'MBC' or cap_p_model == 'BC'): # Pini and Benson, 2017
            swi_hat = swi + 0.01
            swih = (swi_hat-swi)/(1-swi-sgrmax)
            pci = pce_1*(1-swih**0.5)**-1
            pcI = pci*((swfh**-0.5)-1)
            pci2 = pce_2*(1-swih**0.5)**-1
            pcI2 = pci2*((swfh**-0.5)-1)
    
        elif cap_p_model == 'VG':
                
                pcI = pce_1*(((swfh**-1/0.464)-1)**(1/1.864)-1)
                pcI2 = pce_2*(((swfh**-1/0.559)-1)**(1/2.266)-1)

        jsI = pcI*np.sqrt(kx_s/Phi_s)
        jsI[-1] = J_s[0]
        pcI[-1] = pc1[0]
        pcI = np.flip(np.concatenate((zgh,pcI)))
        jsI = np.flip(np.concatenate((zgh,jsI)))
        jshI = pcI2*np.sqrt(kx_sh/Phi_sh)
        jshI[-1] = J_sh[0]
        jshI = np.flip(np.concatenate((zgh,jshI)))
        krgI = np.flip(np.concatenate((zgh,krgf)))
        
        wat_sat[0] = swi

        # Save to txt file using numpy
        zeros = np.zeros(shape = (len(wat_sat),))
        data1 = np.array([wat_sat,krl,zeros,J_s,jsI])
        data2 = np.array([wat_sat,krl,zeros,J_sh,jshI])
        data3 = np.array([sg,krg[::-1],zeros])
        data3[data3 < 0] = 0


        ax0.plot(wat_sat,krg,c=colors[crit_c02_list.tolist().index(snc)],label = labels[crit_c02_list.tolist().index(snc)])
        ax0.set_xlabel('Brine saturation')
        ax0.set_ylabel('Relative permeability')
        ax0.margins(y=0)
        
        ax1.plot(wat_sat,pc1,c=colors[crit_c02_list.tolist().index(snc)],label = labels[crit_c02_list.tolist().index(snc)])
        ax1.set_xlabel('Brine saturation')
        ax1.set_ylabel('Sand Capillary pressure [psi]')
        ax1.margins(y=0)
        ax1.legend()
        ax1.set_ylim(0,20)

        if print_bool == True:
            np.savetxt(path_to_save + '/rel_perm_1_' + str(round(snc,2)) + '.txt',data1.T,fmt=['%.6f','%.6f','%.6f','%.6f','%.6f'])
            np.savetxt(path_to_save + '/rel_perm_2_' + str(round(snc,2)) + '.txt',data2.T,fmt=['%.6f','%.6f','%.6f','%.6f','%.6f'])
            np.savetxt(path_to_save + '/sg_'+ str(snc) + '.txt',data3.T,fmt=['%.6f','%.6f','%.6f'])
        
    ax2.plot(wat_sat,krgI,'k',linestyle='--',label = 'Imbibition relative permeability')
    ax2.plot(wat_sat,krg,'b',label = 'Drainage relative permeability')
    ax2.plot(wat_sat,krl,'b',label = 'Liquid relative permeability')
    ax2.set_xlabel('Brine saturation')
    ax2.set_ylabel('Relative permeability')
    ax2.margins(y=0)
    ax2.legend()

    ax3.plot(wat_sat,pcI,'k',linestyle='--',label = 'Imbibition capillary pressure')
    ax3.plot(wat_sat,pc1,'b',label = 'Drainage capillary pressure')
    ax3.set_xlabel('Brine saturation')
    ax3.set_ylabel('Capillary pressure [psi]')
    ax3.margins(y=0)
    ax3.legend()
    ax3.set_ylim(0,20)


    ax0.plot(wat_sat,krl,'pink',linestyle='--',label = 'Liquid relative permeability')
    ax0.legend()
    plt.tight_layout()
    plt.show()

    return 

        
