"""
A script for reproducing the observed abundances of a star from the solar abundances of Asplund + 2021.
"""

import pandas as pd
import numpy as np
from astropy import constants as const
import matplotlib.pyplot as plt
from tqdm import trange, tqdm
import pkg_resources
import os

pd.options.display.float_format = "{:,.2f}".format


Msun   = const.M_sun.value                  # solar mass (1.9884099×10^30 kg)
Mearth = const.M_earth.value                # terrestrial mass (5.9721679×10^24 kg)
HH     = 1.00794*pow(10, 12)                # atomic weight of Hydrogen 
HHe    = 4.00260*pow(10, 10.93)             # atomic weight of Helium

#configuring plot
plotpar = {'axes.labelsize': 25,
           'xtick.labelsize': 20,
           'ytick.labelsize': 20,
           'font.family': 'serif', 
           'axes.linewidth': 2,
           'text.usetex': True}
plt.rcParams.update(plotpar)

#data path
this_dir, this_filename = os.path.split(__file__)

class pacha(object):
    def __init__(self):
        self.feh         = None
        self.mass        = None
        self.model       = None
        self.chondrites  = None
        self.data_input  = None
        self.Mcon        = None
        self.data_output = None
    
    @classmethod
    def cvmass(self, feh, mass, model=None):
        '''
        Calculate the convective mass of a star given a known [Fe/H] and mass.
        Important to note that Z is scaled to [Fe/H] in the model tables.
        The output is in solar mass units.
        '''
        
        '''
        Args
            [Fe/H] : stellar metallicity
            Mass   : stellar mass
            model  : stellar evolution model
        '''
        
        #choosing model
        if not model:
            model = 'yale.txt'
        if model == 'yale':
            model = 'yale.txt'
        if model == 'lionel':
            model = 'model-lionel.txt'
            
        #reading data
        PATH_model = os.path.join(this_dir, 'data', model)
        data = np.loadtxt(PATH_model)
        
        #double interpolation between mass and metallicity
        for i in range(len(data[0])):
            if feh == data[0][i]:
                Z1 = feh
                Z2 = 0.0
                xp = i
                break
            elif data[0][i-1] < feh < data[0][i]:
                Z1 = data[0][i-1]
                Z2 = data[0][i]
                xp = i

            for j in range(10):
                if mass == data[j][0]:
                    M1 = mass
                    M2 = 0.0
                    yp = j
                    break
                elif data[j-1][0] < mass < data[j][0]:
                    M1 = data[j-1][0]
                    M2 = data[j][0]
                    yp = j
                
        m1 = (Z2 - feh)/(Z2 - Z1)
        m2 = (feh - Z1)/(Z2 - Z1)
        b1 = (M2 - mass)/(M2 - M1)
        b2 = (mass - M1)/(M2 - M1)
        
        if Z1 == feh and M1 == mass:
            CM = data[yp][xp]
        elif Z1 == feh and M1 == data[yp-1][0]:
            CM11 = data[yp-1][xp]
            CM12 = 0
            CM21 = data[yp][xp]
            CM22 = 0
            CM = (m1*CM11 + m2*CM12)*b1 + (m1*CM21 + m2*CM22)*b2
        elif Z1 == data[0][xp-1] and M1 == mass:
            CM11 = data[yp][xp-1]
            CM12 = data[yp][xp]
            CM21 = 0
            CM22 = 0
            CM = (m1*CM11 + m2*CM12)*b1 + (m1*CM21 + m2*CM22)*b2
        elif Z1 == data[0][xp-1] and M1 == data[yp-1][0]:
            CM11 = data[yp-1][xp-1]
            CM12 = data[yp-1][xp]
            CM21 = data[yp][xp-1]
            CM22 = data[yp][xp]
            CM = (m1*CM11 + m2*CM12)*b1 + (m1*CM21 + m2*CM22)*b2
            
        return CM

    
    @classmethod
    def mod_abd(self, feh, mass, data_input=None, model=None, Mcon=None, chondrites=None, data_output=None):
        '''
        Args
            data_input  : table containing the observed abundances
            Mcon        : convective mass
            chondrites  : chondrites abundances
            data_output : name of the star
            
        '''
        
        if not Mcon:
            Mcon = self.cvmass(feh, mass, model)
        else:
            Mcon = Mcon
        
        if not chondrites:
            chondrites = 'CM'
        else:
            chondrites = chondrites
            
        if not data_output:
            data_output = ''
        else:
            data_output = data_output
        
        #solar and earth abundances
        PATH_seab = os.path.join(this_dir, 'data', 'seab_data.csv')
        seab      = pd.read_csv(PATH_seab)
        
        #chondrites
        #the user can choose the chondrites type: CI, CM, CO, CV, H, L, LL, EH, EL
        PATH_met = os.path.join(this_dir, 'data', 'Chondrites.csv')
        met      = pd.read_csv(PATH_met, usecols=[chondrites])
        
        #observed abundances
        tab_ab          = pd.read_csv(data_input)
        carbon          = tab_ab[tab_ab.element == 'C']
        tab_ab['[C/H]'] = tab_ab['[X/H]'] - carbon['[X/H]'].values
        
        #determining abundances of Asplund21
        seab['wat'] = seab['A']*10**(seab['Asplund21']+feh)
        #seab['wat'] = seab['A']*10**(seab['Asplund09']+feh)
        
        seab.at[0,'wat'] = HH
        seab.at[1,'wat'] = HHe
        
        #convective mass
        seab['CMS'] = Mcon*Msun*seab['wat']/np.sum(seab['wat'])
        
        #metheoritic mass
        met['Mass_CM'] = Mearth*met/np.sum(met)
        
        #terrestrial mass
        seab['mass_earth'] = Mearth*seab['Earth']/np.sum(seab['Earth'])
        
        seab['Mass_CM'] = met['Mass_CM']
        
        #getting observed and predicted elements
        merged_tab = pd.merge(tab_ab, seab, how="inner", on=['element'])
        
        total = []
        
        for i in trange(len(merged_tab['Mass_CM']), desc='finding the best solution'):
            eacm1 = []
            TME   = [] # Nmet + Nearth
            for Nmet in np.arange(0., 20, 0.1):
                for Nearth in np.arange(0., 20, 0.1):
                    #equation A.5 in Yana Galarza et al. 2016
                    EACM = np.log10(1 + (Nmet*merged_tab['Mass_CM'][i] + Nearth*merged_tab['mass_earth'][i])/merged_tab['CMS'][i])
                    eacm1.append(EACM)
                    TME.append([Nmet, Nearth])
                           
            #chi2
            chi = (merged_tab['[C/H]'][i] - eacm1)**2/merged_tab['err_[X/H]'][i]**2
            total.append(chi)
            

        chisq  = []
        number = []
        for i in range(len(total[0])):
            soma = []
            for j in range(len(merged_tab['Mass_CM'])):
                soma.append(total[j][i])
            
            chi  = np.sum(soma)
            dof  = len(merged_tab['[C/H]']) - 1
            Xred = chi / dof
            chisq.append(Xred)
            #print (Xred)
        
        TTME  = []
        alpha = []
        beta  = []
        for i in range(len(TME)):
            TTME.append(np.sum(TME[i]))
            alpha.append(TME[i][0])
            beta.append(TME[i][1])
            
        #placing all in data frame
        data = {'chi2':chisq, 'Nmet': alpha, 'Nearth': beta, 'Total_mass': TTME}
        df = pd.DataFrame(data)
        
        #getting the best value given the minimum chi2
        conv = df[df.chi2 == df['chi2'].min()]
        if data_output == '':
            conv.to_csv('chi2_convolution.csv', index=False)
            print (conv.to_string(index=False))
        else:
            conv.to_csv('chi2_convolution_'+data_output+'.csv', index=False)
            print (conv.to_string(index=False))
        
        
        #plotting results
        plt.figure(figsize=(7, 6))
        plt.scatter(conv['Total_mass'], conv['chi2'], facecolors='none', edgecolors='r')
        plt.axhline(y=conv['chi2'].values, c='r', alpha=0.5)
        plt.axvline(x=conv['Total_mass'].values, c='r', alpha=0.5)
        plt.scatter(df['Total_mass'], df['chi2'], s=2)
        plt.xlabel("Total rocky mass (M$_{\oplus}$)")
        plt.ylabel("Chi2")
        plt.ylim(conv['chi2'].values-10, conv['chi2'].values+50)
        plt.xlim(conv['Total_mass'].values-5, conv['Total_mass'].values+5)
        plt.tight_layout()
        if data_output == '':
            plt.savefig('chi2_test.png', dpi=200)
        else:
            plt.savefig('chi2_test_'+data_output+'.png', dpi=200)

                
        #getting abundances with the best solution
        merged_tab['model [C/H]'] = np.log10(1 + (conv['Nmet'].values*merged_tab['Mass_CM'] + conv['Nearth'].values*merged_tab['mass_earth'])/merged_tab['CMS'])
        merged_tab['model [X/H]'] = merged_tab['model [C/H]'] + tab_ab['[X/H]'].values
        header = ['Z', 'Tcond', 'element', '[X/H]', 'err_[X/H]', '[C/H]', 'model [C/H]', 'model [X/H]']
        if data_output == '':
            merged_tab.to_csv('results.csv', columns=header, index=False)
        else:
            merged_tab.to_csv('results_'+data_output+'.csv', columns=header, index=False)
        
        #plotting model vs observed abundances
        plt.figure(figsize=(10, 5))
        plt.scatter(merged_tab['Tcond'], merged_tab['[C/H]'], edgecolors='blue', s=170, marker="s", facecolors='none', linewidths=3., alpha=0.8, zorder= 1, label='Observed abundance')
        plt.errorbar(merged_tab['Tcond'], merged_tab['[C/H]'], yerr=merged_tab['err_[X/H]'], linestyle='None', marker='None', color='blue', elinewidth=2.5, alpha=0.7)
        plt.scatter(merged_tab['Tcond'], merged_tab['model [C/H]'], s=180, c='red', linewidths=1., alpha=0.5, zorder=1, label='Predicted abundance')
        plt.legend(loc=2, numpoints=1, prop={'size':12}, shadow=True)
        plt.xlabel(r'$\mathrm{Condensation\ Temperature (K)}$', color='black')
        plt.ylabel(r'$\Delta \mathrm{[X/C]} \mathrm{(dex)}$')
        plt.tight_layout()
        if data_output == '':
            plt.savefig('model_vs_observed_abundances.png', dpi=200)
        else:
            plt.savefig('model_vs_observed_abundances_'+data_output+'.png', dpi=200)



