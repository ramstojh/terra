import numpy as np
import numpy
import matplotlib as mpl
import matplotlib.pyplot as pyplot
import math
from scipy.optimize import curve_fit
from scipy.stats import chi2

from scipy import stats
from numpy import*

#---------------------------------------------------------------------------------------------
#Input values
Fe     = 0.04                               # stellar metallicity [Fe/H]
m      = 0                                  # differential abundance, 0 for [X/H] and 25 for [X/Fe]
d      = 0                                  # 0 for twin-sun and 1 for sun-twin
ch     = 3                                  # chondrites, CI=2, CM=3, CO=4, CV=5, H=6, L=7, LL=8, EH=9, EL=10
Nmet   = 0.5                                # meteoritic-like mass in units of terrestrial mass
Nearth = 3.5                                # earth-like mass in units of terrestrial mass
Mcon   = 0.017                              # convective mass in units of terrestrial mass
Msun   = 1.989*10**(30)                     # solar mass (kilograms)
Mear   = 5.972*10.0**(24)                   # terrestrial mass (kilograms)
HH     = 1.00794*pow(10, 12)                # atomic weight of Hydrogen 
HHe    = 4.00260*pow(10, 10.93)             # atomic weight of Helium
norabu = -0.019                             # normalize the abundances to some element (almost all cases Carbon)
Norm   = 1                                  # normalize the mean abundance of 11 solar twins (0 yes, 1 no)
#---------------------------------------------------------------------------------------------
# reading data
s1 = "./TABLES/input-data.txt"   
s2 = "./TABLES/chondrites.txt"            
s3 = "abundances_HD45184.txt"            

file1 = open(s1, "r")
anum  = np.genfromtxt(s1, skip_header=3, usecols=[2], dtype="f8")     
earth = np.genfromtxt(s1, usecols=[3], dtype="f8")                    
tcon  = np.genfromtxt(s1, usecols=[4], dtype="f8")                    
abud = np.genfromtxt(s1, skip_header=3, usecols=[5], dtype="f8")     

file2 = open(s2, "r")
CMmet = np.genfromtxt(s2, usecols=[ch], dtype="g")                    

file3 = open(s3, "r")
yourelement   = np.genfromtxt(s3, skip_header=1, usecols=[0], dtype="str")
Z             = np.genfromtxt(s3, skip_header=1, usecols=[1], dtype="f8")
yourabundance = np.genfromtxt(s3, skip_header=1, usecols=[2], dtype="f8")
erroabundance = np.genfromtxt(s3, skip_header=1, usecols=[3], dtype="f8")

#---------------------------------------------------------------------------------------------
element = range(1, 93)
wat = [anum[i]*10.0**(abud[i]+Fe) for i in range(len(abud))] 

SH = sum(wat)
SHT = SH + HH + HHe                                                   

Twat = [HH, HHe]                       
Twat.extend(wat)                                                      

a=Mcon*Msun/SHT
CMS = [a*Twat[i] for i in range(len(Twat))]                           

SMet = sum(CMmet)
Mass_CM  = [Mear*CMmet[i]/SMet for i in range(len(CMmet))]            

searth = sum(earth)
mass_earth = [Mear*earth[i]/searth for i in range(len(earth))]      

EACM = [log10(1 + (Nmet*Mass_CM[i] + Nearth*mass_earth[i])/CMS[i]) for i in range(len(CMmet))]

if d == 0:
  XH = [(EACM[i] - EACM[m]) for i in range(len(EACM))]
elif d == 1:
  XH = [-(EACM[i] - EACM[m]) for i in range(len(EACM))]

XH.insert(43, 0)
XH.insert(61, 0)
XH.insert(84, 0)
XH.insert(85, 0)
XH.insert(86, 0)
XH.insert(87, 0)
XH.insert(88, 0)
XH.insert(89, 0)
XH.insert(91, 0)

difference = list(set(element) - set(Z))
realdiff = [difference[i] - 1 for i in range(len(difference))]
abundred = delete(XH, tuple(realdiff), axis=0)                    
cond = list(tcon)

cond.insert(43, 0)
cond.insert(61, 0)
cond.insert(84, 0)
cond.insert(85, 0)
cond.insert(86, 0)
cond.insert(87, 0)
cond.insert(88, 0)
cond.insert(89, 0)
cond.insert(91, 0)

Tc = delete(cond, tuple(realdiff), axis=0)                        

#---------------------------------------------------------------------------------------------
#solar twins
xvol = np.arange(1, 1251)
xref = np.arange(1253, 1801)
if Norm == 0:
  yvol = (-0.052   + 2.4E-05*xvol)
  yref =  (-0.14937 + 1.01669E-04*xref)
elif Norm == 1:
  yvol = (-0.052   + 2.4E-05*xvol) + 0.052
  yref =  (-0.14937 + 1.01669E-04*xref) + 0.052

#X2-test
#---------------------------------------------------------------------------------------------
yourabundanceN      = [yourabundance[i] - norabu for i in range(len(yourabundance))]
chisqr = sum([(yourabundanceN[k] - abundred[k])**2 / erroabundance[k]**2 for k in range(len(abundred))])
dof = len(abundred) - 1
Xred = chisqr / dof

print '====================================='
print 'The chi squared reduced is:', '%0.2f' %Xred
print '====================================='

#set global settings
#---------------------------------------------------------------------------------------------
def init_plotting():
    pyplot.rcParams['figure.figsize'] = (15, 9)
    pyplot.rcParams['font.size'] = 10           
    pyplot.rcParams['axes.labelsize'] = pyplot.rcParams['font.size']
    pyplot.rcParams['axes.titlesize'] = 1.5*pyplot.rcParams['font.size']
    pyplot.rcParams['legend.fontsize'] = pyplot.rcParams['font.size']
    pyplot.rcParams['xtick.labelsize'] = pyplot.rcParams['font.size']
    pyplot.rcParams['ytick.labelsize'] = pyplot.rcParams['font.size']
    pyplot.rcParams['savefig.dpi'] = 40*pyplot.rcParams['savefig.dpi']
    pyplot.rcParams['xtick.major.size'] = 9
    pyplot.rcParams['xtick.minor.size'] = 5
    pyplot.rcParams['xtick.major.width'] = 1.2
    pyplot.rcParams['xtick.minor.width'] = 1.2
    pyplot.rcParams['ytick.major.size'] = 9
    pyplot.rcParams['ytick.minor.size'] = 5
    pyplot.rcParams['ytick.major.width'] = 1.2
    pyplot.rcParams['ytick.minor.width'] = 1.2
    pyplot.rcParams['legend.frameon'] = True
    pyplot.rcParams['legend.loc'] = 'best'
    pyplot.rcParams['axes.linewidth'] = 2.0

init_plotting()

#---------------------------------------------------------------------------------------------
#plotting results
figure = 'Abundance-vs-Tc_HD45184.pdf'
pyplot.plot(Tc, yourabundanceN, linestyle=" ", marker="o", color="#ec7063", markersize=14, mew=1.3, label=r'$\mathrm{Observed \ abundances}$')
pyplot.errorbar(Tc, yourabundanceN, yerr=erroabundance, linestyle="None", marker="None", color="#ec7063", elinewidth=2.)
pyplot.plot(Tc, abundred, linestyle=" ", marker="^", color="blue", markersize=12, label =r'$\mathrm{ Predicted \ abundances}$')
pyplot.plot(xvol, yvol, linestyle="--", marker=" ", color="black", linewidth = 2.2, label=r'$\mathrm{Twins\ - \ Sun \ (Melendez\ et\ al.\ 2009)}$')
pyplot.plot(xref, yref, linestyle="--", marker=" ", color="black", linewidth = 2.2)
pyplot.legend(loc='best', numpoints=1, prop={'size':19}, shadow=True)
pyplot.xlabel(r'$\mathrm{Condensation\ Temperature (K)}$', color='black', fontsize=28)
pyplot.ylabel(r'$\Delta \mathrm{[X/C]} \mathrm{(dex)}$', fontsize=28)
pyplot.minorticks_on()
pyplot.tick_params(axis='both', labelsize=19)
pyplot.tight_layout()
pyplot.savefig(figure)
pyplot.clf()

#---------------------------------------------------------------------------------------------
# saving results
file4 = open("Output.txt", "w")
file4.write("#   XC   " + "   XO " + "  error  " + "   TC   " + "   element   "  )
file4.write("\n")
for k in range(len(abundred)):
    file4.write("  " +  "%4.4f" %abundred[k] + "  " + "%4.4f" %yourabundance[k] + "  " + "%4.4f" %erroabundance[k] + "  " + "%4.1f" %Tc[k] + "        " + yourelement[k] + "\n")
file4.close()

