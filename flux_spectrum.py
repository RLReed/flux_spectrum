import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import constants
from scipy.integrate import quad
from matplotlib import rc
rc('font',**{'family':'serif'})
from matplotlib import rcParams
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'
rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18
rcParams['lines.linewidth'] = 1.85
rcParams['axes.labelsize'] = 20
rcParams.update({'figure.autolayout': True})

# Boltzmann constant in eV/K
k = constants.value('Boltzmann constant in eV/K')

class Flux():
    """
    This class evaluate neutron flux at energy e. 
    At thermal energy range (e < e1 eV), the flux is approximated by Maxwellian 
    distribution (D&H book Eq.(9-6)). 
    At fast energy range (e2 MeV < e < 20MeV), the flux is approximated by U-235 
    chi-spectrum (D&H book Eq.(2-112)).
    At epithermal energies (e1 eV < e < e2 MeV), flux = 1/e
    
    ratio : thermal-to-fast flux ratio
    """
    def __init__(self, ratio = 2, thermal_temp = 600.0):
        self.e2 = 1e6
        self.thermal_temp = thermal_temp
        
        # Maxwellian distribution, Eq.(9-6)
        self.m = lambda x : x ** 0.5 * np.exp(-x / (k * self.thermal_temp))

        # U235 chi distribution, Eq.(2-112)
        self.chi = lambda x : np.exp(-1.036e-6 * x) * np.sinh((2.29e-6 * x) ** 0.5)     
                         
        # Middle energy range
        self.f = lambda x : 1 / x
        
        # Compute ratio as a function of thermal cutoff
        E = np.logspace(-4, 0.1, 200)
        R = np.array([self.compute_ratio(e1) for e1 in E])
        
        # Create plot of ratio vs cutoff energy
        plt.loglog(R, E)
        plt.ylabel('Thermal Energy Cutoff [eV]')
        plt.xlabel('thermal to non-thermal flux ratio, r')
        plt.savefig('ratio_to_cutoff.pdf')
        plt.clf()
        
        # Compute thermal cutoff for given ratio
        self.e1 = np.interp(ratio, R, E)
        print 'Thermal cutoff is {} eV'.format(self.e1)
        self.ratio = ratio
        
        # Compute constants for each part of the spectrum
        self.c1 = 1.0
        self.c2 = self.m(self.e1) / self.f(self.e1)
        self.c3 = self.c2 * self.f(self.e2) / self.chi(self.e2)
        
    def compute_ratio(self, e1):
        A = quad(self.m, 0, e1)[0]
        C2 = self.m(e1) / self.f(e1)
        C3 = self.f(self.e2) / self.chi(self.e2)
        B = C2 * quad(self.f, e1, self.e2)[0]
        C = C2 * C3 * quad(self.chi, self.e2, 2e7)[0]
        r = A / (B + C)
        return r
        
    def compute_flux(self, e):
        # Evaluate flux at Energy e in eV
        # thermal
        return self.c1 * self.m(e) if e <= self.e1 else (self.c2 / e if e <= self.e2 else self.c3 * self.chi(e))

# ratio of thermal flux to fast flux = 1e-5
# neutron temperature = 600K
f = Flux(1e-5, 600.0)
E = np.logspace(-5, 7, 1000)
flux = np.zeros(E.shape)

# Compute the flux at all energies
for i, e in enumerate(E):
    flux[i] = f.compute_flux(e)
lethergy = E * flux
    
# Plot the flux
plt.loglog(E, flux)
plt.xlabel('E (eV)')
plt.ylabel('flux')
plt.savefig('flux_spectrum.pdf')
plt.clf()

# Plot the lethergy
plt.loglog(E, lethergy)
plt.xlabel('E (eV)')
plt.ylabel('$E\phi(E)$')
plt.savefig('lethergy.pdf')
    
