import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
from scipy.integrate import quad
from scipy.optimize import minimize
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
    At thermal energy range (e < 1eV), the flux is approximated by Maxwellian 
    distribution (D&H book Eq.(9-6)). 
    At fast energy range (0.1MeV < e < 20MeV), the flux is approximated by U-235 
    chi-spectrum (D&H book Eq.(2-112)).
    At epithermal energies (1eV < e < 0.1MeV), flux = 1/e
    
    
    r : thermal-to-fast flux ratio
    """
    def __init__(self, r = 2, thermal_t = 600.0):
        self.e2 = 1e5
        self.thermal_t = thermal_t
        
        # Maxwellian distribution, Eq.(9-6)
        self.m = lambda x : x ** 0.5 * np.exp(-x / (k * self.thermal_t))

        # U235 chi distribution, Eq.(2-112)
        self.chi = lambda x : np.exp(-1.036e-6 * x) * np.sinh((2.29e-6 * x) ** 0.5)     
                         
        # Middle energy range
        self.f = lambda x : 1 / x
        self.r = 0
        E = np.logspace(-5,1,1000)
        R = np.array([self.balance(e) for e in E])
        plt.loglog(E, R)
        plt.xlabel('Energy cutoff [eV]')
        plt.ylabel('thermal to non-thermal flux ratio, r')
        plt.savefig('ratio_to_cutoff.pdf')
        plt.clf()
        self.e1 = np.interp(r, R, E)
        print 'Thermal cutoff is {} eV'.format(self.e1)
        self.r = r
        self.c1 = 1.0
        self.c2 = self.m(self.e1) / self.f(self.e1)
        self.c3 = self.c2 * self.f(self.e2) / self.chi(self.e2)
        
    def balance(self, x):
        A = quad(self.m, 0, x)[0]
        B = self.m(x) / self.f(x) * quad(self.f, x, self.e2)[0]
        C = self.m(x) / self.f(x) * self.f(self.e2) / self.chi(self.e2) * quad(self.chi, self.e2, 2e7)[0]
        Q = A / (B + C) - self.r
        return abs(Q)
        
        
    def compute_flux(self, e):
        """
        Evaluate flux at Energy e.
        e : neutron energy in eV
        """
        # thermal
        if e <= self.e1:
            return self.c1 * self.m(e)
        # epithermal
        elif self.e1 < e <= self.e2:
            return self.c2 / e
        elif e >= self.e2:
            return self.c3 * self.chi(e)

# ratio of thermal flux to fast flux = 1e-5
# neutron temperature = 600K
f = Flux(2, 600.0)
e = np.logspace(-2, 7, 1000)
flux = np.zeros(len(e))
for i in range(len(e)):
    # f.compute_flux can return flux at a specific energy e (eV).
    flux[i] = f.compute_flux(e[i])
    
plt.loglog(e, flux)
plt.xlabel('E (eV)')
plt.ylabel('flux')
plt.savefig('flux_spectrum')
    
