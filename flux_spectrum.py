import numpy as np
from scipy import constants
from scipy.integrate import quad

# Boltzmann constant in eV/K
k = constants.value('Boltzmann constant in eV/K')

class Flux() :
    """
    This class evaluates the neutron spectrum.  The thermal cutoff is treated
    as a variable parameter to ensure a specific fast-to-thermal ratio.
    
    At thermal energy range (e < e1 eV), the flux is approximated by Maxwellian 
    distribution (D&H book Eq.(9-6)). 
    
    At fast energy range (e2 MeV < e < 20MeV), the flux is approximated by U-235 
    chi-spectrum (D&H book Eq.(2-112)).
    
    At epithermal energies (e1 eV < e < e2 MeV), flux = 1/e
    
    
    
    ratio : fast-to-thermal flux ratio
    """
    def __init__(self, ratio = 2, thermal_temp = 600.0):
        self.e2 = 1e6
        self.thermal_temp = thermal_temp
        
        # Maxwellian distribution, Eq.(9-6)
        self.m = lambda x : x ** 0.5 * np.exp(-x/(k*self.thermal_temp))

        # U235 chi distribution, Eq.(2-112)
        self.chi = lambda x : np.exp(-1.036e-6*x)*np.sinh((2.29e-6 * x)**0.5)     
                         
        # Middle energy range
        self.f = lambda x : 1 / x
        
        # Compute ratio as a function of thermal cutoff
        E = np.logspace(-4, 0.1, 200)
        R = np.array([self.compute_ratio(e1) for e1 in E])
        
        # Compute thermal cutoff for given ratio
        self.e1 = np.interp(1.0/ratio, R, E)
        print 'Thermal cutoff is {} eV'.format(self.e1)

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
        
    def evaluate(self, e):
        # Evaluate flux at Energy e in eV
        return  (e<=self.e1)             * self.c1*self.m(e) + \
                (e>self.e1)*(e<=self.e2) * (self.c2 / e)     + \
                (e>self.e2)              * self.c3*self.chi(e)

if __name__ == "__main__" :

    import matplotlib.pyplot as plt
    from multigroup_utilities import *
    from nice_plots import init_nice_plots
    init_nice_plots()
    
    # PWR-like and TRIGA-like spectra
    pwr = Flux(7.0, 600.0)
    triga = Flux(1.0, 600.0)
        
    # Evaluate the flux
    E = np.logspace(-5, 7, 1000)
    phi_pwr = pwr.evaluate(E)
    phi_triga = triga.evaluate(E)
   
    # Collapse the flux and flux per unit lethargy onto WIMS 69-group structure   
    bounds = energy_groups(structure='wims69')
    pwr_mg = collapse(bounds, phi_pwr, E)
    triga_mg = collapse(bounds, phi_triga, E)
    phi_mg_pul = collapse(bounds, E*phi_pwr, E)
    triga_mg_pul = collapse(bounds, E*phi_triga, E)
    
    # Produce step-plot data for each spectrum    
    E_mg, phi_mg = plot_multigroup_data(bounds, pwr_mg) 
    _, triga_mg = plot_multigroup_data(bounds, triga_mg) 
    _, phi_mg_pul = plot_multigroup_data(bounds, phi_mg_pul) 
    _, triga_mg_pul = plot_multigroup_data(bounds, triga_mg_pul) 
    
    plt.figure(1)
    plt.loglog(E, phi_pwr, 'k', E_mg, phi_mg, 'k--', 
               E, phi_triga, 'b', E_mg, triga_mg, 'b:')
    plt.xlabel('$E$ (eV)')
    plt.ylabel('$\phi(E)$')
    
    plt.figure(2)
    plt.loglog(E, E*phi_pwr, 'k', E_mg, phi_mg_pul, 'k--', 
               E, E*phi_triga, 'b', E_mg, triga_mg_pul, 'b:')
    plt.xlabel('$E$ (eV)')
    plt.ylabel('$E\phi(E)$')
    
    plt.show()