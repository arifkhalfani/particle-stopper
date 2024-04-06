#Kinetic Energy vs Position Graph of a Particle Traversing Graphite

import numpy as np
import math
import matplotlib.pyplot as plt

#conversion
MeV = 1.6*10**-13
eV = 1.6*10**-19

#constants
c = 299792458                       #speed of light
e = 1.6 * 10**-19                   #elementary charge
E_e = 0.5110*MeV                    #electron rest energy
Z = 6                               #carbon atomic number
A = 12                              #carbon relative mass
Na = 6.022 * 10**23                 #Avogadro's constant
eps_0 = 8.85 * 10**-12              #vacuum permittivity
rho = 2260                          #density of graphite
I = 78*eV                           #carbon mean excitation energy
alpha = 1/137                       #fine structure constant
c1 = (e**2/(4*math.pi*eps_0))**2

#settings
E0 = 3000*MeV                       #energy of one particle (3000 for electron beam, 15000 for proton beam)
m = 0.511*MeV/(c**2)                #mass of particle (0.511 for electrons, 938.27 for protons) 
z = 1                               #absolute charge of particle (in e)
x = np.linspace(0, 0.2, 10000)      #distance travelled by particle (in m)

#initial values
E_p = m*c**2                        #rest energy of particle
T0 = E0 - E_p                       #initial kinetic energy
gamma0 = T0/(E_e) + 1               #lorentz factor
beta0 = math.sqrt(1-1/gamma0**2)    #v/c

T = T0
gamma = gamma0
beta = beta0

T_values= np.zeros_like(x)

for i, n in enumerate(x):
    dx = x[1] - x[0]

    #Coulomb interactions
    dE_c = -1*rho*c1*(2*math.pi*(z**2)*Na*Z*rho/(E_e*A*beta**2))*(
            math.log(T*(T+E_e**2)**2*beta**2/(2*I**2*E_p)) + (1-beta**2) -
            (2*math.sqrt(1-beta**2)-1+beta**2)*math.log(2) + 1/8*(1-math.sqrt(1-beta**2))**2)
    
    #Bremsstrahlung
    dE_b = -1*rho*c1*(alpha*Z**2*Na*(T+E_e)*rho/(E_e**2*A))*(
            4*math.log(2*(T+E_e**2)/E_p) - 4/3)
    
    #Total dE/dx
    dE_total = dE_c + dE_b

    T += dE_total * dx
    gamma = T/(E_e) + 1
    beta = math.sqrt(1-1/gamma**2)
    T_values[i] = T/MeV  #Change units to MeV
    
#Plotting
plt.figure(figsize=(10, 4))
plt.plot(x, T_values)
plt.xlabel('Position (m)')
plt.ylabel('Kinetic Energy (MeV)')
plt.title('Kinetic Energy vs Position')
plt.grid()
plt.show()
