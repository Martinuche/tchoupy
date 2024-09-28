""" List of many values about physics constants and various orders of magnitude.
"""

print('BAD MODULE, DO NOT USE YET: prior to revamping of the module the '
      'dimension was indicated, and the minus sign in the unit was implicit.'
      ' I have removed the dimension however I did not reintroduce the minus '
      'sign when necessary, so many variables are describing 1 over what they '
      'are actually supposed to describe.')

__authors__ = ("Martin Teuscher",)
__contact__ = ("teuscher.edu@gmail.com",)
__version__ = "1.0"
__date__ = "2024/09"

import numpy as np
import scipy.constants as const
from ..systems import HEP,SI,Imperial,Cosmo
# =============================================================================
# SOURCES
# =============================================================================

"""Source of scipy.constants : CODATA 2018
https://docs.scipy.org/doc/scipy/reference/constants.html
wall chart of constants accessible at https://physics.nist.gov/cuu/Constants/index.html

OTHER SOURCES:
"""

__SOURCES = {
    0 : 'https://docs.scipy.org/doc/scipy/reference/constants.html \nSee also: https://physics.nist.gov/cuu/Constants/index.html',
    1 : 'https://web.archive.org/web/20131110215339/http://asa.usno.navy.mil/static/files/2014/Astronomical_Constants_2014.pdf',
    2 : 'https://en.wikipedia.org/wiki/Sun',
    3 : 'https://en.wikipedia.org/wiki/Earth',
    4 : 'https://en.wikipedia.org/wiki/Milky_Way',
    5 : 'my Dark Matter course from M2 ICFP ENS',
    6 : 'https://arxiv.org/pdf/1807.06209.pdf',
    7 : 'my Orders of Magnitude notebook',
    8 : 'https://en.wikipedia.org/wiki/List_of_refractive_indices',
    9 : 'https://arxiv.org/abs/1109.3350',
    10 : 'https://en.wikipedia.org/wiki/Higgs_boson',
    11 : 'https://en.wikipedia.org/wiki/Speed_of_sound',
    12 : 'https://fr.wikipedia.org/wiki/Proxima_Centauri',
    13 : 'https://en.wikipedia.org/wiki/Millimetre_of_mercury',
    14 : 'https://arxiv.org/pdf/1801.04268.pdf'
}

#Function used to create the __doc__ of all constants listed in the file.
def source(identifier):
    return f'\n\nSource: {__SOURCES[identifier]}'

# =============================================================================
# UNITS
# =============================================================================

###dim -1 units###
#1 second [0]
s = HEP(1., 's', doc = 'One second'+source(0))
#1 minute [0]
minute = HEP(1., 'min', doc = 'One minute'+source(0))
#1 hour [0]
hr = HEP(1., 'hr', doc = 'One hour'+source(0))
#1 year [0]
yr = HEP(1., 'yr', doc = 'One year'+source(0))
#1 meter [0]
m = HEP(1., 'm', doc = 'One meter'+source(0))
#1 foot [0]
ft = Imperial(1., 'ft', doc = 'One foot'+source(0))
#1 mile  [0]
mile = Imperial(1., 'mile', doc = 'One mile'+source(0))
#1 angstrom [0]
angstrom = HEP(1., 'angstrom', doc = 'One angstrom'+source(0))
#1 astronomical unit [0]
AU = HEP(1., 'AU', doc = 'One astronomical unit'+source(0))
#1 lightyear [0]
lyr = HEP(1., 'lyr', doc = 'One light-year'+source(0))
#1 parsec [0]
pc = HEP(1., 'pc', doc = 'One parsec'+source(0))
#1 megaparsec [0]
Mpc = HEP(1., 'Mpc', doc = 'One mega parsec'+source(0))
      
###dim 1 units###  
#1 electron volt [0]   
eV = HEP(1., 'eV', doc = 'One electronvolt'+source(0))  
#1 giga electron volt [0]
GeV = HEP(1., 'GeV', doc = 'One giga electronvolt'+source(0)) 
#1 joule [0]
J = HEP(1., 'J', doc = 'One joule'+source(0))
#1 erg [0]
erg = HEP(1., 'erg', doc = 'One erg'+source(0))    
#1 kilogram [0]
kg = HEP(1., 'kg', doc = 'One kilogram'+source(0))
#1 gram per mol [0]
gpermol = HEP(1., 'gpermol', doc = 'One gram per mole'+source(0))
#1 kelvin [0]
K = HEP(1., 'K', doc = 'One kelvin'+source(0))
#1 ampère [0]
A = HEP(1., 'A', doc = 'One ampere'+source(0))
#1 volt [0]
V = HEP(1., 'V', doc = 'One volt'+source(0))
#1 hertz [0]
Hz = HEP(1., 'Hz', doc = 'One hertz'+source(0))

###dim -2 units###
#1 square meter [0]
m2 = HEP(1., 'm2', doc = 'One square meter'+source(0))

###dim 2 units###
#1 square giga electronvolt [0]
eV2 = HEP(1., 'eV2', doc = 'One square electronvolt'+source(0))
#1 newton [0]
N = HEP(1., 'N', doc = 'One newton'+source(0))
#1 watt [0]
W = HEP(1., 'W', doc = 'One watt'+source(0))
#1 tesla [0]
T = HEP(1., 'T', doc = 'One tesla'+source(0))
#1 gauss [0]
Gauss = HEP(1., 'G', doc = 'One gauss'+source(0))

###dim -3 units###
#1 cubic meter [0]
m3 = HEP(1., 'm3', doc = 'One cubic meter'+source(0))

###dim 3 units###
#1 cubic giga electronvolt [0]
eV3 = HEP(1., 'eV3', doc = 'One cubic electronvolt'+source(0))
#1 Jansky [0]
Jy = HEP(1., 'Jy', doc = 'One jansky'+source(0))

###dim 4 units###
#1 quartic giga electronvolt [0]
eV4 = HEP(1., 'eV4', doc = 'One quartic electronvolt'+source(0))
#1 gram per cubic centimeter [0]
gpercm3 = HEP(1., 'g1cm_3', doc = 'One gram per cubic centimeter'+source(0))
#1 gram per cubic meter [0]
gperm3 = HEP(1., 'g1m_3', doc = 'One gram per cubic meter'+source(0))
#1 pascal [0]
Pa = HEP(1., 'Pa', doc = 'One pascal'+source(0))
#1 bar [0]
bar = SI(1., 'bar', doc = 'One bar'+source(0))

###dim 0 units###
#1 ohm (dimensionless) [0]
Ohm = HEP(1., 'ohm', doc = 'One ohm'+source(0))
#1 coulomb (dimensionless) [0]
C = HEP(1., 'C', doc = 'One coulomb'+source(0))
#1 mol (dimensionless) [0]
mol = HEP(1., 'mol', doc = 'One mole'+source(0))
#1 meter per second [0]
mpers = HEP(1., 'mpers', doc = 'One meter per second'+source(0))
#1 kilometer per hour [0]
kmpers = HEP(1., 'kmpers', doc = 'One kilometer per second'+source(0))
#1 mile per hour [0]
mileperh = Imperial(1., 'mileperh', doc = 'One mile per hour'+source(0))

# =============================================================================
# FUNDAMENTAL CONSTANTS
# =============================================================================

#Reduced Planck constant [0]
hbar = HEP(1., 'Planck', doc = 'Reduced Planck constant (Planck constant / (2 pi)'+source(0))
#Planck constant [0]
h = HEP(2.*const.pi, 'Planck', doc = 'Planck constant'+source(0))
#Speed of light [0]
lightspeed = HEP(1., 'Planck', doc = 'Speed of light in vacuum'+source(0))
c = HEP(1., 'Planck', doc = 'Speed of light in vacuum'+source(0))
#Boltzmann constant [0]
kB = HEP(1., 'Planck', doc = 'Boltzmann constant'+source(0))
boltzmann = HEP(1., 'Planck' , doc = 'Boltzmann constant'+source(0))
#Gravitational constant [0]
G = HEP(1., 'MPl2', doc = 'Newton gravitational constant'+source(0))
Grav = HEP(1., 'MPl2', doc = 'Newton gravitational constant'+source(0))
#Elementary charge [0]
e = HEP(const.e, 'C', doc = 'Elementary charge'+source(0)) #C is Coulomb
#Fine structure constant [0]
#alpha = HEP(const.alpha, 'one', doc = 'Fine structure constant'+source(0))
"""WHY np.sqrt(alpha) == e.Planck WHILE THERE SHOULD BE A FACTOR sqrt(4pi) OF DIFFERENCE ???"""
alpha = const.alpha
#Fermi coupling constant [0]
GF = HEP(const.value('Fermi coupling constant'), 'GeV2', doc = 'Fermi coupling constant'+source(0))

"""The following two need some units to be defined first"""
#Vacuum permittivity [0] #Other units : F.m^-1 , A^2.s^4.kg^-1.m^-3 , C^2.N^-1.m^-2 , C.V^-1.m^-1 
eps_0 = HEP(const.epsilon_0, 'kg_1m_3s4A2', doc = 'Vacuum dielectric permittivity'+source(0))
#Vacuum permeability [0] #Other units : H.m^-1, N.A^-2
mu_0 = HEP(const.mu_0 * N.eV2 * A.eV**(-2), 'one', doc = 'Vacuum magnetic permeability'+source(0))

#Planck mass [0]
m_Pl = HEP(1., 'MPl', doc = 'Planck mass'+source(0))
#Reduced Planck mass [0]
reduced_m_Pl = HEP(1., 'reducedMPl', doc = 'Reduced Planck mass (Planck mass / sqrt(8 pi))'+source(0))
#Planck time [0]
t_Pl = HEP(1., 'tPl', doc = 'Planck time'+source(0))

# =============================================================================
# NON-FUNDAMENTAL CONSTANTS
# =============================================================================

#Avogadro constant [0]
N_A = const.N_A
#Ideal gas constant  [0]
R = HEP(const.R, 'J1K_1mol_1', doc = 'Ideal gas constant'+source(0))
#1 faraday (charge of a mole of electrons) [7]
#faraday = HEP(N_A*e.C, 'C')
faraday = HEP(const.N_A*e.C, 'C1mol_1', doc = 'One faraday'+source(7))
#1 debye (unit for dipolar moment) [7]
debye = HEP(1e-21 * C.one / c.mpers, 'm', doc = 'One debye'+source(7))


# =============================================================================
# ATOMIC PHYSICS
# =============================================================================

#Electron mass [0]
m_e = HEP(const.m_e, 'kg', doc = 'Electron mass'+source(0))
#Proton mass [0]
m_p = HEP(const.m_p, 'kg', doc = 'Proton mass'+source(0))
#Neutron mass [0]
m_n = HEP(const.m_n, 'kg', doc = 'Neutron mass'+source(0))
#Atomic mass unit [0]
u = HEP(const.u, 'kg', doc = 'Atomic mass unit'+source(0))
#Higgs boson mass [10] #125.25±0.17GeV from https://arxiv.org/pdf/2402.00727.pdf
m_H = HEP(125.0, 'GeV', doc = 'Higgs boson mass'+source(10))
#Rydberg constant [0]
rydberg = HEP(const.Rydberg, 'm', doc = 'Rydberg constant'+source(0))

#Ionization energy of the hydrogen atom [7]
E_ioniz = HEP(13.6, 'eV', doc = 'Ionization energy of the hydrogen atom'+source(7))

#Electron velocity on the 1st Bohr orbit [7]
v_e_Bohr = HEP(1e-2, 'Planck', doc = 'Electron velocity on the 1st Bohr orbit'+source(7))

#Magnetic moment of the hydrogen atom [7]
moment_H = HEP(1e-23 * J.eV / T.eV2, 'eV', doc = 'Magnetic moment of the hydrogen atom'+source(7))

#Ammoniac molecule oscillation frequency [7]
f_ammoniac = HEP(20., 'GHz', doc = 'Ammoniac molecule oscillation frequency'+source(7))

# =============================================================================
# ASTROPHYSICS
# =============================================================================

#Sun
#mass [1]
M_sun = HEP(1., 'Msun', doc = 'Sun mass'+source(1))
#equatorial radius [2]
R_sun = HEP(696342, 'km', doc = 'Sun equatorial radius'+source(2))
#mean density [2]
rho_sun = HEP(1.408, 'g1cm_3', doc = 'Sun mean mass density'+source(2))
#volume [2]
volume_sun = HEP(1.41e18 * 1e9, 'm3', doc = 'Sun volume'+source(2))
#luminosity [2]
L_sun = HEP(3.828e26, 'W', doc = 'Sun luminosity'+source(2))
#velocity of the Sun around the Miky Way [2]
v_sun = HEP(251, 'kmpers', doc = 'Velocity of the Sun around the Milky Way'+source(2))
#distance from Sun to Milky Way center [2]
sun_to_milkyway = HEP(26.66, 'klyr', doc = 'Distance from Sun to Milky Way center'+source(2))

#Earth
#mass [1]
M_earth = HEP(1/332946.0487, 'Msun', doc = 'Earth mass'+source(1))
#equatorial radius [1]
R_earth = HEP(6378.1366, 'km', doc = 'Earth equatorial radius'+source(1))
#mean density [3]
rho_earth = HEP(5.5134, 'g1cm_3', doc = 'Earth mean mass density'+source(3))
#volume [3]
volume_earth = HEP(1.08321e12 * 1e9, 'm3', doc = 'Earth volume'+source(3))
#velocity of the Earth around the Sun [3]
v_earth = HEP(29.7827, 'kmpers', doc = 'Velocity of the Earth around the Sun'+source(3))

#Moon
#mass [1]
M_moon = HEP(1.23000371e-2/332946.0487, 'Msun', doc = 'Moon mass'+source(1))
#equatorial radius [1]
R_moon = HEP(1737.4, 'km', doc = 'Moon equatorial radius'+source(1))

#Mercury
#mass [1]
M_mercury = HEP(1/(6.0236e6), 'Msun', doc = 'Mercury mass'+source(1))
#equatorial radius [1]
R_mercury = HEP(2439.7, 'km', doc = 'Mercury equatorial radius'+source(1))
#Venus
#mass [1]
M_venus = HEP(1/(4.08523719e5), 'Msun', doc = 'Venus mass'+source(1))
#equatorial radius [1]
R_venus = HEP(6051.8, 'km', doc = 'Venus equatorial radius'+source(1))
#Mars
#mass [1]
M_mars = HEP(1/(3.09870359e6), 'Msun', doc = 'Mars mass'+source(1))
#equatorial radius [1]
R_mars = HEP(3396.19, 'km', doc = 'Mars equatorial radius'+source(1))
#Jupiter
#mass [1]
M_jupiter = HEP(1/(1.047348644e3), 'Msun', doc = 'Jupiter mass'+source(1))
#equatorial radius [1]
R_jupiter = HEP(71492, 'km', doc = 'Jupiter equatorial radius'+source(1))
#Saturn
#mass [1]
M_saturn = HEP(1/(3.4979018e3), 'Msun', doc = 'Saturn mass'+source(1))
#equatorial radius [1]
R_saturn = HEP(60268, 'km', doc = 'Saturn equatorial radius'+source(1))
#Uranus
#mass [1]
M_uranus = HEP(1/(2.290298e4), 'Msun', doc = 'Uranus mass'+source(1))
#equatorial radius [1]
R_uranus = HEP(25559, 'km', doc = 'Uranus equatorial radius'+source(1))
#Neptune
#mass [1]
M_neptune = HEP(1/(1.941226e4), 'Msun', doc = 'Neptune mass'+source(1))
#equatorial radius [1]
R_neptune = HEP(24764, 'km', doc = 'Neptune equatorial radius'+source(1))
#Pluto
#mass [1]
M_pluto = HEP(1/(1.36566e8), 'Msun', doc = 'Pluto mass'+source(1))
#equatorial radius [1]
R_pluto = HEP(1195, 'km', doc = 'Pluto equatorial radius'+source(1))

#Io mass [1]
M_io = HEP(4.704e-5*M_jupiter.Msun, 'Msun', doc = 'Io (one of Jupiter galilean moons) mass'+source(1))
#Ganymede mass [1]
M_ganymede = HEP(7.805*M_jupiter.Msun, 'Msun', doc = 'Ganymede (one of Jupiter galilean moons) mass'+source(1))
#Callisto #mass [1]
M_callisto = HEP(5.667e-5*M_jupiter.Msun, 'Msun', doc = 'Callisto (one of Jupiter galilean moons) mass'+source(1))
#Europa mass [1]
M_europa = HEP(2.528e-5*M_jupiter.Msun, 'Msun', doc = 'Europa (one of Jupiter galilean moons) mass'+source(1))
#Titan mass [1]
M_titan = HEP(2.366e-4*M_saturn.Msun, 'Msun', doc = 'Titan (Saturne biggest moon) mass'+source(1))
#Titania mass [1]
M_titania = HEP(4.06e-4*M_uranus.Msun, 'Msun', doc = 'Titania (one of Uranus moons) mass'+source(1))
#Oberon #mass [1]
M_oberon = HEP(3.47e-5*M_uranus.Msun, 'Msun', doc = 'Oberon (one of Uranus moons) mass'+source(1))
#Triton mass [1]
M_triton = HEP(2.089e-4*M_neptune.Msun, 'Msun', doc = 'Triton (one of Neptune moons) mass'+source(1))

#Milky Way
#mass [4]
M_milkyway = HEP(1.15e12, 'Msun', doc = 'Milky Way galaxy mass'+source(4))
#radius [4]
R_milkyway = HEP(26.8/2., 'kpc', doc = 'Milky Way galaxy disk radius'+source(4))
#thickness [4]
width_milkyway = HEP(2.6, 'kpc', doc = 'Milky Way galaxy disk thickness'+source(4))
#velocity of the Milky Way in the local galaxy cluster [5]
v_milkyway = HEP(650, 'kmpers', doc = 'Velocity of the Milky Way galaxy in the Local Galaxy Cluster'+source(5))
#distance to Andromeda galaxy [5] (consistent with [7])
milkyway_to_andromeda = HEP(780., 'kpc', doc = 'Distance from Milky Way galaxy to Andromeda Galaxy'+source(7))

#Solar System radius [7]
R_solarsystem = HEP(1e13, 'm', doc = 'Solar System radius'+source(7))
#Distance to Proxima Centauri [12]
sun_to_proxima = HEP(4.244, 'lyr', doc = 'Distance of Sun to Proxima Centauri star'+source(12))

# =============================================================================
# GEOPHYSICS
# =============================================================================

#Earth gravity [0]
earthgravity = HEP(1., 'earthgravity', doc = 'Earth gravitational acceleration'+source(0))

#differential tide effect [7]
"""CONTEXT ???"""
diff_tide = HEP(1e-6 * earthgravity.earthgravity, 'earthgravity', doc = 'Differential tide effect on Earth [MISSING CONTEXT FOR PROPER DEFINITION]'+source(7))

#angular velocity of Earth proper rotation [7]
f_Earth = HEP(7e-5 / (2*const.pi), 'Hz', doc = 'Frequency of Earth proper rotation (angular velocity / (2 pi))'+source(7))



# =============================================================================
# HIGH ENERGY PHYSICS
# =============================================================================

#1 barn [0]
barn = HEP(1., 'barn', doc = 'One barn'+source(0))

#Maximum number of relativistic degrees of freedom [5]
gstar = 106.75

#Usual mass scale M* for conformal time [9]
Mstar = HEP(np.sqrt(90/(8*(const.pi**3)*gstar)), 'MPl', doc = 'Standard mass scale to define conformal time in cosmology. Mstar = sqrt(90/(8 pi^3 g_*)) M_Planck with g_* = 106.75.'+source(9))

#typical kinetic energy of a proton [5]
"""CONTEXT ???"""
proton_kinetic = HEP(1., 'keV', doc = 'Typical proton kinetic energy [MISSING CONTEXT FOR PROPER DEFINITION]'+source(5))

#electroweak cross section (typical neutrino cross section, ~ g^2 /M_W^4) [5]
"""CONVERSION ISSUE, cf dark matter course"""
EW_sigma = HEP(4e-33, 'm2', doc = 'Electroweak cross section (typical neutrino cross section) [VALUE MAY BE ERRONEOUS]'+source(5))

#limit on dark matter cross section [5]
darkmatter_sigma_max = HEP(1e-42, 'm2', doc = 'Experimental bound on dark matter cross section [MISSING CONTEXT FOR PROPER DEFINITION]'+source(5))

# =============================================================================
# COSMOLOGY
# =============================================================================

#Current Hubble parameter IN KILOMETER PER SECOND PER MEGAPARSEC [6]
H0 = HEP(67.37, 'km1s_1Mpc_1', doc = 'Hubble constant (Hubble parameter at present time)'+source(6))
smallh = H0.km1s_1Mpc_1 / 100
#Current cosmological parameter : dark energy [6]
Omega0_lambda = 0.685
#Current cosmological parameter : matter [6]
Omega0_m = 0.315
#Current cosmological parameter : radiation [6]
#Omega0_r = 0.0
#Current cosmological parameter : curvature [6]
Omega0_k = 0.0007
#Temperature of the Cosmic Microwave Background [5]
T_CMB = Cosmo(1., 'zplusone', doc = 'Temperature of the Cosmic Microwave Background (CMB)'+source(5))
#Temperature of Neutrino Microwave Background [5]
T_nuCMB = Cosmo(T_CMB.K * ((4/11)**(1/3)), 'K', doc = 'Temperature of the Neutrino Microwave Background'+source(5))
#Redshift at reionization [6]
z_reioniz = Cosmo(1+7.68, 'zplusone', doc = 'Redshift at reionization'+source(6))
#Redshift at recombination [6]
z_recomb = Cosmo(1+1089.8, 'zplusone', doc = 'Redshift at recombination'+source(6))
#Redshift at radiation-matter equality [6]
z_eq = Cosmo(1+3387., 'zplusone', doc = 'Redshift at radiation-matter equality'+source(6))

#age of the Universe [5]
t_universe = HEP(13.8, 'Gyr', doc = 'Age of the Universe'+source(5))
#size of observable Universe [5]
R_universe = HEP(93, 'Glyr', doc = 'Size of observable Universe (i.e. Size of a Hubble bubble) [MISSING CONTEXT : SIZE IS RADIUS OR DIAMETER ?]'+source(5))
#time of recombination [5]
t_recomb = HEP(380000, 'yr', doc = 'Cosmic time at recombination'+source(5))
#time of Lambda-dominated era [5]
"""A VERIFIER"""
t_lambda = HEP(10., 'Gyr', doc = 'Starting time of Lambda-dominated era [VALUE MAY BE ERRONEOUS]'+source(5))

#mean density of photons in the Universe [5]
n_gamma = HEP(411e6, 'm3', doc = 'Mean number density of photons in the Universe'+source(5))
#mean density of neutrinos in the Universe (per flavour) [5]
n_nu = HEP(112e6, 'm3', doc = 'Mean number density of neutrinos in the Universe (per flavor)'+source(5))
#mean mass density of protons in Earth [5]
rho_proton_earth = HEP(1., 'g1cm_3', doc = 'Mean mass density of protons in the Earth'+source(5))
#mean density of protons around the Sun (or in the Solar System) [5]
n_proton_solar = HEP(1e6, 'm3', doc = 'Mean number density of protons in the Solar System'+source(5))
#mean density of protons in the Milky Way [5]
n_proton_milkyway = HEP(1e3, 'm3', doc = 'Mean number density of protons in the Milky Way'+source(5))
#mean density of protons between galaxies [5]
n_proton_galaxy = HEP(1., 'm3', doc = 'Mean number density of protons between galaxies'+source(5))
#mean density of baryons in the Universe (mostly protons) [5]
n_baryon_universe = HEP(0.24, 'm3', doc = 'Mean number density of baryons in the Universe'+source(5))
#mean energy density of dark matter around the Sun (like in the Solar System) [5]
rho_darkmatter_solar = HEP(0.3 * GeV.g, 'g1cm_3', doc = 'Mean mass density of Dark Matter in the Solar System'+source(5))
#mean energy density of dark matter in the Universe [5]
rho_darkmatter_universe = HEP(1e-6 * GeV.g, 'g1cm_3', doc = 'Mean mass density of Dark Matter in the Universe'+source(5))

#number of protons in the visible Universe [5]
number_protons_universe = HEP(R_universe.m**3 * n_baryon_universe.m3, 'one', doc = 'Number of protons in the visible Universe'+source(5))

#mass matter universe [5]
M_universe_matter = HEP(number_protons_universe.one, 'Mp', doc = 'Mass of matter in the visible Universe'+source(5))

#Current effective number of entropic degrees of freedom [14]
g_S_now = 3.91

# =============================================================================
# ELECTROMAGNETIC FIELDS
# =============================================================================

#vacuum impedance [7]
z_vacuum = HEP(120*const.pi, 'ohm', doc = 'Vacuum impedance'+source(7))
#conductivity of copper [7]
gamma_copper = HEP(6e7 / (Ohm.one * m.m), 'm', doc = 'Electric conductivity of copper'+source(7))
#density of electrons in a good conductor [7]
n_conductor = HEP(1e28, 'm3', doc = 'Number density of electrons in a good conductor'+source(7))
#density of electrons in a semiconductor [7]
n_semiconductor = HEP(1e22, 'm3', doc = 'Number density of electrons in a semiconductor'+source(7))
#typical free path time of an electron in a conductor [7]
tau_conductor = HEP(1e-14, 's', doc = 'Typical free path time of an electron in a conductor'+source(7))
#typical speed of an electron in a conductor [7]
electron_speed = HEP(1., 'mmpers', doc = 'Typical speed of an electron in a conductor'+source(7))
#electric field during an electric arc in dry air [7]
E_arc = HEP(3e6, 'V1m_1', doc = 'Electric field during an electric arc in dry air'+source(7))
#number of loops per unit length in a coil [7]
n_loops_coil = HEP(1e4, 'm', doc = 'Typical number of loops per unit length in a coil'+source(7))

#Horizontal component of Earth magnetic field [7]
B_earth_h = HEP(0.2, 'G', doc = 'Horizontal component of Earth magnetic field'+source(7))
#Vertical component of Earth magnetic field [7]
B_earth_v = HEP(0.4, 'G', doc = 'Vertical component of Earth magnetic field'+source(7))
#Magnetic field of a standard magnet [7]
B_magnet = HEP(1e2, 'G', doc = 'Magnetic field of a standard magnet'+source(7))
#Magnetic field of RMI [7]
B_RMI = HEP(1., 'T', doc = 'Magnetic field inside a RMI'+source(7))
#Magnetic field of a superconductive coil [7]
B_supercoil = HEP(10., 'T', doc = 'Magnetic field in a superconductive coil'+source(7))
#Magnetic field of a tokamak [7]
B_tokamak = HEP(100., 'T', doc = 'Magnetic field in a tokamak'+source(7))

#Energy stored in a usual capacitor [7]
E_capacitor = HEP(50e-6, 'J', doc = 'Energy stored in a usual capacitor'+source(7))
#Energy stored in a usual coil [7]
E_coil = HEP(50e-3, 'J', doc = 'Energy stored in a usual coil'+source(7))

#Copper skin depth for various signals [7]
#in usual AC plugs (f = 50Hz)
skin_depth_plug = HEP(1e-2, 'm', doc = 'Skin depth of copper at AC plugs frequency (50Hz)'+source(7))
#in laboratory (f = kHz)
"""CONTEXT ???"""
skin_depth_lab = HEP(1e-3, 'm', doc = 'Skin depth of copper at lab frequency (kHz) [MISSING CONTEXT FOR PROPER DEFINITION]'+source(7))
#for wifi (f = 2.45 GHz)
skin_depth_wifi = HEP(1e-6, 'm', doc = 'Skin depth of copper at wifi frequency (2.45GHz)'+source(7))

#Ionosphere critical frequency [7]
f_ionosphere = HEP(10., 'MHz', doc = 'Ionosphere critical frequency'+source(7))
#Ion number density in the ionosphere [7]
n_ionosphere = HEP(1e12, 'm3', doc = 'Number density of ions in the ionosphere'+source(7))

#Typical section of a conductor wire [7]
section_wire = HEP(1e-6, 'm2', doc = 'Typical section of a conductor wire'+source(7))

#Current density in a usual conductor [7]
j_conductor = HEP(1e6, 'A1m_2', doc = 'Current density in a usual conductor'+source(7))

# =============================================================================
# OPTICS
# =============================================================================

#Radio wavelengths [7]
l_radio_max = HEP(1., 'km', doc = 'Upper bound on radio wavelengths'+source(7))
l_radio_min = HEP(1., 'm', doc = 'Lower bound on radio wavelengths'+source(7))
#Microwave wavelengths [7]
l_micro_max = HEP(0.1, 'm', doc = 'Upper bound on microwave wavelengths'+source(7))
l_micro_min = HEP(1., 'mm', doc = 'Lower bound on microwave wavelengths'+source(7))
#Infrared wavelengths [7]
l_infrared_max = HEP(1., 'mm', doc = 'Upper bound on infrared wavelengths'+source(7))
l_infrared_min = HEP(700., 'nm', doc = 'Lower bound on infrared wavelengths'+source(7))
#Red wavelengths [7]
l_red_max = HEP(700., 'nm', doc = 'Upper bound on red wavelengths'+source(7))
l_red_min = HEP(600., 'nm', doc = 'Lower bound on red wavelengths'+source(7))
#Yellow wavelengths [7]
l_yellow_max = HEP(600., 'nm', doc = 'Upper bound on yellow wavelengths'+source(7))
l_yellow_min = HEP(550., 'nm', doc = 'Lower bound on yellow wavelengths'+source(7))
#Green wavelengths [7]
l_green_max = HEP(550., 'nm', doc = 'Upper bound on green wavelengths'+source(7))
l_green_min = HEP(500., 'nm', doc = 'Lower bound on green wavelengths'+source(7))
#Blue wavelengths [7]
l_blue_max = HEP(500., 'nm', doc = 'Upper bound on blue wavelengths'+source(7))
l_blue_min = HEP(450., 'nm', doc = 'Lower bound on blue wavelengths'+source(7))
#Purple wavelengths [7]
l_purple_max = HEP(450., 'nm', doc = 'Upper bound on purple wavelengths'+source(7))
l_purple_min = HEP(400., 'nm', doc = 'Lower bound on purple wavelengths'+source(7))
#Ultraviolet wavelengths [7]
l_ultraviolet_max = HEP(400., 'nm', doc = 'Upper bound on ultraviolet wavelengths'+source(7))
l_ultraviolet_min = HEP(10., 'nm', doc = 'Lower bound on ultraviolet wavelengths'+source(7))
#X-ray wavelengths [7]
l_xray_max = HEP(10., 'nm', doc = 'Upper bound on x-ray wavelengths'+source(7))
l_xray_min = HEP(10e-3, 'nm', doc = 'Lower bound on x-ray wavelengths'+source(7))
#Gamma ray wavelengths [7]
l_gammaray_max = HEP(10e-3, 'nm', doc = 'Upper bound on gamma-ray wavelengths'+source(7))
l_gammaray_min = HEP(10e-6, 'nm', doc = 'Lower bound on gamma-ray wavelengths'+source(7))
#Sodium doublet wavelength 1 [7]
l_sodiumdoublet1 = HEP(589.0, 'nm', doc = 'First wavelength of sodium doublet (lowest one)'+source(7))
l_sodiumdoublet2 = HEP(589.6, 'nm', doc = 'Second wavelength of sodium doublet (highest one)'+source(7))

#Wifi frequency [7]
f_wifi = HEP(2.45, 'GHz', doc = 'wifi frequency'+source(7))
#frequency at AC power plug [7]
f_plug = HEP(50., 'Hz', doc = 'frequency of AC current in a power plug'+source(7))

#photon energy of visible light [5]
E_visible = HEP(10., 'eV', doc = 'Energy of a photon in visible light domain'+source(5))
#photon energy of X-rays [5]
E_xray = HEP(1., 'keV', doc = 'Energy of a photon in x-ray domain'+source(5)) 
#photon energy of gamma rays [5]
E_gammaray = HEP(1., 'MeV', doc = 'Energy of a photon in gamma-ray domain'+source(5))

#refraction index of vacuum [0]
nvacuum = 1.
#refraction index of air [8]
nair = 1.000293
#refraction index of liquid helium [8]
nliquidhelium = 1.025
#refraction index of water [8]
nwater = 1.333
#refraction index of glass [8]
nglass = 1.458
#refraction index of diamond [8]
ndiamond = 2.417

#laser aperture angle [7]
laser_angle = HEP(1e-3, 'one', doc = 'Aperture angle of a typical laser'+source(7))

#response time for light detection [7]
#of human eye
tau_eye = HEP(0.1, 's', doc = 'Response time for light detection of a human eye'+source(7))
#of photoresistance
tau_photoresistance = HEP(1., 'ms', doc = 'Response time for light detection of a photoresistance'+source(7))
#of photodiode
tau_photodiode = HEP(1., 'mus', doc = 'Response time for light detection of a photodiode'+source(7))
#of photomultiplicator
tau_photomultiplicator = HEP(1., 'ns', doc = 'Response time for light detection of a photo multiplicator'+source(7))


#Solar light power in Southern France [7]
Pi_sun = HEP(1e3, 'W1m_2', doc = 'Power of solar light at Southern France latitude'+source(7))
#Red laser power [7]
Pi_laser = HEP(1e3, 'W1m_2', doc = 'Power of a red laser'+source(7))

#Coherence times [7]
#of natural light 
coherence_sunlight = HEP(1e-15, 's', doc = 'Coherence time of natural light'+source(7))
#of sodium spectral lamp
coherence_Na_lamp = HEP(1e-11, 's', doc = 'Coherence time of a sodium spectral light'+source(7))
#of HeNe laser
coherence_HeNe_laser = HEP(1e-9, 's', doc = 'Coherence time of a helium-neon laser'+source(7))
#of monomode laser
coherence_monomode_laser = HEP(1e-5, 's', doc = 'Coherence time of a monomode laser'+source(7))

#Bandwidth of usual color filters [7]
bandwidth_colorfilter = HEP(10, 'nm', doc = 'Wavelength bandwidth of usual color filters'+source(7))

# =============================================================================
# THERMODYNAMICS
# =============================================================================

#millimeter of mercury [13]
mmHg = HEP(133.322387415, 'Pa', doc = 'One millimeter of mercury'+source(13))

#Air
#molar mass [7]
M_air = HEP(29, 'gpermol', doc = 'Air molar mass'+source(7))
#volumic mass (under normal conditions) [7]
rho_air = HEP(1.33e-3, 'g1cm_3', doc = 'Air volumic mass under normal atmospheric conditions'+source(7))
#speed of sound at 20°C [11]
c_sound_air = HEP(343, 'mpers', doc = 'Speed of sound in air at 20 celsius degrees'+source(11))
#Isopressure thermal expansion coefficient of air (normal conditions) [7]
alpha_air = HEP(1e-2, 'K', doc = 'Isopressure thermal expansion coefficient of air under normal atmospheric conditions'+source(7))
"""This isn't available at the moment since it's expressed in inverse pascal,
which have dimension -4, and i didn't create a -4 dim units dictionary because
for god's sake it's pointless
#Isothermal compressibility coefficient of air (normal conditions) [7]
chi_air = HEP(1e-5, 'Pa', doc = 'Isothermal compressibility coefficient of air under normal atmospheric conditions'+source(7))
"""
#mean quadratic speed of air molecules [7]
v_air = HEP(500, 'mpers', doc = 'Typical mean quadratic speed of air molecules under normal atmospheric conditions'+source(7))

#Water
#molar mass [7]
M_water = HEP(18, 'gpermol', doc = 'Water molar mass'+source(7))
#volumic mass (under normal conditions) [7]
rho_water = HEP(1., 'g1cm_3', doc = 'Liquid water volumic mass under normal pressure conditions'+source(7))
#speed of sound at 20°C [11]
c_sound_water = HEP(1481, 'mpers', doc = 'Speed of sound in liquid water at 20 celsius degrees')
#massic thermic capacity of water [7]
C_water = HEP(4.18e3, 'J1K_1kg_1', doc = 'Massic thermic capacity of liquid water'+source(11))
#massic thermic capacity of ice [7]
C_ice = HEP(2.1e3, 'J1K_1kg_1', doc = 'Massic thermic capacity of ice'+source(7))


#Atmospheric pressure [7]
atm = SI(101325., 'Pa', doc = 'One atmospheric pressure (standard unit)'+source(7))
#Lunar atmospheric pressure [7]
moon_atm = SI(1e-9, 'Pa', doc = 'Moon atmospheric pressure'+source(7))
#Pressure in a vacuum tube [7]
P_vacuum_tube = SI(1e-6, 'Pa', doc = 'Pressure inside a vacuum tube'+source(7))
#Pressure in an incandescent light bulb [7]
P_lightbulb = SI(10., 'Pa', doc = 'Pressure inside an incandescent lightbulb'+source(7))
#Pressure of CO2 in a champagne bottle [7]
P_champagne = SI(5., 'bar', doc = 'carbon dioxyde pressure in a champagne bottle'+source(7))
#Air pressure in a diving cylinder [7]
P_diving = SI(200., 'bar', doc = 'Air pressure inside a diving cylinder'+source(7))

#Mean free path in air at 1 bar pressure [7]
mfp_1bar = HEP(5e-5, 'm', doc = 'Mean free path of air molecules at a pressure of one bar'+source(7))
#Mean free path in air at 1e-8 bar [7]
mfp_10nbar = HEP(50., 'm', doc = 'Mean free path of air molecules at a pressure of 1e-8 bar'+source(7))


#Thermal conductivity at 298 Kelvin [7]
#of copper 
k_copper = HEP(400., 'W1K_1m_1', doc = 'Thermal conductivity of copper'+source(7))
#of typical metal
k_metal = HEP(1e2, 'W1K_1m_1', doc = 'Thermal conductivity of a typical metal'+source(7))
#of concrete
k_concrete = HEP(1., 'W1K_1m_1', doc = 'Thermal conductivity of concrete'+source(7))
#of glass
k_glass = HEP(1., 'W1K_1m_1', doc = 'Thermal conductivity of glass'+source(7))
#of water
k_water = HEP(0.6, 'W1K_1m_1', doc = 'Thermal conductivity of liquid water'+source(7))
#of wood
k_wood = HEP(0.1, 'W1K_1m_1', doc = 'Thermal conductivity of wood'+source(7))
#of air
k_air = HEP(2.5e-2, 'W1K_1m_1', doc = 'Thermal conductivity of air'+source(7))

#Diffusion coefficient [7]
#of air
D_air = HEP(1e-5 * m2.s2, 's', doc = 'Diffusion coefficient of air'+source(7))
#of glass
D_glass = HEP(1e-7 * m2.s2, 's', doc = 'Diffusion coefficient of glass'+source(7))
#of typical metal
D_metal = HEP(1e-5 * m2.s2, 's', doc = 'Diffusion coefficient of a typical metal'+source(7))

#Coefficient conducto-convectif (coefficient de convection thermique) [7]
#of solid/gas
h_solid_gas = HEP(10., 'W1K_1m_2', doc = 'Coefficient conducto-convectif pour un interface solide/gaz [ENGLISH TRANSLATION REQUIRED]'+source(7))
#of solid/liquid
h_solid_liquid = HEP(100., 'W1K_1m_2', doc = 'Coefficient conducto-convectif pour un interface solide/liquide [ENGLISH TRANSLATION REQUIRED]'+source(7))


# =============================================================================
# SOLID MECHANICS
# =============================================================================

#Coefficient of dynamic friction steel vs steel [7]
fd_steel_steel = 0.2
#Coefficient of dynamic friction rubber vs asphalt [7]
fd_rubber_asphalt = 0.6

#Young modulus of iron [7]
E_iron = HEP(190e9, 'Pa', doc = 'Young modulus of iron'+source(7))
#Young modulus of steel [7]
E_steel = HEP(210e9, 'Pa', doc = 'Young modulus of steel'+source(7))
#Young modulus of carbon nanotube cable [7]
E_nanotube = HEP(1000e9, 'Pa', doc = 'Young modulus of carbon nanotube cable'+source(7))

#Steel volumic mass [7]
rho_steel = HEP(7.87e6, 'g1m_3', doc = 'Volumic mass of steel'+source(7))

# =============================================================================
# CHEMISTRY
# =============================================================================

#Oxydoreduction constant  R*T*ln(10)/F at T=298 Kelvin [7]
oxred_constant = HEP(const.R * 298 * np.log(10) / faraday.C, 'V', doc = 'Oxydoreduction constant R*T*ln(10)/F at temperature T=298 Kelvin with R ideal gas constant, ln natural logarithm and F faraday constant'+source(7))

#Standard electrode potential of O2 / H2O [7]
E0_O2_H2O = HEP(1.23, 'V', doc = 'Standard electrode potential of O2 / H2O couple'+source(7))
#Standard electrode potential of Fe3+ / Fe2+ [7]
E0_Fe3_Fe2 = HEP(0.77, 'V', doc = 'Standard electrode potential of Fe3+ / Fe2+ couple'+source(7))
#Standard electrode potential of Fe2+ / Fe [7]
E0_Fe2_Fe = HEP(-0.44, 'V', doc = 'Standard electrode potential of Fe2+ / Fe couple'+source(7))
#Standard electrode potential of saturated calomel electrode [7]
E0_calomel = HEP(0.246, 'V', doc = 'Standard electrode potential of saturated calomel electrode'+source(7))

#Surtension cathodique [7]
#on a platinum electrode
eta_platinum = HEP(-0.01, 'V', doc = 'Surtension cathodique sur electrode de platine [ENGLISH TRANSLATION REQUIRED]'+source(7))
#on a iron electrode
eta_iron = HEP(-0.4, 'V', doc = 'Surtension cathodique sur electrode de fer [ENGLISH TRANSLATION REQUIRED]'+source(7))
#on a zinc electrode
eta_zinc = HEP(-0.75, 'V', doc = 'Surtension cathodique sur electrode de zinc [ENGLISH TRANSLATION REQUIRED]'+source(7))

#Standard enthalpy of usual reactions [7]
DeltaHr0 = HEP(1e5 * J.g, 'gpermol', doc = 'Standard enthalpy of usual chemical reactions'+source(7))
#Standard entropy of usual reactions [7]
DeltaSr0 = HEP(1e2, 'J1K_1mol_1', doc = 'Standard entropy of usual chemical reactions'+source(7))

# =============================================================================
# MISCELLANEOUS
# =============================================================================

#human eye angular resolution [7]
human_eye_resolution = 1e-4 #RADIANS

#room temperature [7]
room_temp = HEP(298, 'K', doc = 'Conventional temperature of a lab room'+source(7))

#Size of a small transistor [7]
transistor_size = HEP(1e-8, 'm', doc = 'Size of small transistor'+source(7))
#Size of a virus [7]
virus_size = HEP(1e-7, 'm', doc = 'Size of a virus'+source(7))
#Size of a Compact Disk alveolus [7]
CD_alveolus_size = HEP(1e-6, 'm', doc = 'Size of a Compact Disk (CD) alveolus. If anyone in the future sees this doc, they will probably laugh very hard at how OUTDATED this value is. COMPACT DISK, REALLY ?! That probably sounds as antiquated to them as floppy disks sound to me.'+source(7))
#Size of a red blood cell [7]
redbloodcell_size = HEP(1e-6, 'm', doc = 'Size of a red blood cell'+source(7))
#width of a hair [7]
hair_width = HEP(1e-4, 'm', doc = 'Width of a hair'+source(7))

#Power delivered by a nuclear plant [7]
P_nuclear = HEP(1., 'GW', doc = 'Power delivered by one nuclear plant'+source(7))