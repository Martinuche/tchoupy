"""Prebuilt conversion dictionary dedicated to High Energy Physics.

Properties:
- Only one fundamental dimension, energy, with electronvolt (eV) as its primary
unit.
- Conversions of other units to powers of electronvolt are based on:
    - c = 1 (speed of light.)
    - hbar = 1 (reduced Planck constant)
    - k_B = 1 (Boltzmann constant)
    - q_p = 1 (Planck charge. The Planck charge is defined as 
               q_p = sqrt(4 pi epsilon_0 hbar c). It is chosen so that the
               electrostatic repulsion of two objects with Planck charge and 
               Planck mass balances the Newtonian attraction between them.
               There is big caveat: while values of hbar and c are exact in the
               International System of Units, epsilon_0 is measured experimentally.
               Thus, the link between Planck units and SI units for dimensions
               involving charge or current contains an uncertainty. This is the
               same for the gravitational constant, whose uncertainty is 
               incidentally much bigger than the one on epsilon_0: the link 
               between e.g. the Planck length and a meter is burdened by
               an uncertainty.
               Also, note that q_p = 1 <==> mu_0 = 4 pi in this system of units.
               However, according to the CODATA 2018, in the SI system
               mu_0 = (1 + 5.5e-10)*4 pi*1e-7 N/A^2. Thus, a perhaps better
               choice of normalization would be q_p = 1+5.5e-10 rather than
               1. As of this version, normalization to 1 is preferred.)
"""

__authors__ = ("Martin Teuscher",)
__contact__ = ("teuscher.edu@gmail.com",)
__version__ = "1.0"
__date__ = "2024/04"

import numpy as np
import scipy.constants as const

from .. import ConversionSystem


######
# Sun mass in KILOGRAMS (source: https://en.wikipedia.org/wiki/Earth or
# https://web.archive.org/web/20131110215339/http://asa.usno.navy.mil/static/files/2014/Astronomical_Constants_2014.pdf)
# IF YOU CHANGE KILOGRAMS TO ANOTHER UNIT HERE, IMPERATIVELY CHANGE IT TOO IN
# THE DEFINITION OF THE "Msun" UNIT BELOW
SUN_MASS = 1.9884e30 # /!\ KILOGRAMS
######

HEP = ConversionSystem('HEP',
                'eV',
                [
                  ('J', 1. / const.e),
                  
                  ('g', (1./const.kilo) * const.c**2 / const.e),
                  # Or: const.value['Planck mass energy equivalent in eV']
                  ('MPl', np.sqrt(const.c**5 * const.hbar / const.G) / const.e),
                  
                  ('K', const.k / const.e),
                  
                  ('V', const.e),
                ],
                
                'eV0',
                [
                 ('Planck', 1.),
                 # Fake unit for pure numbers like const.alpha or const.pi
                 ('one', 1.),
        
                 ('mol', const.N_A),
                 
                 # Definition of the Planck charge: see beginning of the file
                 ('C', np.sqrt(const.alpha) / const.e),
                ],
                
                'eV_1',
                [
                  ('s', const.e / const.hbar),
                 
                  ('m', 1. / const.c * const.e / const.hbar),
                 
                  ('tPl', (np.sqrt(const.hbar * const.G / const.c**5) 
                           * const.e / const.hbar)),
                ],
                
                SI_prefixes = ['T', 'G', 'M', 'k', 'c', 'm', 'mu', 'n'],
                custom_prefixes = [],
                
                ### Composite units ###
                
                erg = (const.erg, 'J'),
                 
                Msun = (SUN_MASS, 'kg'), # assume SUN_MASS is in KILOGRAMS
                reducedMPl = (1./np.sqrt(8.*const.pi), 'MPl'),
                Me = (const.m_e, 'kg'),
                Mp = (const.m_p, 'kg'),
                Mn = (const.m_n, 'kg'),
                u = (const.u, 'kg'),
                 
                Hz = 's_1',
                 
                # Definition of the Planck charge: see beginning of the file
                A = (1., 'C1s_1'),
                 
                gpermol = 'g1mol_1',
                 
                earthgravity = (const.g, 'm1s_2'),
                 
                lightspeed = 'one',
                mpers = 'm1s_1',
                mperh = (1. / const.hour, 'm1s_1'),
                
                ohm = (1., 'V1A_1'),   
                
                lPl = 'tPl',
                min = (const.minute, 's'),
                hr = (const.hour, 's'),
                yr = (const.Julian_year, 's'),
                 
                pc = (const.parsec , 'm'),
                lyr = (const.c * const.Julian_year, 'm'), # = const.light_year
                AU = (const.au, 'm'),
                angstrom = ( const.angstrom, 'm'),
                 
                N = (1., 'kg1m1s_2'),
                T = (1., 'V1s1m_2'),
                G = (1e-4, 'T'),
                W = (1., 'J1s_1'),
                 
                barn = (1e-28, 'm2'),
                 
                Jy = (1e-26, 'W1m_2Hz_1'),
                 
                Pa = (1., 'N1m_2')
                 )