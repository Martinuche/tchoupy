"""Prebuilt conversion system dedicated to units of the International System.

Properties:
- Seven fundamental dimensions with seven fundamental units: time (s),
length (m), mass (kg), temperature (K), electrical current (A), 
amount of substance (mol), luminous intensity (cd). 
"""

__authors__ = ("Martin Teuscher",)
__contact__ = ("teuscher.edu@gmail.com",)
__version__ = "1.0"
__date__ = "2024/04"

import numpy as np
import scipy.constants as const

from .. import ConversionSystem


SI = ConversionSystem('SI',
                         's',
                         [],
                         's0',
                         [
                             # Fake unit for pure numbers like const.alpha or const.pi
                             ('one', 1.),
                             ('rad', 1/(2*np.pi))
                         ],
                         'm',
                         [],
                         'g',
                         [],
                         'K',
                         [],
                         'cd',
                         [],
                         'mol',
                         [],
                         'A',
                         [],
                         
                         SI_prefixes = '_ALL',
                         custom_prefixes = [],
                         
                         Hz = 's_1',
                         C = 'A1s1',
                         
                         N = 'kg1m1s_2',
                         J = 'N1m1',
                         W = 'J1s_1',
                         V = 'W1A_1',
                         ohm = (1., 'V1A_1'), 
                         S = (1., 'ohm_1'), 
                         T = (1., 'V1s1m_2'),
                         
                         min = (const.minute, 's'),
                         hr = (const.hour, 's'),
                         day = (24.0, 'hr'),
                         yr = (const.Julian_year, 's'),
                         
                         mpers = 'm1s_1',
                         mperh = 'm1hr_1',
                         
                         Pa = (1., 'N1m_2'),
                         bar = (1e5, 'Pa'),
                         atm = (1013.0, 'hPa'),
                         
                         Me = (const.m_e, 'kg'),
                         Mp = (const.m_p, 'kg'),
                         Mn = (const.m_n, 'kg'),
                         u = (const.u, 'kg'),
                         angstrom = (const.angstrom, 'm'),
                         AU = (const.au, 'm'),

                         gpermol = 'g1mol_1',
                         earthgravity = (const.g, 'm1s_2'),
                         )
# SI.change_repr('g','kg')