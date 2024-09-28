"""Prebuilt conversion system dedicated to cosmology.

Properties: z+1 (redshift + 1) as its fundamental  unit.
- Only one fundamental dimension, with redshift+1 (zplusone) as its primary
unit.
"""

__authors__ = ("Martin Teuscher",)
__contact__ = ("teuscher.edu@gmail.com",)
__version__ = "1.0"
__date__ = "2024/09"

import numpy as np
import scipy.constants as const

from .. import ConversionSystem


######
# Cosmic Microwave Background (CMB) temperature in KELVIN
# IF YOU CHANGE KELVINS TO ANOTHER UNIT HERE, IMPERATIVELY CHANGE IT TOO IN
# THE DEFINITION OF THE "TCMB" UNIT BELOW
TEMPERATURE_CMB = 2.725 # /!\ KELVINS
######


Cosmo = ConversionSystem('Cosmo',
                         'zplusone',
                         [],
                         # assume TEMPERATURE_CMB is in KELVIN
                         K = (1. / TEMPERATURE_CMB, 'zplusone'),
                         eV = (const.e / const.k, 'K')
                         )