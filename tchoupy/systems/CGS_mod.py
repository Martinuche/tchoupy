"""Prebuilt conversion system dedicated to units of the CGS System.

Properties:
- Seven fundamental dimensions with seven fundamental units: time (s),
length (cm), mass (g), temperature (K), electrical current (A), 
amount of substance (mol), luminous intensity (cd). 
"""

__authors__ = ("Martin Teuscher",)
__contact__ = ("teuscher.edu@gmail.com",)
__version__ = "1.0"
__date__ = "2024/09"

import numpy as np
import scipy.constants as const

CGS = lambda x : x+2