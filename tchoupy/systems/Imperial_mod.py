"""Prebuilt conversion system dedicated to the metric and imperial systems.

Properties:
- Two fundamental dimensions, length and time, whose respective primary units
are meter (m) and second (s).
"""

__authors__ = ("Martin Teuscher",)
__contact__ = ("teuscher.edu@gmail.com",)
__version__ = "1.0"
__date__ = "2024/09"

import numpy as np
import scipy.constants as const

from .. import ConversionSystem

Imperial = ConversionSystem('Imperial',
                         'm',
                         [],
                         's',
                         [],
                         ft = (const.foot, 'm'),
                         mile = (const.mile, 'm'),
                         hr = (const.hour, 's'),
                         mileperh = 'mile1hr_1',
                         )