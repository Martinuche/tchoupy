"""Toolbox building methods for manipulating physical quantities and
convert them in arbitrary units. This is not to be used by itself but inside
the classes defined in conversionsystem_metaclass."""

__authors__ = ("Martin Teuscher",)
__contact__ = ("teuscher.edu@gmail.com",)
__version__ = "1.0"
__date__ = "2024/04"


from fractions import Fraction
from functools import reduce, singledispatchmethod
import numpy as np
# regex was used as the re module does not support multiple lookbehinds, see
# https://stackoverflow.com/questions/58310184/python-regex-error-look-behind-requires-fixed-width-pattern
import regex

# See __pow__ method for the use case below. Basically, if an object representing
# a dimensional quantity is taken to a float power, this power will be
# approximated by the closest rational with denominator less than or equal to
# MAX_DENOMINATOR. This ensures that the dimensions of physical quantities being
# manipulated (e.g. 2 in cm^2) remain exact (i.e. rational) and do not suffer
# from floating point errors. This is also why the fractions module is used.
# See README.md file for more information.
MAX_DENOMINATOR = 100

def split(unit):
    """
    Split a rational power of a unit into its name, numerator and denominator.

    Parameters
    ----------
    unit : string
        Assumed to match the following pattern :
        - starts by a nonzero number of alphabetical characters, except for
          a single 'o'
        - optionally followed by (
             - an optional underscore (if present, it denotes a negative sign
               preceding the numerator)
             - a nonzero number of digit characters : the numerator of the power
             - optionally followed by (
                 - the letter 'o', denoting a fraction bar (short for 'over')
                 - a nonzero number of digit characters : the denominator
              )
          )
    Said differently, accepted patterns are e.g. cm, cm2, cm_2, cm1o2, cm_1o2,
    which stand for centimeter, (cm)^2, (cm)^(-2), (cm)^(1/2), (cm)^(-1/2).
    Notice that the prefix centi is included in the power as naturally expected.
    
    Returns
    -------
    3-tuple of string, int, int
        first element is the name of the unit (e.g. cm), second is the numerator
        (1 if not specified), third is the numerator (1 if not specified)

    """
    
    split = regex.split('(?<=[A-Za-z]{2}|(?!o)[A-Za-z])(?=\d)|_|(?<=\d)o', unit)
    return (
            #name of the base unit
            split[0],
            #integer numerator
            1 if len(split) < 2 else int(split[1]) * (-1 if '_' in unit else 1),
            #integer denominator
            1 if len(split) < 3 else int(split[2])
           )

def split_composite(composite_unit):
    """
    Split a succession of power of units into names and digits, using split().

    Parameters
    ----------
    composite_unit : string
        Assumed to match a nonzero repetition of the pattern detailed in the
        split() method docstring.

    Returns
    -------
    list of tuples of string, int, int
        list containing the outputs of the split() method applied to each
        pattern present in the composite_unit string.

    """
    
    # ?: is here to avoid mess, but I'm unsure of what it does, see:
    # https://docs.python.org/3/library/re.html#re.findall
    individual_units = regex.findall('[A-Za-z]+(?:_?\d+(?:o\d+)?)?', composite_unit)
    
    return list(map(split, individual_units))


def _init__(self, value, unit,
             from_dimension = None,
             doc = ''):
    """
    Instantiate an object representing a physical number.
    
    A physical number consists of a value (float) and a unit, associated to a
    certain physical dimension.
    
    The class of this object is generated by the ConversionSystem factory
    in the conversionsystem_metaclass module. When this class is created, its 
    __init__ method is set equal to the hereby _init__ function.

    Parameters
    ----------
    value : float
        value of the physical quantity when expressed in the unit 'unit'.
    unit : string
        Must match one of the following patterns, same as in split_composite()
        1) 'abc' or any sequence of at least one alphabetic character,
           except for a single 'o'
        2) 'abc52' or any 1)-case followed by any sequence of at least
           one digit
        3) 'abc_52' or any 1)-case followed by '_' (denoting a minus sign)
           followed by any sequence of at least one digit
        4) 'abc52o33' or any 2)-case followed by 'o' (denoting a fraction bar)
           followed by any sequence of at least one digit
        5) 'abc_52o33' or any 3)-case followed by 'o' (denoting a fraction bar)
           followed by any sequence of at least one digit
        6) a repeted sequence of any of the above cases, e.g.
           'abc52de_3f1g_3o5hi_1'.
        Physically speaking this corresponds to the product of units
            (abc)^52 * (de)^(-3) * f * g^(-3/5) * (hi)^(-1)
        where abc, de, f, g, hi would be valid unit names like kg, mol, mm...
    from_dimension : None or numpy array, optional
        Used for internal calls only when not None. Default is None.
    doc : string, optional
        Describe the meaning of the physical object (e.g. 'Sun Mass').
        Default is ''.

    Raises
    ------
    ValueError
        The unit does not match the pattern mentioned above.

    Returns
    -------
    None.

    """
    cls = self.__class__
    
    # Internal use case: the object is already expressed in irreducible units
    # (This is to gain time when multiplying or adding objects of the class)
    if from_dimension is not None:
        self._val = value
        self.dim = from_dimension
    else:
        
        # Check whether unit matches the pattern explained in the docstring
        # flags=regex.IGNORECASE is not used because a single "o" must be
        # refused while a single capital "O" is valid 
        if not bool(regex.match('^([A-Za-z]+(_?\d+(o\d+)?)?)+$', unit)):
            raise ValueError('Misspelled unit.')
        else:
            unitlist = split_composite(unit)
        
            # Commpute the dimension of unit and conversion value of value*unit
            # into the (possibly compound) primary unit of same dimension.
            # tup is an element of unitlist, of the form (name, numer, denom)
            # cls[unit_name] is a tuple ( dimension of unit_name , conversion
            # of unit_name to (possibly compound) primary unit of same dimension ).
            
            self.dim = reduce(lambda array, tup : (
        array + cls[tup[0]][0] * Fraction(tup[1], tup[2]).limit_denominator(
                                    max_denominator = cls.MAX_DENOMINATOR)
                                        ),
                              unitlist,
                              np.full(cls.nb_total_dims, Fraction(0,1))
                             )
            
            self._val = reduce(lambda val, tup : val * (cls[tup[0]][1] ** (
                                                            tup[1]/tup[2])), 
                               unitlist,
                               value )
            
            self.__doc__ = doc
    
    # Compute string representing the primary unit of the instance,
    # based on the primary units of its class
    self._unit_repr = ''
    for (fraction, primary_unit) in zip(self.dim, cls.fund_dims):
        if fraction == 0:
            pass
        else:
            # Numerator of power, even if equal to 1, is always displayed in
            # order to clearly separate two units. Else, a confusion could occur
            # between e.g. "millikelvin" (mK) and "meter times kelvin" 
            # (here displayed as m1K1 and not mK).
            self._unit_repr += (
                primary_unit
                + ('_' if fraction.numerator < 0 else '')
                + f'{abs(fraction.numerator)}'
                + ('' if fraction.denominator == 1 else f'o{fraction.denominator}')
                               )
    
def _getattr__(self, unit):
    """
    Return the value of self (a physical quantity) expressed in the unit unit.

    Parameters
    ----------
    unit : string
        Must match the pattern described in split_composite() or __init__.

    Raises
    ------
    ValueError
        unit does not match the pattern.

    Returns
    -------
    float
        Value of self after conversion in unit "unit".
    """
    
    # flags=regex.IGNORECASE is not used because a single "o" must be
    # refused while a single capital "O" is valid 
    if not bool(regex.match('^([A-Za-z]+(_?\d+(o\d+)?)?)+$', unit)):
        if unit.startswith('__') or unit.endswith('__'):
            #TODO : This happens for internal calls when the class is called
            #with weird builtins, like __abstractmethods__ or when you try to
            #use it with numpy and it calls __array_struct__
            #This is probably not good use to pass, but it's just that I want
            #the triggered error to be consistent with what the user is trying
            #to do, not to raise "misspelled unit" when they're actually asking
            #for e.g. np.sqrt(x) where x is an instance of self
            #It still fails, but the error displayed is more consistent.
            pass
        else:
            raise ValueError('Misspelled unit.')
    else:    
        cls = self.__class__
        
        unitlist = split_composite(unit)
        
        # c.f. _init__ method above where the same computation is performed
        asked_dim = reduce(lambda array, tup : array + (
                                    cls[tup[0]][0] * Fraction(tup[1], tup[2]) ),
                           unitlist,
                           np.full(cls.nb_total_dims, Fraction(0,1))
                           )
        if np.any(self.dim != asked_dim):
            raise ValueError(f'Quantity {self} cannot be expressed in [{unit}]'
                             f' within {cls.__name__} system. Dimensions '
                             'do not match.')
        
        # c.f. _init__ method above where the same computation is performed
        return self._val / reduce(lambda val, tup : val * (cls[tup[0]][1] ** (
                                                                tup[1]/tup[2])), 
                           unitlist,
                           1. )
            
# https://www.pythontutorial.net/python-oop/python-__repr__/
def _repr__(self):
    """
    Display a physical quantity as a float value followed by a unit.
    
    The unit is always a product of powers of the primary units chosen
    in the object's class.
    """
    
    return str(self._val) + '[' + self._unit_repr + ']'


#=============================
# All the following are meant ot override basic operations + - * / **, allowing
# physical quantities to be manipulated as numeric types.
# https://docs.python.org/3/reference/datamodel.html?highlight=rmul#emulating-numeric-types
#=============================

def _add__(self, other):
    if not np.all(self.dim == other.dim):
        raise ValueError(self.__class__._common_error_msg + 
        'adding two numbers of different dimensions is not allowed.')
    return self.__class__(self._val + other._val,
                          '',
                          from_dimension = self.dim,)
_radd__ = _add__

def _sub__(self, other):
    if not np.all(self.dim == other.dim):
        raise ValueError(self.__class__._common_error_msg + 
        'substracting two numbers of different dimensions is not allowed.')
    return self.__class__(self._val - other._val,
                          '',
                          from_dimension = self.dim,)

def _rsub__(self, other):
    if not np.all(self.dim == other.dim):
        raise ValueError(self.__class__._common_error_msg + 
        'substracting two numbers of different dimensions is not allowed.')
    return self.__class__(other._val - self._val,
                          '',
                          from_dimension = self.dim,)
 
# See conversionystem_metaclass for why this functionality is being used
# https://stackoverflow.com/questions/24063788/python3-singledispatch-in-class-how-to-dispatch-self-type 
@singledispatchmethod
def _mul__(self, other):
    return NotImplemented

# The other registered function is implemented in conversionsystem_metaclass
@_mul__.register(int)
@_mul__.register(float)
def _scalar_mul(self, scalar):
    return self.__class__(self._val * scalar,
               '',
               from_dimension = self.dim
               )

_rmul__ = _mul__

# See conversionsystem_metaclass for why this functionality is being used
# https://stackoverflow.com/questions/24063788/python3-singledispatch-in-class-how-to-dispatch-self-type
@singledispatchmethod
def _truediv__(self, other):
    return NotImplemented

# The other registered function is implemented in conversionsystem_metaclass
@_truediv__.register(int)
@_truediv__.register(float)
def _scalar_div(self, scalar):
    return self.__class__(self._val / scalar,
               '',
               from_dimension = self.dim
               )

def _rtruediv__(self, scalar):
    return self.__class__(scalar / self._val,
               '',
               from_dimension = -self.dim
               )

def _pow__(self, scalar):
    #Fraction module automatically does gcd and sign simplification
    frac = Fraction(scalar).limit_denominator(
                            max_denominator = self.__class__.MAX_DENOMINATOR)
    
    return self.__class__(self._val ** frac,
               '',
               from_dimension = self.dim * frac
               )

@property
def inverse_(self):
    """
    Alternative way to write 1/self as self.inverse_.

    Returns
    -------
    same type as self
        1/self

    """
    return self.__class__(1./self._val,
               '',
               from_dimension = -self.dim
               )


