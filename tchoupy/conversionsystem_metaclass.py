"""Module containing a metaclass used to generate unit systems for
conversions and dimensional quantity computations."""

__authors__ = ("Martin Teuscher",)
__contact__ = ("teuscher.edu@gmail.com",)
__version__ = "1.0"
__date__ = "2024/04"


from fractions import Fraction
from functools import reduce
import numpy as np
# regex was used as the re module does not support multiple lookbehinds, see
# https://stackoverflow.com/questions/58310184/python-regex-error-look-behind-requires-fixed-width-pattern
import regex
import warnings

# Other module from the package
# https://stackoverflow.com/questions/8953844/import-module-from-subfolder
from . import qty


class ConversionSystem(type):
    """Class factory (metaclass) for systems of unit conversions.

    Take any number of unit names and conversion rules between some of these and
    return a class. Any instance of this class will then be a physical quantity
    that can be expressed within this system of units and easily converted
    into any valid unit.
    Moreover, each generated class overrides most of arithmetic operation
    builtins (e.g. __add__) so that instances can easily be manipulated as
    numbers with +-*/** (although with an associated physical dimension).

    See README.md file for more information on following attributes and methods.
   
    Metaclass attributes
    --------------------
    _SI_prefixes : dict
        list of the international system of units' prefixes (dict keys,
        string type) and their corresponding values (dict values, float type) 
     
    Class (instance of metaclass) attributes
    ----------------------------------------
    prefixes : list
    MAX_DENOMINATOR : int
    nb_total_dims : int
    fund_dims : list of string
    _conversions : dict
    _common_error_msg : string
        
    Instance (instance of class) attributes
    ---------------------------------------
    inverse_
        @property method returning 1/instance
    
    Metaclass methods
    -----------------
    
    Class methods
    -------------
    __getitem__(key)
        return information about the unit key, if available within this class
    __setitem__(key, val)
        add unit key to the set of units available for this class
    change_repr(old_unit, new_unit)
        replace a primary unit of the system by another unit
    _add_prefixes
        internal method
    
    Instance methods
    ----------------
    __getattr__
    __repr__
    __add__
    __radd__
    __sub__
    __rsub__
    __mul__
    __rmul__
    __truediv__
    __rtruediv__
    __pow__
    """

    
    _SI_prefixes = {  'Y':1e24,
                      'Z':1e21,
                      'E':1e18,
                      'P':1e15,
                      'T':1e12,
                      'G':1e9,
                      'M':1e6,
                      'k':1e3,
                      'h':1e2,
                      'da':1e1,
                      '':1.,
                      'd':1e-1,
                      'c':1e-2,
                      'm':1e-3,
                      'µ':1e-6,
                      'mu':1e-6,
                      'n':1e-9,
                      'p':1e-12,
                      'f':1e-15,
                      'a':1e-18,
                      'z':1e-21,
                      'y':1e-24,
                   }
    
    def __new__(metacls, clsname, 
                *units_and_subunits,
                SI_prefixes = [],
                custom_prefixes = [],
                MAX_DENOMINATOR = qty.MAX_DENOMINATOR,
                **composite_units):
        """Create an instance of the ConversionSystem metaclass.
    
        Create a new instance, set part of its namespace with the prefixes (SI
        and non-SI) that will later be used to name units, then pass it to the
        __init__ constructor.
        
        Parameters
        ----------
        
        See README.md file for more information.
       
        metacls : type type
            refer to ConversionSystem class
        clsname : string
            name of the class to be created by the class factory. Good practice
            recommends it is the same as the name of the variable it is assigned
            to, e.g. HighEnergy = ConversionSystem('HighEnergy', ...).
            It can be different, but there is no need to complicate things.
        *units_and_subunits : even number of arguments; even-numbered ones
                              are strings and odd-numbered ones are lists
        SI_prefixes : iterable of strings, optional
            Sublist of International System of Units prefixes that the user 
            wants to use to name units. They must be non-empty keys of
            ConversionSystem._SI_prefixes.
            Default is [].
            If set to '_ALL', all SI prefixes will be integrated to the class.
        custom_prefixes : iterable of 2-tuple, optional
            Prefixes that the user wants to use to name units, but that are
            not part of the SI prefixes. Each tuple must be of the form
            (string giving the prefix' name, float giving the prefix' value).
            Default is [].
        MAX_DENOMINATOR : int, optional
            Maximal irreducible denominator tolerated by the fractions module.
            See the qty module for further information.
            Default is qty.MAX_DENOMINATOR.
        **composite_units : any number of keyword arguments of physics_string type
            Additional unit names to be registered in the class 

        Raises
        ------
        ValueError
            One of the prefixes is not made of alphabetical characters.

        Returns
        -------
        ConversionSystem object, which is a class itself
        """
        
        if SI_prefixes == '_ALL':
            SI_prefixes = list(metacls._SI_prefixes.keys())
        if not np.all(np.array(list(map(
                lambda string : string.isalpha(),
                map(lambda tup : tup[0], custom_prefixes)
                )))):
            raise ValueError('Prefixes must be alphabetical.')
        
        prefixes_dict = {
            # prefix '' ensures the unit names will also be present without prefix
            '': 1.,
            # select desired prefixes from built-in International System prefixes
            **dict(map(lambda key : (key, metacls._SI_prefixes[key]), SI_prefixes)),
            # add user-defined prefixes
            **dict(custom_prefixes),
           }
        
        # super() refers to 'type', which is subclassed.
        # Return a new class, instance of the ConversionSystem metaclass,
        # and generate part of its namespace, mainly bound methods that will be
        # used to perform basic operations +-*/ on instances of this new class. 
        return super().__new__(metacls, clsname, (), {
            'prefixes' : prefixes_dict,
            'MAX_DENOMINATOR' : MAX_DENOMINATOR,
            '__init__' : qty._init__,
            '__getattr__' : qty._getattr__,
            '__repr__' : qty._repr__,
            '__add__' : qty._add__,
            '__radd__' : qty._radd__,
            '__sub__' : qty._sub__,
            '__rsub__' : qty._rsub__,
            '__mul__' : qty._mul__,
            '__rmul__' : qty._rmul__,
            '__truediv__' : qty._truediv__,
            '__rtruediv__' : qty._rtruediv__,
            '__pow__' : qty._pow__,
            'inverse_' : qty.inverse_
            })
        
        # Note for my own understanding :
        # https://stackoverflow.com/questions/114214/class-method-differences-in-python-bound-unbound-and-static    
        # object.method(*) is translated to class.method(object, *).
        #So there is no memory waste for bound methods, there are redefined each
        # time an instance is created, although object.method is class.method
        # will return False. object.method is still calling class.method, but
        # fixes its first argument.
        # On the contrary, there is no such translation for attributes, so 
        # `object.attr is class.attr` will return True if attr was defined in
        # the class body.
    
    # customary "self" has been replaced by "cls" for the sake of clarity,
    # since instances of ConversionSystem are classes 
    # (this is codewise genuine)
    def __init__(cls, clsname, 
                 *units_and_subunits,
                 SI_prefixes = [],
                 custom_prefixes = [],
                 MAX_DENOMINATOR = qty.MAX_DENOMINATOR,
                 **composite_units):
        """
        Instantiate the ConversionSystem object returned by the __new__ method.
        
        See README.md file for information about the arguments.
        
        Parameters
        ----------
        clsname : string
            name of the class to be created by the class factory. Good practice
            recommends it is the same as the name of the variable it is assigned
            to, e.g. HighEnergy = ConversionSystem('HighEnergy', ...).
            It can be different, but there is no need to complicate things.
        *units_and_subunits : even number of arguments; even-numbered ones
                              are strings and odd-numbered ones are lists
        SI_prefixes : iterable of strings, optional
            Sublist of International System of Units prefixes that the user 
            wants to use to name units. They must be non-empty keys of
            ConversionSystem._SI_prefixes.
            Default is [].
            If set to '_ALL', all SI prefixes will be integrated to the class.
        custom_prefixes : iterable of 2-tuple, optional
            Prefixes that the user wants to use to name units, but that are
            not part of the SI prefixes. Each tuple must be of the form
            (string giving the prefix' name, float giving the prefix' value).
            Default is [].
        MAX_DENOMINATOR : int, optional
            Maximal irreducible denominator tolerated by the fractions module.
            See the qty module for further information.
            Default is qty.MAX_DENOMINATOR.
        **composite_units : any number of keyword arguments of physics_string type
            Additional unit names to be registered in the class 

        Raises
        ------
        ValueError
            Some unit name contains non-alphabetical characters or consists
            of a single 'o'. A single 'o' is protected because it is used to
            denote a fraction bar in attribute names. For instance if x is an
            instance of cls representing 1 sqrt(Hz), x.Hz1o2 calls for the
            value of x in sqrt(Hz) (/ is not a valid character for an
            attribute name, and _ is already used to denote a negative power,
            e.g. (x**-2).Hz_2 ).
        KeyError
            A compound unit refered to a unit that had not been defined yet.
            Recall that keyword arguments can be ordered in a unique way in
            python, as dictionaries can as long as their size is not modified.
            https://stackoverflow.com/questions/58227855/best-way-performance-wise-to-iterate-over-a-dictionary

        Returns
        -------
        None

        """
        
        super().__init__(clsname, (), {})
        
        cls._conversions = {}
        cls.fund_dims = []
        
        # Find base units written before digits in the unit names given by the
        # user. A unit that contains no number will necessarily be added to this set.
        indep_units = reduce(
            lambda unitset, newunit : (
                    unitset | set(regex.split('_?\d+', newunit, 1)[:1]) 
                                       ), 
            units_and_subunits[::2],
            set()
            )
        
        cls.nb_total_dims = len(indep_units)
        dim_ids = {}
        counter = 0
        
        try:       
            for (unit, secondary_units) in (
                    zip(units_and_subunits[::2], units_and_subunits[1::2])):
                  
                base_unit, power_num, power_den = qty.split(unit)
                
                if base_unit == 'o':
                    raise ValueError('One of your units is: \"o\". This name is'
                                     'restricted as a keyword to denote a ' 
                                     'fraction, so it cannot be used to name'
                                     'a unit. Underscores cannot be contained'
                                     'in the name of a unit either.') 
                    
                if power_den > cls.MAX_DENOMINATOR:
                    warnings.warn('An attempt to use rational powers with large' 
                                  ' denominator was detected. Please '
                                  'consider changing the value of MAX_DENOMINATOR '
                                  'to an appropriate value, otherwise some '
                                  'conversions may be crudely approximated by the '
                                  'program.', category = UserWarning)
                
                #If the key exists in the dict, return the value, else return
                #the specified default object
                unit_id = dim_ids.get(base_unit, 
                                      np.full(cls.nb_total_dims, Fraction(0,1))
                                      )
                
                if np.all(unit_id == 0):
                    
                    # I do this rather than unit_id[counter]=1 to preserve the
                    # Fraction dtype
                    unit_id[counter] += 1
                    counter += 1
                    dim_ids[base_unit] = unit_id
                    # cls.fund_dims[i] contains the fundamental (primary) unit of
                    # dimension (0,...,0,1,0...0) with 1 at position i
                    # (starting by 0)
                    cls.fund_dims.append(base_unit)
                    
                    cls._add_prefixes(base_unit, unit_id, 1.)
                        
                unit_id = Fraction(power_num, power_den) * unit_id         
            
                # concatenate the lists is slow but append would modify the input
                for (secondary_unit, conversion_value) in secondary_units:
                    cls._add_prefixes(secondary_unit, unit_id, conversion_value)
            
            for (unit, definition) in composite_units.items():
                # either the kwarg is e.g. g = (9.81, 'm1s_2')
                if isinstance(definition, tuple):
                    #Summing up this to cls[unit] = definition is much less readable
                    conversion_value, sequence_of_units = definition
                    cls[unit] = (conversion_value, sequence_of_units)
                # or, if the constant is 1., it is defined through shortcut,
                # e.g. N = 'kg1m1s_2'
                elif isinstance(definition, str):
                    cls[unit] = (1., definition)
                else:
                    raise ValueError('Composite units must be defined either as '
                                     'a string of units or a 2-tuple of (value, ' 
                                     'string of units). This is not the case of '
                                     f'{unit}.')
                    
        except KeyError as e:
            raise KeyError(f'Unit {e.args[0]} referenced before assignment.')
    
        # https://stackoverflow.com/questions/2058925/how-can-i-break-up-this-long-line-in-python
        # Start of an error message which is then completed depending on the
        # situation. See qty module.
        cls._common_error_msg = ('The space of dimensional numbers is a graded'
                                  ' algebra: you can add two numbers of the '
                                  'same dimension, or multiply a number by a '
                                  'dimensionless scalar or another number, but '
                                  )       
        
        # Unfortunately for some reason what follows has to go OUTSIDE the qty
        # module (probably because it depends on the instance) : magic method
        # bound to class instance are not looked up. I have not understood every
        # thing and I do not understand what that last sentence means, refer to
        # https://stackoverflow.com/questions/24063788/python3-singledispatch-in-class-how-to-dispatch-self-type
        # and
        # https://docs.python.org/3/reference/datamodel.html#special-method-lookup
        
        # To differentiate the following operations :
        # a) scalar (float) * dimensional number
        # b) dimensional number * dimensional number,
        # we use functools.singledispatchmethod:
        # several functions are registered under the same name, then type
        # inference is performed. See qty module.
        @qty._mul__.register(cls)
        def _qty_mul(self, other : cls) -> cls:
            return cls(self._val * other._val,
                       '',
                       from_dimension = self.dim + other.dim
                       )
            """
            Multiply a dimensional number by another.
            
            qty._mul__ overrides __mul__.
            
            Parameters
            ----------
            other : instance of cls (the class-instance of ConversionSystem)
                dimensional number
    
            Returns
            -------
            instance of cls
                dimensional number product of self and other.
    
            """
            
        # To differentiate the following operations :
        # a) scalar (float) * dimensional number
        # b) dimensional number * dimensional number
        # the functionalities of functools.singledispatchmethod is used :
        # several functions are registered under the same name, then type
        # inference is performed. See qty module.
        @qty._truediv__.register(cls)
        def _qty_div(self, other : cls) -> cls:
            """
            Divide a dimensional number by another.
            
            qty._truediv__ overrides __truediv__.
            
            Parameters
            ----------
            other : instance of cls (the class-instance of ConversionSystem)
                dimensional number

            Returns
            -------
            instance of cls
                dimensional number quotient of self by other.

            """
            return cls(self._val / other._val,
                       '',
                       from_dimension = self.dim - other.dim
                       )
    
    def _add_prefixes(cls, alphabetical_unit, unit_id, conversion_value):
        """
        Add to cls all the unit names composed of a single unit preceded by all
        the registered prefixes.

        All these unit names are stored in the dict cls._conversions.
        
        Parameters
        ----------
        alphabetical_unit : string
            unit name (eV, kg, mol, ...).
        unit_id : numpy array
            dimension of alphabetical_unit in the space of units that is proper
            to cls. If this space has k primary units, unit_id[k] gives the 
            dimension in the k-th direction.
        conversion value : float
            conversion factor between alphabetical_unit and the primary unit
            of same dimension, defined as some product of powers of the primary
            units of cls.

        Raises
        ------
        ValueError
            Unit name contains non-alphabetical characters or consists
            of a single 'o'. A single 'o' is protected because it is used to
            denote a fraction bar in attribute names. For instance if x is an
            instance of cls representing 1 sqrt(Hz), x.Hz1o2 calls for the
            value of x in sqrt(Hz) (/ is not a valid character for an
            attribute name, and _ is already used to denote a negative power,
            e.g. (x**-2).Hz_2 ).

        Returns
        -------
        None.

        """
        if not alphabetical_unit.isalpha() or alphabetical_unit == 'o':
            raise ValueError('Unit names must contain only letters and be '
                             'different from a single \"o\". This is not '
                             f'the case of {alphabetical_unit}.')
        for prefix in cls.prefixes:
            prefixed_unit = prefix + alphabetical_unit
            if prefixed_unit in cls._conversions:
                raise ValueError(f'Unit {prefixed_unit} is defined in '
                                'too many ways.')
            #The setitem method cls[prefixed_unit] is not used here
            #because it needs the primary units to be created first.
            #Basically the setitem method is here for the user only, to be
            #used after instance creation.
            cls._conversions[prefixed_unit] = (unit_id, conversion_value * cls.prefixes[prefix])
    
    def __getitem__(cls, key):
        """
        Syntactic sugar to retrieve information on some unit through cls[key].
        
        key refers to the unit the user wants to look up.

        Parameters
        ----------
        key : string
            A valid unit name (can include a prefix), present in cls.

        Raises
        ------
        KeyError
            The desired unit cannot not found in cls.

        Returns
        -------
        tuple, size 2
            The result of cls[unit], providing key is a valid unit in cls, is a 
            2-tuple consisting of
                - an array representing the dimension of key in the space of
                  units that is proper to cls. If this space has k primary units,
                  array[k] gives the dimension in the k-th direction.
                - the conversion factor between key and the primary unit of same
                  dimension, defined as some product of powers of the primary
                  units of cls. 
        """
        
        # try/except VS if time cost: same efficiency when exception doesnt occur 
        # https://stackoverflow.com/questions/2522005/cost-of-exception-handlers-in-python
        try:
            return cls._conversions[key]
        except KeyError as e:
            raise KeyError(f'Unit {e.args[0]} not recognized in'
                           f' {cls.__name__} system.')
    

    def __setitem__(cls, key, val):
        """
        Syntactic sugar to add new units to cls even after instantiation.

        cls[key] = (factor, unit) will add the unit named key to cls such that
        key = factor*unit. unit can be a compound unit here e.g. m1s_1.
        Prefixes will automatically be added in front of key.
        
        Parameters
        ----------
        key : string
            Must only contain alphabetical characters and differ from a single 'o'.
        val : 2-tuple
            See above.

        Returns
        -------
        None.

        """
        pre_conversion_value, composite_unit = val
        new_object = cls(pre_conversion_value, composite_unit)
        cls._add_prefixes(key, new_object.dim, new_object._val)
     
    def change_repr(cls, old_unit, new_unit):
        """
        [TODO] Replace one primary unit of cls by another chosen by the user.
        
        The primarym units of cls are originally inferred at instantiation.
        It is possible to change them one by one, providing the replacement unit
        still allow to span all the space of units (this is checked by the 
        program). The new unit must already exist inside cls.        
        
        This function does not ensure compatibility of conversions before and 
        after the change of primary units. Changing the primary unit
        at instantiation and consider it unmutable should be preferred.
        TODO : assert everything that shall be modified by this function
        indeed is. (typically I initially forgot to modify cls.fund_dims).
        
        Parameters
        ----------
        old_unit : string
            unit present in cls.fund_dims.
        new_unit : string
            unit present in cls._conversions.

        Raises
        ------
        ValueError
            Requirements described above are not met.

        Returns
        -------
        None.

        """
        
        dim_index = np.where(np.array(cls.fund_dims) == old_unit)[0]
        
        if dim_index.size == 0:
            raise ValueError('This function is made to replace a primary '
                             f'unit of this system of conversions. {old_unit} '
                             'is not one of them.')
        dim_index = dim_index[0]
        try:
            dim_vector, conversion = cls[new_unit]
        except KeyError:
            raise ValueError(f'{new_unit} must be made of letters only.')
        pivot = dim_vector[dim_index]
        if pivot == 0:
            raise ValueError(f'The change from {old_unit} to {new_unit} does '
                             'not allow to cover all the dimensional space. The'
                             ' vectors are not free.')
        if pivot > cls.MAX_DENOMINATOR:
           warnings.warn('An attempt to use rational powers with large '
                         'denominator was detected. Please '
                         'consider changing the value of MAX_DENOMINATOR '
                         'to an appropriate value, otherwise some '
                         'conversions may be crudely approximated by the '
                         'program.', category = UserWarning) 
        for (unit, (dim_array, value)) in cls._conversions.items():
            #/!\/!\ Parentheses are very important to conserve a Fraction type !!
            # vector * (dim_array[dim_index]/pivot) would turn everything into floats !!
            new_dim_array = dim_array - (dim_vector * dim_array[dim_index]) / pivot
            new_dim_array[dim_index] = dim_array[dim_index] / pivot
            new_val = value * conversion**(- dim_array[dim_index] / pivot)
            
            cls._conversions[unit] = (new_dim_array, new_val)
        
        cls.fund_dims[dim_index] = new_unit
        
        warnings.warn('Be careful not to combine any two quantities whose one '
                      'is defined and the other after the call to change_repr. '
                      'Values of quantities defined before are left unchanged '
                      'and so are meaningless after a call to change_repr.',
                      category = UserWarning)

