# tchoupy
# Toolbox for Conversions (H ?) and Operations and Units with PYthon
# Time-saving tool for unit conversions and dimensional quantity computations

# ===INSTALLATION===

The project is written in **Python3** and makes use of various scientific
libraries which you have to install. To do so, we recommend using pip.
A good practice in python is to set up a virtual environment
```bash
python -m venv my_env
```
Now before any utilisation, please activate the virtual environment
```bash
source my_env/bin/activate
```
Install the requirements using pip.
Note: if pip is unable to find some librairies, it may be that they are already
builtins of your python installation. This often happens with e.g. the warnings
package.
```bash
pip install -r requirements.txt
```

# ===QUICKSTART===

The core of this package allows you to

1) create variables representing physical quantities: x = SI(3, 'km') represents
a length of 3 kilometers ('SI' is explained in the next section)

2) easily convert these variables into any compatible unit: x.m equals 3000
(m stands for meter), x.cm equals 100000, x.mile equals 1.864, etc

3) manipulate these quantities as numbers : x*x gives another physical quantity
equal to 9 (km)², 1/x**(1/2) can be expressed in meter^(-1/2) etc

4) perform more advanced functionalities described in what follows

## Different systems for unit conversions

A crucial thing to understand is that two units may or may not be convertible 
into one another, depending on what conversions you are allowing. For example,
In the International System of units (SI), 1 second cannot be expressed in
meters, while in High Energy Physics (HEP) it can, because the speed of light
is set to 1, implying that 1 second = 299792458 meters. You could also be a
meteorologist for whom height of water on the ground (in millimeters) is
equivalent to a volume of rainwater (in liters), thus wish to convert liters 
into millimeters and vice-versa; but obviously this conversion would not make
any sense outside of your field.

To account for this, you need to specify which set of rules ('system') you are 
using for conversions  when creating an object representing a physical quantity. 
Python-wise, this object will be an instance of the conversion system you want 
to use (which itself is a class, hence a callable). The syntax is very easy. 
Just call the system's name with a value and a unit as its two arguments:
```python
from tchoupy.systems import HEP, SI
>>> time1 = HEP(60, 's')   #60 seconds within High Energy Physics system
>>> time2 = SI(60, 's')   #60 seconds within International System of units
```
The value of the created object in a given unit is then retrieved by putting
the unit (its abbreviatied form) as an attribute of the object:
```python
>>> time1.min   # min stands for minute
1.0
>>> time2.min
1.0
>>> time1.m   # m stands for meter
17987547480.0
>>> time2.m
ValueError: Quantity 60.0[s1] cannot be expressed in [m] within SI system.
Dimensions do not match.
```
The latter error signals an incompability between seconds and meters in SI.


To provide as much versatility as possible, several systems are already
implemented in the package. As of version 1.0 you have access to :

- SI : conversion system according to the International System of units.
There are seven base units that are independent of each other (meaning, they
cannot be converted into one another) : second (s), meter (m), gram* (g),
kelvin (K), ampere (A), candela (cd), mole (mol). 
Additional units are: Hz, C, N, J, W, V, ohm, S, T, min, hr, day, yr, mpers,
mperh, Pa, bar, atm, Me, Mp, Mn, u, angstrom, AU, gpermol, earthgravity, 
optionally preceded by any of the following prefixes:
Y, Z, E, P, T, G, M, k, h, da, d, c, m, µ, mu, n, p, f, a, z, y. 

- HEP : conversion system for High Energy Physics, based on c = k_B = hbar = 1
and the Planck charge q_p = sqrt(4 pi epsilon_0 hbar c) = 1.
There is only one base unit : electronvolt (eV). Therefore, all units have a
dimension equal to some power of eV and any two units with the same dimension
can be converted into one another.
Additional units are J, g, MPl, K, V, Planck, one, mol, C, s, m, tPl, erg, 
Msun, reducedMpl, Me, Mp, Mn, u, Hz, A, gpermol, earthgravity, lightspeed,
mpers, mperh, ohm, lPl, min, hr, yr, pc, lyr, AU, angstrom, N, T, G, W, barn,
Jy, Pa, optionally preceded by any of the following prefixes: 
T, G, M, k, c, m, mu, n.

- Cosmo : TODO.

To import any of them, just look in the 'systems' submodule:
```python
>>> from tchoupy.systems import SI
```

## How to write a valid unit name

So far we have only dealt with quantities whose units are straightforward to
express (min, s, km). But what if the quantity's unit would be made of several
'subunits', e.g. meter per second, joule per kelvin, cubic meter and so forth?

Let us break this down with a few examples, then explain the rationale behind:

```python
>>> speed = SI(36, 'm1s_1')   # 36 meter per second
>>> speed.km1h_1   # x in kilometer per hour
10.0
>>> bar = SI(1e5, 'kg1s_2m_1')   # 1 bar = 1e5 pascal = 1e5 kg/m/s^2
>>> bar.Pa   # 1 bar in pascals
100000.0
>>> bar.kg1s_2m_1
```
So, to write a compound unit, concatenate the subunits, separated by the number
at which each subunit is powered. Note that whenever this number is negative,
an underscore (_) is used in place of a minus sign (-). The rationale is
visible in the examples: as an unit name shall always be a valid python
attribute name, it cannot contain anything else than letters, numbers or _.
A power of 1 cannot be omitted, because the program would wrongly decipher the
string: 'ms_1' would be read millisecond^(-1) and not meter * (second^(-1)).

You may also encounter the case when the power is a fractional number, but
in a similar fashion the character '/' cannot be part of an attribute name. You
must use 'o' instead (alias for 'over'), as in
```python
>>> length = SI(1, 'ha1o2')   # a length expressed in hectare^(1/2)
>>> length.m    
100.   # the side of a square with surface 1 hectare is 100 meter long
```
Finally, note that any prefix in the unit name is automatically included inside
the power, as naturally expected in the common language. For instance:
```python
>>> density = SI(1000, 'm_3')   # a number density of 1000 per cubic meter
>>> density.cm_3   # y per cubic centimeter
0.001   # and obviously not 10.
```

With all of the above you can express quantities in rather complex units. A last 
famous example from cosmology is the Hubble constant:
```python
>>> H0 = HEP(1.44e-42, 'GeV')
>>> H0.km1s_1Mpc_1   # H0 in kilometer per second per megaparsec
67.506760860605
```

Extra note: what if you ask for the object itself, without specifying any
argument? The object's builtin method __repr__ has been modified to return 
its value expressed in the primary units of the system, under a string
form:
```python
>>> bar = SI(1e5, 'kg1s_2m_1')   # 1 bar = 1e5 pascal = 1e5 kg/m/s^2
>>> bar
100000.0[s_2m_1kg1]   # 1e5 per squared second per meter times kg kilogram
```

## Operations on quantities

The last powerful feature of tchoupy is to help simplify your computations when
handling with physical quantities, so that you don't have to express all your
variables in the same unit or worry about unit conversion. It allows you
to add, substract, multiply, divide and take the power of the objects you define
(it overrides builtin methods like __add__ or __truediv__), then returns another
object whose unit is consistent with your computation. Consider the following:
```python
>>> length = SI(300,'m')   # 300 meters
>>> length2 = SI(5, 'km')   # 5 kilometers
>>> speed = SI(5, 'm1s_1')   # 5 meter per second
>>> volume = length**3
>>> type(volume)
conversionsystem_metaclass.SI
>>> volume.dam3   # volume in cubic decameter
27000.0
>>> length + length2
5300.0[m1]   # 5300 meters
>>> length / speed
60.0[s1]   # it takes 60 seconds to travel 300m at 5m/s
>>> time = length / speed
>>> time.min   # time in minutes
1.0
>>> length + time
ValueError: The space of dimensional numbers is a graded algebra: you can add
two numbers of the same dimension, or multiply a number by a dimensionless
scalar or another number, but adding two numbers of different dimensions
is not allowed.
```
Note that the last error message would not have been triggered if HEP system
had been chosen over SI system, as time and length have the same dimension
from the HEP point of view.

To maintain consistency, it is only possible to mix quantities defined from
the same system.

Remark: floating powers of physical quantities (e.g. length**3.45) are
approximated by the closest rational number (with sufficiently low denominator),
in order to only keep track of rational dimensions. This is motivated by 
physics considerations, as there is to our knowledge no example of irrational
dimensions (e.g. meter**sqrt(2)). See advanced features for more details.

# ===ADVANCED FEATURES===

## Create one's own conversion system

To increase the versatility of possible conversions, the user can define their
own conversion systems. They can create any unit and any conversion rules that
are suitable for their purpose. Let's see how.

To create a new system of units, you must call the ConversionSystem class (the
returned object will itself be a class, as HEP or SI above, hence
ConversionSystem - CoSy henceforth - is a so-called "metaclass").
Suppose for instance that we are working on a hiking website, and want to
convert between lengths and durations of hiking paths. We could define a
conversion system as follows:
```python
from tchoupy import ConversionSystem
Hikes = ConversionSystem('Hikes', 

                     'hundredmeter',
                     [('minute', 0.7), ('mile', 16.09), ('step', 1/140)],
                     
                     'gram',
                     [('ounce', 28.3)],
                     
                     'hundredmeter2',
                     [('hectare', 1.)],
                     
					 SI_prefixes = ['k', 'c', 'm'],
                     custom_prefixes = [('big', 2.)],
                     MAX_DENOMINATOR = 100,
					 
                     meter = (0.01, 'hundredmeter'),
                     hour = (60., 'minute'),
                     second = (1/60, 'minute'),
                     heartbeat = (0.8, 'second'),
                     mileph = 'mile1hour_1',
                     kilostepph = 'kstep1hour_1',
                     imperialhectare = 'mile2',
                     
                     )
```
1) 1st argument to CoSy is the name of your unit conversion system. In principle
it can be different from the name of the variable to which it is assigned,
but there is no reason to complicate things. So, Hikes is a system as SI
or HEP before.

2) Following the name come (any even number of) positional arguments, working
by pair. The first argument of a pair is a unit. The second is a list of other
units that can be converted into the primary unit you just specified, together
with the conversion ratio. 

2a) Here, we have defined our first fundamental unit, 'hundredmeter'.
We then defined a minute, equal to 0.7 hundredmeter (based on
the assumption that someone walks 70 meters in one minute), a mile 
(1609 meters) and a step (a 140th of 100meter or 71cm). We emphasize that
'meter' or 'cm' do not exist in the conversion system, we simply
mention them for the sake of clarity. Any two of these units can then be 
converted into one another:
```python
>>> short_hike = Hikes(3, 'mile')
>>> short_hike.minute
68.95714285714286   # It takes arounds 69 minutes to walk 3 miles
```

2b) The primary unit of the second pair is 'gram', and we define what an ounce
is (28.3 grams). Each individual primary unit is treated as independent from 
the others, so conversion between hundredmeter and gram are supposed to be
impossible.

2c) The third pair refers to a previously defined primary unit ('hundredmeter')
that has been raised to some power. Contrary to 2b), no new primary unit is 
defined here, but additional unit "hectare" is created. As a consequence,
units equivalent to hundredmeter, when raised to the 
appropriate power, will automatically be compatible with hectare:
```python
>>> other_hike = Hikes(30, 'minute')   # A 30-minute long hike
>>> other_hike.hundredmeter
21.0
>>> (other_hike**2).hectare   # squaring a quantity equivalent to length 
441.0
```
All in all, the number of independant dimensions here is 2 (hundredmeter and
grams).

3) Keyword arguments:

3a) In SI_prefixes, list the common prefixes you want to make available in your
system. They will be attached to unit names. Here, since we included 'k'
(kilo), units khundredmeter, kminute, kstep, kgram, kounce and khectare are
thus available. Same is true for 'c' (0.01) and 'm' (0.001):
```python
>>> afewsteps = Hikes(2000.0, 'step')
>>> afewsteps.kstep
2.0
>>> surface = Hikes(1.0, 'khundredmeter2')
>>> surface.hundredmeter2
1000000.0
```
IMPORTANT: in the second example, note that the power following the unit ALSO
applies to the prefix value (1000^2 = 1e6). This will always be the case:
cm3 (cubic centimeter) = (cm)^3 = 1e-6 (m^3).

IMPORTANT: if you want ALL common prefixes from yotta to yocto to be included,
you can set `SI_prefixes = "_ALL"` instead of having to write the whole list.


3b) In custom_prefixes you can generate new (alphabetical characters only)
prefix abbreviations that are not part of the classical list of prefixes. 
The rest is similar to 3a).

3c) MAX_DENOMINATOR: see below.

4) Extra keywords arguments:
This feature makes it practical to define a lot of extra units, as soon as they
can be defined by a composition of other previously defined units, with a 
optional prefactor. As you can see, the syntax is:
    `new_unit_name = (prefactor, composite_unit_string)`
where a valid string for composite_unit_string value has been defined in the 
Quickstart section above. A few important remarks:
- the syntax `new_unit_name = composite_unit_string` is also valid and is a 
    shortcut for `new_unit_name = (1.0, composite_unit_string)`
- as dictionaries have a canonical order in python, you can use a new unit 
    defined this way inside the composite_unit_string of another, as long as 
    the latter is defined after the former was. See the definition of heartbeat
    or mileph above: they are based on units second and hour, themselves defined
    as a keyword argument.
- new_unit_name must only be made of alphabetical characters
- Prefixes will also be attached to the units defined this way.
 
### MAX_DENOMINATOR

Best use: do not touch this argument.

The possibility of powering physical dimensions (cubic meter, second squared,...)
has already been mentioned. However, to our knowledge, physics always deals 
with powers being rational numbers (for an example of a non-integer power, 
hertz^(1/2) is a unit sometimes used in signaling theory. But more simply, we 
would like (x^3)^(1/3) to return x, so allowing for integer powers immediately 
calls for the need of fractional powers anyway). Mathematically speaking, the 
space of physical dimensions is a Q-vector space (or even, a Q-graded algebra). 

Ensuring compatibility between units and quantities therefore implies to keep 
track of these rational dimensions: if x can be expressed in cm, x**3 must be 
expressible in cm3, but not in cm2. As we chose to represent unit names as
strings, yhis requirement is naturally in conflict with python floating point
errors, since e.g. x**(0.1+0.2) should a priori be expressed in 
`cm0.30000000000000004` !
To bypass this problem, we need to round up floating point errors to determine 
what power was "intented" by the user" (here 0.3 was intended although python 
returned 0.30000000000000004). The rounding precision is then an arbitrary 
parameter, left to the user. The following rule has been retained:
When a physical quantity is raised to some float power, the program will 
replace the float power by the closest rational with denominator (in 
irreducible form) lower than or equal to MAX_DENOMINATOR. This way, the 
dimension remains a rational fraction of the original dimension.

The default is MAX_DENOMINATOR=100, as it is assumed that fractional powers 
with denominator bigger than 100 very rarely occurs in physics.

All of this is achieved thanks to the Fractions module, which is reputably slow.
Future versions might leverage this dependency.

## Add custom units to a preexisting conversion system

You can add a unit to an already existing system this way:
```python
>>> Masses = ConversionSystem('Masses',
                                  'gram',
                                  [('ton', 1e6), ('ounce', 28.3)]
                                  )
>>> Masses['kilo'] = (1000.0, 'gram')
```
The syntax is similar to item 4) in the last subsection:
    `System[new_unit_name] = (prefactor, composite_unit_string)`
and the types & restrictions are the same. Prefixes included in the system 
(see 3a) of the last subsection) will also be attached to this new unit.

This can be very handy when you need a system similar to a builtin one,
except for a few units, but don't want to rebuild it entirely from scratch:
```python
>>> from tchoupy.systems import HEP as HEP_but_with_speed_of_sound
>>> HEP_but_with_speed_of_sound['soundspeed'] = (343.0, 'm1s_1')
>>> HEP_but_with_speed_of_sound(1., 'soundspeed').lightspeed   
1.1441248465296617e-06   # speed of sound in speed of light unit
```

## Change the fundamental unit of a (preexisting) system

TODO. This is already available through the change_repr(old, new) method, but
is not safe as of the first version.

# ===COMMON PITFALLS===

Below we list a few common mistakes that are worth knowing about.

1) Be careful not to use an unappropriate primary units for your system
(e.g. electronvolt while all your computations are about miles and kilometers).
As the program internally converts every value into primary units, very
large conversion factors may quickly appear, and even be enhanced by squaring 
or cubing your units. This can potentially lead to large floating point errors
or even overflows. To avoid this behaviour, create your own system with an
appropriate primary unit, or use the change_repr method of a builtin dict
to change one of its fundamental units (FEATURE TO COME).

2) Because of prefixes, overlaps in unit definitions may happen when defining
your own systems - or adding new units to an existing one. For instance,
 "giga-al" (Gal) and gallons (Gal) cannot be present in the same system. 
This is checked by the program and it will raise an error if that happens.  

