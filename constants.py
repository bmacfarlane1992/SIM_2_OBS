'''

constants.py

Definitions of units, and conversion factors.

Author: Benjamin MacFarlane
Contact: bmacfarlane@uclan.ac.uk
Date: 27/02/2018

'''


# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #
# --------------- General constants/conversions --------------- #
# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #


### - - - - - Size of axis/title text and ticks in plots - - - - - ###

fontsize = 18

### - - - - - Size of legend text and ticks in plots - - - - - ###

leg_fontsize = 12

# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #
# --------------- General constants/conversions --------------- #
# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #


### - - - - - Gravitational constant (SI) - - - - - ###

G_grav = 6.674e-11

### - - - - - Speed of light (SI) - - - - - ###

c = 2.998e8

### - - - - - Speed of light (CGS) - - - - - ###

c_cgs = 2.998e10

### - - - - - Planck constant (SI) - - - - - ###

h = 6.626e-34

### - - - - - Proton mass (SI) - - - - - ###

m_p = 1.673e-27

### - - - - - Stefan-Boltzmann constant (SI) - - - - - ###

k_b = 1.381e-23

### - - - - - Mean molecular density - - - - - ###

mu = 2.35

### - - - - - Atomic mass unit (CGS) - - - - - ###

kg_per_amu = 1.661e-27

### - - - - - Atomic mass unit (CGS) - - - - - ###

g_per_amu = 1.661e-30

### - - - - - g -> kg conversion (SI) - - - - - ###

kg_per_g = 1.000e-3

### - - - - - g -> kg conversion (SI) - - - - - ###

g_per_kg = 1.00 / kg_per_g

### - - - - - H2 mass (SI) - - - - - ###

kg_per_H2 = 3.35e-27

### - - - - - H2 mass (cgs) - - - - - ###

g_per_H2 = kg_per_H2 * g_per_kg

### - - - - - m s^{-1} -> cm s^{-1} conversion (CGS) - - - - - ###

cms_per_ms = 1.000e2

### - - - - - km s^{-1} -> cm s^{-1} conversion (CGS) - - - - - ###

cms_per_kms = 1.000e5

### - - - - - m s^{-1} -> km s^{-1} conversion - - - - - ###

kms_per_ms = 1.000e-3

### - - - - - km s^{-1} -> m s^{-1} conversion (SI) - - - - - ###

ms_per_kms = 1.0 / kms_per_ms

### - - - - - Stefan-Boltzmann constant (SI) - - - - - ###

sigma_sb = 5.670e-8

### - - - - - Stefan-Boltzmann constant (CGS) - - - - - ###

sigma_sb_cgs = 5.670e-5

### - - - - - m -> cm conversion (CGS) - - - - - ###

cm_per_m = 1.000e2

### - - - - - cm -> m conversion (CGS) - - - - - ###

m_per_cm = 1.00 / cm_per_m

### - - - - - m^2 -> cm^2 conversion (CGS) - - - - - ###

cm2_per_m2 = cm_per_m**2.00

### - - - - - m^2 -> cm^2 conversion (CGS) - - - - - ###

m2_per_cm2 = 1.00 / cm2_per_m2

### - - - - - Micron -> m conversion (SI) - - - - - ###

m_per_micron = 1.000e-6

### - - - - - Micron -> cm conversion - - - - - ###

cm_per_micron = 1.000e-4

# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #
# --------------- Astronomical constants/conversions --------------- #
# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #


### - - - - - Solar radius (SI) - - - - - ###

m_per_rsol = 6.955e+8

### - - - - - Solar radius (CGS) - - - - - ###

cm_per_rsol = 6.955e+10

### - - - - - Solar radius (CGS) - - - - - ###

rsol_per_cm = 1.0 / cm_per_rsol

### - - - - - Assumed dust-to-gas ratio (dimensionless) - - - - - ###

dust_to_gas = 0.010

### - - - - - Assumed dust-to-gas ratio (dimensionless) - - - - - ###

gas_to_dust = 1.0 / dust_to_gas

### - - - - - Solar mass (SI) - - - - - ###

kg_per_msol = 1.998e+30

### - - - - - Solar mass (CGS) - - - - - ###

g_per_msol = 1.998e+33

### - - - - - Solar mass (CGS) - - - - - ###

msol_per_g = 1.0 / g_per_msol

### - - - - - Astronomical units (AU - SI) - - - - - ###

m_per_au = 1.496e11

### - - - - - Astronomical units (AU - SI) - - - - - ###

au_per_m = 1.0 / m_per_au

### - - - - - Astronomical units (CGS) - - - - - ###

cm_per_au = 1.496e13

### - - - - - cm -> AU conversion - - - - - ###

au_per_cm = 1.0 / cm_per_au

### - - - - - Parsec -> AU conversion - - - - - ###

au_per_pc = 206265.000

### - - - - - Parsec -> cm conversion - - - - - ###

cm_per_pc = au_per_pc * cm_per_au

### - - - - - Parsec -> cm conversion - - - - - ###

m_per_pc = cm_per_pc * m_per_cm

### - - - - - Parsec -> cm conversion - - - - - ###

rsol_per_pc = cm_per_pc * rsol_per_cm

### - - - - - Jansky (SI) - - - - - ###

Jy_SI = 1.000e-26

### - - - - - Jansky (cgs) - - - - - ###

Jy_cgs = 1.000e-23

### - - - - - erg (SI) - - - - - ###

erg_per_J = 1.000e7

### - - - - - erg (SI) - - - - - ###

J_per_erg = 1.00 / (erg_per_J)

### - - - - - Solar Luminosity (SI) - - - - - ###

Lsol_SI = 3.828e26

### - - - - - Solar Luminosity (cgs) - - - - - ###

Lsol_cgs = Lsol_SI * erg_per_J
