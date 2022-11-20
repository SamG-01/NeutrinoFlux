# Initialization

## Standard Scientific Computing Packages
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.constants as C

## SciPy Numerical Methods
from scipy.optimize import curve_fit
from scipy.integrate import quad, nquad
from scipy.interpolate import InterpolatedUnivariateSpline, RectBivariateSpline

## Uncertainty Handling
from uncertainties import ufloat
import uncertainties.unumpy as unp

## MCEq
###import solver related modules
from MCEq.core import MCEqRun
import mceq_config
###import primary model choices
import crflux.models as pm

# Constants
V_eff = (C.kilo)**3 # m**3; effective detector volume

R_E = 6371 * C.kilo # m; radius of the earth

nucleon_mass = 1.67e-27 * C.kilo # gram
nucleons_per_gram = 1/nucleon_mass # used to convert mass density to number density

seconds_per_year = 31556925.2160000 # used to convert from events/second to events/year
number_density_water = C.N_A / (C.centi)**3 # / m**3; water equivalent for N_A

# Loads Pickled Fits
import os
import dill as pickle

path = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))

with open(path + "/Fits/cross_sections", "rb") as f:
    cross_sections = pickle.load(f)
with open(path + "/Fits/attenuation_function", "rb") as f:
    attenuation_function = pickle.load(f)
with open(path + "/Fits/atmo_flux", "rb") as f:
    flux_funcs = pickle.load(f)
