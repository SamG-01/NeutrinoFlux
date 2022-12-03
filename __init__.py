# Initialization

## Standard Scientific Computing Packages
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.constants as C

## SciPy Numerical Methods
from scipy.integrate import quad, nquad
from scipy.interpolate import InterpolatedUnivariateSpline, RectBivariateSpline

## MCEq (for atmo flux)
###import solver related modules
from MCEq.core import MCEqRun
import mceq_config
###import primary model choices
import crflux.models as pm

# Constants
mass_density_ice = 0.917 / (C.centi)**3 # g/m^3

R_E = 6371 * C.kilo # m; radius of the earth

nucleon_mass = 1.67e-27 * C.kilo # gram
nucleons_per_gram_earth = 1/nucleon_mass # used to convert mass density to number density
electrons_per_gram_earth = 0.5 * nucleons_per_gram_earth # number of electrons = number of protons, which make up half the nucleon mass

nucleons_per_gram_water = nucleons_per_gram_earth
electrons_per_gram_water = 10/18 * nucleons_per_gram_water # 10 electrons per 18 nucleons

seconds_per_year = 31556925.2160000 # used to convert from events/second to events/year

# Loads Pickled Fit Functions
import os
import dill as pickle

path = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))

with open(path + "/Fits/cross_sections", "rb") as f:
    cross_sections = pickle.load(f)
with open(path + "/Fits/attenuation_parameter", "rb") as f:
    attenuation_parameter = pickle.load(f)
with open(path + "/Fits/atmo_flux", "rb") as f:
    atmo_flux_funcs = pickle.load(f)
