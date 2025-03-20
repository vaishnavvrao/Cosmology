# Look through and make sure you have all these packages installed before running! 
import sys
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import ascii,fits
import astropy.units as u
import astropy.constants as const
from astropy.wcs import WCS
from astropy.table import Table
from astropy.coordinates import SkyCoord,match_coordinates_sky,Angle
from astropy.coordinates import match_coordinates_sky
from scipy.interpolate import UnivariateSpline
from scipy.io import readsav
import matplotlib.patches as patches
import matplotlib.path as path
from scipy.spatial import cKDTree
from scipy import stats
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as PathEffects
from astroquery.simbad import Simbad
import seaborn as sns #For KDE plot
from astropy.visualization import wcsaxes
from sklearn.neighbors import KernelDensity #for KDE
from scipy.stats import poisson

#For making postage stamps
from astropy.visualization import make_lupton_rgb
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clipped_stats
from astropy.visualization.wcsaxes import add_scalebar
from astropy.visualization.wcsaxes import SphericalCircle

# These are my personal preferences for the default aesthetics of Matplotlib figures
# You can change these to whatever you want! 
params = {'figure.figsize': (8,8),
          'text.usetex': True,
          'font.family': 'serif',
          'font.size': 20,
          'xtick.direction': 'in',
          'xtick.top': True,
          'xtick.labelsize': 16,
          'ytick.direction': 'in',
          'ytick.right': True,
          'ytick.labelsize': 16,
          'axes.titlesize': 20,
          'axes.labelsize': 20,
          'legend.fontsize': 18,
          'axes.spines.left': True,
          'axes.spines.bottom': True,
          'axes.spines.top': True,
          'axes.spines.right': True,
          'axes.grid.axis': 'both',
          'axes.grid.which': 'both'}

pub_params = {'figure.figsize': (7.2,7.2),
          'text.usetex': True,
          'font.family': 'serif',
          'font.size': 10,
          'xtick.direction': 'in',
          'xtick.top': True,
          'xtick.labelsize': 8,
          'ytick.direction': 'in',
          'ytick.right': True,
          'ytick.labelsize': 8,
          'axes.titlesize': 10,
          'axes.labelsize': 10,
          'legend.fontsize': 9,
          'axes.spines.left': True,
          'axes.spines.bottom': True,
          'axes.spines.top': True,
          'axes.spines.right': True,
          'axes.grid.axis': 'both',
          'axes.grid.which': 'both'}
mpl.rcParams.update(params)


# Quiet all warnings...because they add clutter and are usually unhelpful
import warnings
warnings.filterwarnings('ignore')