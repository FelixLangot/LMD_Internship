import adsfunc as ads
import scipy.stats as sp
from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib import font_manager as fm, rcParams
import matplotlib.pyplot as plt
import os
from numba import njit, jit
import numpy as np
from tqdm.notebook import tqdm as tqdmn  # progress bar
import jax.numpy as jnp

import seaborn as sns
sns.set(color_codes=True)
sns.set_palette(sns.color_palette("Blues_r", 3)+['#D3D3D3', '#808080'])

mpl.rcParams['axes.linewidth'] = 0.2  # set the value globally
mpl.rcParams['xtick.major.width'] = 0.3
mpl.rcParams['xtick.minor.width'] = 0.3
mpl.rcParams['ytick.major.width'] = 0.3
mpl.rcParams['ytick.minor.width'] = 0.3


fpath = os.path.join(
    rcParams["datapath"], "/Users/felixlangot/Library/Fonts/Cronos-Pro-Light_12448.ttf")
Csprop = fm.FontProperties(fname=fpath)
Cstitleprop = Csprop.copy()
Cslabelprop = Csprop.copy()
Cstitleprop.set_size(30)
Csprop.set_size(15)
Cslabelprop.set_size(25)

fpath2 = os.path.join(
    rcParams["datapath"], "/Users/felixlangot/Library/Fonts/MinionPro-Regular.otf")
Mpprop = fm.FontProperties(fname=fpath2)
Mptitleprop = Mpprop.copy()
Mplabelprop = Mpprop.copy()
Mptitleprop.set_size(30)
Mpprop.set_size(25)
Mplabelprop.set_size(15)

fpath3 = os.path.join(
    rcParams["datapath"], "/Users/felixlangot/Downloads/mnsymbol/tex/MnSymbol.sty")
Tgprop = fm.FontProperties(fname=fpath3)
Tgtitleprop = Tgprop.copy()
Tglabelprop = Tgprop.copy()
Tgtitleprop.set_size(30)
Tgprop.set_size(25)
Tglabelprop.set_size(15)
