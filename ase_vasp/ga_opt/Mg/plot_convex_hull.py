#!/usr/local/bin/python

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import make_interp_spline, interp1d
from scipy.optimize import curve_fit
import ase
from ase.io import read, write

plt.rcParams.update({'font.size': 18, 'font.family': 'Arial', 'figure.figsize':(9,6)})
ion = 'Mg'
data = pd.read_csv(f'{ion}_ga_energy.csv')

fig, ax = plt.subplots()
ax.scatter(data[f'{ion}_count']/32, data['Form Energy'], color = 'grey')


ax.set_ylabel('Formation energy (eV)', fontweight="bold")
ax.set_xlabel(f'x in {ion}$_x$'+'Fe$_2$'+'(CN)$_6$', fontweight="bold")
ax.set_ylim([-2.5,0])
ax.tick_params(axis='x', labelsize=16)  # Change x-axis tick label size
ax.tick_params(axis='y', labelsize=16)  # Change y-axis tick label size


plt.savefig(f'convexhull_{ion}.pdf')
plt.show()
