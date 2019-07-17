# Header template to use ASE.

import os, sys, re
import numpy as np

import ase
from ase import Atoms
from ase.db import connect
from ase.io import read
from ase.optimize.pcobfgs import PCOBFGS
from kinbot import readers.py

label = '{label}'
atom = {atom}
geom = {geom}

mol = Atoms(symbols = atom, positions = geom)

ppn = {ppn}
method = '{method}'
basis = '{basis}'
mult = {mult}
charge = {charge}
chk = '{chk}'
guess = '{guess}'
integral = '{integral}'
dummy = {dummy}
sella = {sella}
order = {order}
freq = {freq}
task = '{task}'
irc_maxpoints = {irc_maxpoints}
irc_stepsize = {irc_stepsize}
qc = '{qc}'
fix = {fix}
change = {change}
maxattempt = {maxattempt}
qc_command = '{qc_command}'

