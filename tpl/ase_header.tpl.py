"""
Header template to use ASE.
"""

import os, sys, re
import numpy as np

import ase
from ase import Atoms
from ase.db import connect
from ase.io import read
from ase.optimize.pcobfgs import PCOBFGS

label = '{label}'
atom = {atom}
geom = {geom}

mol = Atoms(symbols = atom, positions = geom)

