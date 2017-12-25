# If you just have a stupid sense of humor, you'll never run out of things to laugh at
# "The larger our ignorance, the stronger the magnetic field." - Woltier
"""
Exoplanets Final Project: N-Body Simulation
@author: Jonas Powell
October 2017
"""


## PACKAGES ###
import datetime

import astropy.constants as c
#from twilio.rest import Client


### CONSTANTS ###
from barnes_hut_gridding import node, tree
from brute_force import systemBF

G = -c.G.value                                              # m3 kg-1 s-2
AU = c.au.value                                             # m
mSol = c.M_sun                                              # kg
rSol = c.R_sun                                              # m
mEarth = c.M_earth                                          # kg
rEarth = c.R_earth                                          # m


# This is the only place dt is defined (not passed as an argument anywhere). JSYK

# For stellar scales (v2434)
#dt = 30000
#axlims = 220 * AU

dt = 1800                                                   # duration of one timestep (s) (1800s = half hour)
axlims = 2 * AU                                             # AU
theta = 0.0                                                 # for BH Approx
steps = []                                                  # the output file

# Automatically make good new file names.
now = datetime.datetime.now()
months = ['jan', 'feb', 'march', 'april', 'may', 'june', 'july', 'aug', 'sep', 'oct', 'nov', 'dec']
today = months[now.month - 1] + str(now.day)
logfileName = today + 'run.log'
gifOutName = today + 'run.gif'
trajOutName = 'z' + today + 'run.png'

gifNames = []
elapsedTimeBF = []
durationBF = []
elapsedTimeBH = []
durationBH = []


n = node('root', [-axlims, axlims], [-axlims, axlims], None)
bh = tree(n)
bh.visualizeTree()
bh.run_BarnesHut(30000, 0)

bf = systemBF()
bf.run_BruteForce(100, 0)

# The End
