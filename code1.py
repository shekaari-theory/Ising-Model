#!/usr/bin/python3
# This code is free & open-source, 
# distributed under the terms of GNU General Public 
# License (Version 3, 29 June 2007), 
# https://www.gnu.org/licenses/gpl-3.0.txt .
#   
#   Programmer: Ashkan Shekaari
#
#      Contact: shekaari.theory@gmail.com
#               shekaari@email.kntu.ac.ir
#
# Social Links: https://rt.academia.edu/AshkanShekaari
#               https://orcid.org/0000-0002-7434-467X
#  
# >>> PLEASE, cite the main paper in case of using this code, as follows:
#     A. Shekaari, M. Jafari (2021). "Theory and simulation of the Ising model". arXiv:2105.00841v1 [cond-mat.stat-mech]. doi:10.48550/arXiv.2105.00841.
#
# PURPOSE: Ising model in 2 dimensions.
#
import random
import numpy as np
from numpy.random import rand
#-------------------------------------------
L = 5           # Size of the square lattice  
mcstp = 100     # Number of MC sweep
eqstp = 100     # Number of MC sweep for equilibration
D = 2           # Lattice dimension
N = L**D        # Number of spins
NN = 4          # Number of nearest-neighbors
T_0 = 0.5       # Initial temperature
T_f = 10.0      # Final temperature
dT = 0.1        # Temperature step
#----------------------------------
c = 0           # Counter
T = T_0         # Temperature initialization
norm = 1./(mcstp*N*NN)
nr2 = norm/mcstp
#---------------
def INIT(L):
	""" Generates initial confuguration (cfg): all spins up  """
	cfg = np.random.randint(1, size = (L, L)) + 1
	return cfg
#--------------------------
def METROPOLIS(spin, beta):
	"""The Metropolis algorithm"""
	for i in range(L): 
		for j in range(L):
			x = np.random.randint(0, L) # X coordinate
			y = np.random.randint(0, L) # Y coordinate
			s = spin[x, y]
			# Periodic boundary conditions --------------
			R = spin[(x + 1)%L, y] + spin[x, (y + 1)%L] \
			+ spin[(x - 1)%L, y] + spin[x, (y - 1)%L]
			#--------------------------------------------
			dE = 2*s*R # Energy difference
			if dE < 0.:
				s *= -1 # Flips the spin
			elif rand() < np.exp(-dE*beta):	# Throws the die
				s *= -1 # Flips the spin
			spin[x, y] = s # Returns the spin without flipping
	return spin
#--------------------------
def CALCULATE_ENERGY(spin):
	"""Computes energy"""
	En = 0. # Energy
	for x in range(len(spin)):
		for y in range(len(spin)):
			S = spin[x, y]
			# Periodic boundary conditions --------------
			R = spin[(x + 1)%L, y] + spin[x, (y + 1)%L] \
				+ spin[(x - 1)%L, y] + spin[x, (y - 1)%L]
			#--------------------------------------------			
			En -= S*R
	return En
#---------------------------------
def CALCULATE_MAGNETIZATION(spin):
	"""Computes magnetization"""
	mgnt = np.sum(spin) # Magnetization
	return mgnt
#-----------------------------------
# Opens the output file "output.txt" 
#  containing beta (B= 1/k_{B}T), 
#  temperature (T), and the mean values of energy (E), 
#  magnetization (M), specific heat (C_v), 
#  and susceptibility (Chi) of the system 
#  at every temperature step.
g = open('output.txt', 'w')
g.write('# T^{-1}        T           E           M           C_v         Chi\n')
g.write('#=======    ========     ========    ========    ========    ========\n')
#--------------------------------------------------
while T <= T_f:		# Temperature loop (main loop)
	T += dT
	summ = [0, 0, 0, 0]	# Accumulators
	spin = INIT(L)
	c = c + 1
	B = 1./T; B2 = B*B
	for sweep in range(eqstp):  # Equilibration loop
		METROPOLIS(spin, B)
#---------------------------------------------
	for sweep in range(mcstp):  # Main MC loop
		METROPOLIS(spin, B)
		E = CALCULATE_ENERGY(spin)
		M = CALCULATE_MAGNETIZATION(spin)
		#---------------------------------
		summ[0] += E      # E accumulator
		summ[1] += E*E    # E^2 accumulator
		summ[2] += M      # M accumulator
		summ[3] += M*M    # M^2 accumulator
#-------------------------------------------------
# Thermodynamic averages at each temperature step
	mean_E = summ[0]*norm  # Mean E
	mean_M = summ[2]*norm  # Mean M
	C_v = (norm*summ[1] - nr2*summ[0]*summ[0])*B2  # Specific heat (C_v)
	chi = (norm*summ[3] - nr2*summ[2]*summ[2])*B   # Susceptibility (Chi)
	g.write('%f    %f    %f    %f    %f    %f\n' \
	%(B, T, mean_E, mean_M, C_v, chi))
	print('>>> Temperature step %d' %(c))
