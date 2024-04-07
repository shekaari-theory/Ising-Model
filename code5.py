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
# (\Delta E is calculated here according 
#  to Eq. 10 of the main paper).
#
import random
import numpy as np
from numpy.random import rand
#-------------------------------------------
L = 5           # Size of the square lattice  
mcstp = 1000    # Number of MC sweeps
eqstp = 100     # Number of MC sweeps for equilibration
D = 2           # Lattice dimension
N = L**D        # Number of spins
NN = 4          # Number of nearest-neighbors
T = 2.5         # Temperature
B = 1./T        # Beta = 1/k_{B}T
B2 = B*B        # Beta^2
norm = 1./(mcstp*N*NN)
nr2 = norm/mcstp
#----------------
def INIT(L):
	""" Generates initial confuguration (cfg): all spins up  """
	cfg = np.random.randint(1, size = (L, L)) + 1
	return cfg
#--------------
def bw(dE):
	'''Computes Boltzmann weights'''
	if dE == -8 :
		return 24.532530197     # exp(-dE/T)=exp(8/2.5)
	elif dE == -4 :
		return 4.953032424      # exp(-dE/T)=exp(4/2.5)
	elif dE == 0 :
		return 1.0              # exp(-dE/T)=exp(0)
	elif dE == 4 :
		return 0.201896518      # exp(-dE/T)=exp(-4/2.5)
	elif dE == 8 :
		return 0.040762204      # exp(-dE/T)=exp(-8/2.5)
#-------------------------------------------------------
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
			dE = int(2*s*R) # Energy difference
			if dE < 0.:
				s *= -1 # Flips the spin
			elif rand() < bw(dE):	# Throws the die
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
#--------------------------------------
summ = [0, 0, 0, 0, 0]	# Accumulators
spin = INIT(L)
B = 1./T; B2 = B*B
# Equilibration loop -----
for sweep in range(eqstp):  
	METROPOLIS(spin, B)
#------------------------------------------------------
# Storing initial configuration in the file "spin.xyz"
# The file will contain all spins at every sweep,
# and will be used as the input file for XMakemol.
f = open('spin.xyz', 'w')
f.write("%d\n\n" %(N))
for x in range(L):
	for y in range(L):
		if spin[x, y] == 1:
			f.write("O %d %d\n" %(x, y))
               # O means oxygen: XMakemol parameter
		elif spin[x, y] == -1:
			f.write("Na %d %d\n" %(x, y))
               # Na means sodium: XMakemol parameter
#----------------------------------------------------
for sweep in range(mcstp):  # Main MC loop
	METROPOLIS(spin, B)
	E = CALCULATE_ENERGY(spin)
	M = CALCULATE_MAGNETIZATION(spin)
	#---------------------------------
	summ[0] += E      # E accumulator
	summ[1] += E*E    # E^2 accumulator
	summ[2] += M      # M accumulator
	summ[3] += M*M    # M^2 accumulator
	summ[4] += abs(M) # |M| accumulator
	#-----------------------------------
	f = open('spin.xyz', 'a')
	for x in range(L):
		for y in range(L):
			if spin[x, y] == 1:
				f.write("O %d %d\n" %(x, y))
			elif spin[x, y] == -1:
				f.write("Na %d %d\n" %(x, y))
#--------------------------------------------
# Computes thermodynamic averages
mean_E = summ[0]*norm      # <E>
mean_E2 = summ[1]*norm     # <E^2> 
mean_M = summ[2]*norm      # <M> 
mean_M2 = summ[3]*norm     # <M^2> 
mean_abs_M = summ[4]*norm  # <|M|> 
C_v = (mean_E2 - nr2*summ[0]*summ[0])*B2  # Specific heat
chi = (mean_M2 - nr2*summ[2]*summ[2])*B   # Susceptibility
print('\n -----------------------------------------------------------------')
print(' >>> THERMODYNAMIC AVERAGES:\n')
print('                              Lattice size  =  %d*%d' %(L, L))
print('                           Temperature [T]  =  %f' %(T))
print('                         Mean energy [<E>]  =  %f J' %(mean_E))
print('                                     <E^2>  =  %f J^2' %(mean_E2))
print('                   Mean manetization [<M>]  =  %f' %(mean_M))
print('                                     <M^2>  =  %f' %(mean_M2))
print('                                     <|M|>  =  %f' %(mean_abs_M))
print('        Energy fluctuation [<E^2> - <E>^2]  =  %f J^2' %(mean_E2 - nr2*summ[0]*summ[0]))
print(' Magnetization fluctuation [<M^2> - <M>^2]  =  %f' %(mean_M2 - nr2*summ[2]*summ[2]))
print('               Specific heat [<C_v>/k_{B}]  =  %f' %(C_v))
print('                    Susceptibility [<Chi>]  =  %f' %(chi))
print(' -----------------------------------------------------------------\n')
