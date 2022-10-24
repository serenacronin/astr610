#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 19:56:40 2022

@author: Serena A. Cronin

This program calculates the black hole spin parameter as a function of 
black hole mass. Homework 3 for ASTR 680: High Energy Astrophysics.
"""

### ---- MY FUNCTIONS ---- ###

### ---- OTHER FUNCTIONS & PACKAGES ---- ###
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import time
from matplotlib.offsetbox import AnchoredText
startTime = time.time()  # TIME THE SCRIPT

##### --- DEFINE FUNCTIONS --- #####

def spec_ang_mom(direction, r, M, j):
	
	"""
	The specific angular momentum (i.e., per unit rest mass) of a 
	particle in a circular geodesic at radius r around a black hole 
	of mass M and spin parameter a = jM.
	
	
	direction : string; options: 'prograde' or 'retrograde'
	r : float; circular geodesic radius r around a black hole
	M : float; black hole mass
	j : float; dimensionless spin parameter
	
	Returns: u_phi
	
	"""
	
	a = j*M  # spin parameter
	
	if direction == 'prograde':
		
		num = np.sqrt(M*r) * (r**2 - 2*a*np.sqrt(M*r) + a**2)
		denom = r*(r**2 - 3*M*r + 2*a*np.sqrt(M*r))**(1/2)
		u_phi = num / denom
		
	elif direction == 'retrograde':
		
		num = np.sqrt(M*r) * (r**2 + 2*a*np.sqrt(M*r) + a**2)
		denom = r*(r**2 - 3*M*r - 2*a*np.sqrt(M*r))**(1/2)
		u_phi = -(num / denom)
		
	return u_phi
	
	
def spec_energy(direction, r, M, j):
	
	"""
	The specific energy (i.e., per unit rest mass) of a 
	particle in a circular geodesic at radius r around a black hole 
	of mass M and spin parameter a = jM.
	
	
	direction : string; options: 'prograde' or 'retrograde'
	r : float; circular geodesic radius r around a black hole
	M : float; black hole mass
	j : float; dimensionless spin parameter
	
	Returns: -u_t
	
	"""
	
	a = j*M  # spin parameter
	
	if direction == 'prograde':
		
		num = r**2 - 2*M*r + a*np.sqrt(M*r)
		denom = r*(r**2 - 3*M*r + 2*a*np.sqrt(M*r))**(1/2)
		u_t = num / denom
		
	elif direction == 'retrograde':
		
		num = r**2 - 2*M*r - a*np.sqrt(M*r)
		denom = r*(r**2 - 3*M*r - 2*a*np.sqrt(M*r))**(1/2)
		u_t = num / denom
		
	return u_t


def ISCO(direction, M, j):
	
	"""
	The location of the innermost stable circular orbit (ISCO)
	
	
	direction : string; options: 'prograde' or 'retrograde'
	M  : float; black hole mass
	j : float; dimensionless spin parameter
	
	Returns: r_isco
	
	"""
	
	# define quantities to be used in the ISCO calculation
	Z1 = 1 + (1 - j**2)**(1/3) * ((1 + j)**(1/3) + (1 - j)**(1/3))
	Z2 = (3*j**2 + Z1**2)**(1/2)
	
	
	if direction == 'prograde':
		r_isco = M*(3 + Z2 - ((3 - Z1)*(3 + Z1 + 2*Z2))**(1/2))
		
	elif direction == 'retrograde':
		r_isco = M*(3 + Z2 + ((3 - Z1)*(3 + Z1 + 2*Z2))**(1/2))
		
	return r_isco


def spin_param(J_old, u_phi, m_part, M_BH_old, u_t):
	
	"""
	
	Calculate the spin parameter as accretion of particles happens.
	This was not given in the problem and was calculated by hand.
	
	"""
	
	J_new = J_old + u_phi*m_part
	M_BH_new = M_BH_old + u_t*m_part
	
	j_new = J_new / M_BH_new**2
	
	return j_new


def run_jnew(M0, m_part, j0, j_upp, direction):
	
	"""
	
	Runs the calculations for a series of upper limits on j, black hole
	initial masses, mass of particle, and initial j
	
	"""
	
	M_list = []
	j_list = []
	j = j0
	M = M0
	
	if direction == 'prograde':
		while j < j_upp:
			
			# calculate ISCO, u_phi, and u_t
			r = ISCO(direction, M, j)
			u_phi = spec_ang_mom(direction, r, M, j)
			u_t = spec_energy(direction, r, M, j)
			
			# calculate the angular momentum and the new spin parameter
			J_old =  j * M**2
			j_new = spin_param(J_old, u_phi, m_part, M, u_t)
			
			# save to a list
			M_list.append(M)
			j_list.append(j_new)
			
			# update new black hole mass and spin parameter
			M = M + m_part
			j = j_new
			
	elif direction == 'retrograde':
		while j > j_upp:
			
			# calculate ISCO, u_phi, and u_t
			r = ISCO(direction, M, j)
			u_phi = spec_ang_mom(direction, r, M, j)
			u_t = spec_energy(direction, r, M, j)
			
			# calculate the angular momentum and the new spin parameter
			J_old =  j * M**2
			j_new = spin_param(J_old, u_phi, m_part, M, u_t)
			
			# save to a list
			M_list.append(M)
			j_list.append(j_new)
			
			# update new black hole mass and spin parameter
			M = M + m_part
			j = j_new
		
	return M, j, M_list, j_list

##### --- SET UP PLOT --- #####
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["axes.linewidth"] = 2.5
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.variant"] = "small-caps"

##### --- DEFINE VARIABLES --- #####

run_partA = True
run_partB = True

##### -------------- PART A --------------- #####

if run_partA == True:
	direction = 'prograde'
	M0 = 1
	m_part = 0.001  # want the calculation within 0.1%
	j0 = 0
	
	# run for 0.5, 0.9, 1.0
	j_upps = [0.5, 0.9, 0.999999, 0.99999999]
	
	fig = plt.figure(figsize=(13,13))
	p = 1
	for j_upp in j_upps:
		M, j, M_list, j_list = run_jnew(M0, m_part, j0, j_upp, direction)
		
		ax = plt.subplot(2, 2, p)
		ax.plot(M_list, j_list, color='tab:pink', lw=5)
		ax.set_xlabel('M', fontsize=20)
		ax.set_ylabel('j', fontsize=20)
		ax.set_title('j = %s' % j_upp, fontsize=20)
		ax.tick_params(axis='both', which='both',direction='in',
	                          width=2.5, labelsize=16, length=7)
		
		at = AnchoredText('M_final = %s' % round(M, 3), prop=dict(size=17), 
	                            frameon=True, loc='lower right')
		at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		ax.add_artist(at)
		p = p+1
		
	plt.show()

##### -------------- PART B --------------- #####
if run_partB == True:
	
	direction = 'retrograde'
	M0 = 1
	m_part = 0.001  # want the calculation within 0.1%
	j0 = 0.999
	
	# run for 0.5, 0.9, 1.0
	j_upps = [0.99, 0.9, 0.5, 0]
	
	fig = plt.figure(figsize=(13,13))
	sp = 1
	for j_upp in j_upps:
		M, j, M_list, j_list = run_jnew(M0, m_part, j0, j_upp, direction)
		
		ax = plt.subplot(2, 2, sp)
		ax.plot(M_list, j_list, color='tab:cyan', lw=5)
		ax.set_xlabel('M', fontsize=20)
		ax.set_ylabel('j', fontsize=20)
		ax.set_title('j = %s' % j_upp, fontsize=20)
		ax.tick_params(axis='both', which='both',direction='in',
	                          width=2.5, labelsize=16, length=7)
		
		at = AnchoredText('M_final = %s' % round(M, 3), prop=dict(size=17), 
	                            frameon=True, loc='upper right')
		at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		ax.add_artist(at)
		sp = sp+1
		
	plt.show()

##### --- TIME THE SCRIPT --- #####
outfile = open('Time_.txt', 'w')
executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))
outfile.close()