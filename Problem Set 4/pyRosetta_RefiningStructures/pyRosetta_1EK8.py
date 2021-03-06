import os
import sys
import time
import pprint
import subprocess
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import random
from math import exp
from numpy.random import normal

# import all functions from PyRosetta
from rosetta import *
init()

# this the standard score function which you should use to calculate the energy
# of any structures
scorefxn = get_fa_scorefxn()

# this is a function provided to you and called to make a simple plot of the
# energy at each iteration
def make_energy_plot(energy_list, out_file,\
		# all of the following are options which may or may not be
		# modified
		x_label = "Iteration",\
		title = "", plot_dots = True,\
		xticks_list = [], xticks_labels_list = []):
	"""
	- Makes a simple plot of energies contained in energy_list (for example,
		energies of a structure at each iteration)
	- Saves the plot as out_file in the current directory
	- Can optionally pass in an x-axis label, title, whether to include
		dots or not, and where to plot x-tick labels and what those
		labels should be
	NOTE: must have matplotlib installed (it is on Athena dialup)
	"""
	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	# plot the energy_list
	xs = range(len(energy_list))
	if (plot_dots == True):
		plt.plot(xs, energy_list, 'bo-')
	else:
		plt.plot(xs, energy_list, 'b-')

	# set title and axis labels
	ax.set_title(title)
	if (xticks_list != []):
		ax.set_xticks( xticks_list )
		ax.set_xticklabels(xticks_labels_list, size=12, rotation = 45)
	ax.set_xlabel(x_label, fontsize=14)
	ax.set_ylabel("Energy", fontsize=14)

	# this will save the figure in the directory from which this script is
	# run
	fig.savefig(out_file)
	print "\n\tPlot saved as: ", out_file



def part_a_energy_minimization():
	"""
	- Implement a simple energy minimization (gradient descent-like) algorithm:
		- For each residue's phi and psi angles, calculate the energy
			of the structure with that angle perturbed by + and - 1
			degree (4 possible changes per residue)
		- Accept the one change with the lowest energy (select 1 of the 
			(4*number_of_residue) structures)
		- Continue iterating until convergence, defined here as an energy
			change of less than 1.
	- Returns: the tuple (energy_list, energy_min_pose), where energy_list is a
		list that contains the starting energy followed by energy at each
		iteration, and energy_min_pose is a Pose of the final structure.

	"""
	
	energy_list = []
	
	energy_min_pose = pose_from_pdb("1EK8.rotated.pdb")
	########## ADD YOUR CODE HERE ############
	scorefxn = get_fa_scorefxn()
	energy_list.append(scorefxn(energy_min_pose))
	print 'starting energy: ', energy_list[0]
	energydiff = 2
	while energydiff > 1:
		residue_energies = []
		residue_angles = []
		residue_angle_vals = []
		for i in range(1, energy_min_pose.total_residue() + 1):
			phi = energy_min_pose.phi(i)
			psi = energy_min_pose.psi(i)
			energy_min_pose.set_phi( i, phi+1 )
			min_energy = scorefxn( energy_min_pose )
			min_angle = 'phi'
			min_angle_val = 1			

			energy_min_pose.set_phi( i, phi-1 )
                	energy = scorefxn( energy_min_pose )
                	if (energy < min_energy):
				min_energy = energy
				min_angle_val = -1
			energy_min_pose.set_phi(i, phi)
						
			energy_min_pose.set_psi( i, psi+1 )
                	energy = scorefxn( energy_min_pose )
                	if energy < min_energy:
				min_energy = energy
				min_angle_val = 1
				min_angle = 'psi'
			energy_min_pose.set_psi( i, psi-1 )
                	energy = scorefxn( energy_min_pose )
                	if energy < min_energy:
				min_energy = energy
				min_angle_val = -1
				min_angle = 'psi'
			residue_energies.append(min_energy)
			residue_angles.append(min_angle)
			residue_angle_vals.append(min_angle_val)
			energy_min_pose.set_psi(i, psi)
		min_residue = min(residue_energies)
		min_index = residue_energies.index(min_residue)
		energydiff = abs(energy_list[-1] - min_residue)
		energy_list.append(min_residue)
		if residue_angles[min_index] == 'phi':
			energy_min_pose.set_phi(min_index, energy_min_pose.phi(min_index)+residue_angle_vals[min_index])
		else:
			energy_min_pose.set_psi(min_index, energy_min_pose.psi(min_index)+residue_angle_vals[min_index])
		
	print 'final energy: ', energy_list[-1]
	print 'iterations: ', len(energy_list) - 1			
	########## END YOUR CODE HERE ############
	return (energy_list, energy_min_pose)


def part_c(kT):
	"""
	- Input: A float, kT (the temperature)
	- This function calls onesetofMetropolismoves(), which performs 100*x
		random perturbations (x = # of residues), accepting each
		according to the Metropolis criterion - see that function for
		more specific details
	"""	

	# load the rotated pose
	rotated_pose = pose_from_pdb("1EK8.rotated.pdb")

	# run one set (100*x) of Metropolis moves, where x is the number of
	# residues in the pose
	energy_list, new_pose = onesetofMetropolismoves(rotated_pose, kT)

	return energy_list



def onesetofMetropolismoves(pose, kT):
	"""
	- Implements one "set" of Metropolis moves (=100*x moves) on the pose,
		where x is the number of residues in the pose, and kT is a
		float of the temperature
	
	- You should make multiple calls to the function
		make_backbone_change_metropolis(), which
		makes one move according to the Metropolis acceptance criterion
	
	- Returns: a tuple of (energy_list, pose), where energy_list is a list
		of the energies after each of the moves, and pose is the
		final structure after the 100*x moves
	"""
	
	energy_list = []
	
	########## ADD YOUR CODE HERE ############
	for i in range(pose.total_residue()*100):
		pose, energy_of_pose = make_backbone_change_metropolis(pose, kT)
		energy_list.append(energy_of_pose)	
	
	########## END YOUR CODE HERE ############

	return energy_list, pose




def make_backbone_change_metropolis(pdb_pose, kT):
	"""
	- Performs one move according to the Metropolis acceptance criterion
	- Takes in a PyRosetta pose (pdb_pose) and float temperature (kT) and:
		1. Randomly selects one of its residues
		2. Randomly chooses either the phi or psi angle
		3. Chooses a perturbation angle drawn from a normal distribution
			with mean 0 and standard deviation 20
		4. Accepts the random angle change according to the
			Metropolis criterion
	- Returns: a tuple of (pose, energy_of_pose) for whichever pose is
		selected by the Metropolis criterion
	"""

	energy_of_pose = 0.

	########## ADD YOUR CODE HERE ############
	scorefxn = get_fa_scorefxn()
	energy_of_pose = scorefxn(pdb_pose)
	rand_residue = random.randint(1, pdb_pose.total_residue())
	rand_angle = random.randint(0, 1)
	rand_degrees = random.normalvariate(0,20)
	psi = pdb_pose.psi(rand_residue)
	phi = pdb_pose.phi(rand_residue)
	if rand_angle == 0:
		pdb_pose.set_phi(rand_residue, phi+rand_degrees)
	else:
		pdb_pose.set_psi(rand_residue, psi+rand_degrees)
	energy = scorefxn(pdb_pose)
	if energy < energy_of_pose:
		return (pdb_pose, energy)
	else:
		prob = exp(-(energy-energy_of_pose)/float(kT))
		rand_var = random.random()
		if rand_var < prob:
			return (pdb_pose, energy)
		else:
			if rand_angle == 0:
				pdb_pose.set_phi(rand_residue, phi)
			else:
				pdb_pose.set_psi(rand_residue, psi)
			return (pdb_pose, energy_of_pose)
	
	########## END YOUR CODE HERE ############


def part_d(starting_kT, number_of_kT_halvings):
	"""
	- Implements an annealing schedule starting at high temperature
		(starting_kT) and lowering slowly
	- At each temperature, perform 100*x moves according to the Metropolis
		criterion (hint: the onesetofMetropolismoves() function you
		implemented in part (C) should be helpful here!)
	- Decrease the temperature by a factor of 2 and repeat
	- Continue so that the last iteration performed is the first such that
		kT is less than 1 (so number_of_kT_halvings+1 total iterations)
	- Returns: a tuple of (full_energy_list, kT_byround_list, rotated_pose)
		- full_energy_list: a list of the energy at each iteration,
			where each onesetofMetropolismoves() consists of
			100*number_of_residues iterations
		- kT_byround_list: a list of kT used for each
			onesetofMetropolismoves()
		- rotated_pose: a Pose of the final structure after all rounds
			of Simulated Annealing
	- NOTE: this should take upwards of an hour when run for the full
		number_of_kT_halvings. Therefore, it's recommended that you run
		your code with number_of_kT_halvings=0 or 1 to make sure there
		are no errors, and then change number_of_kT_halvings to the
		value you've calculated and let it run to completion
	"""
	# load the rotated pose
	rotated_pose = pose_from_pdb("1EK8.rotated.pdb")

	# Lists to keep track of the energy at each iteration and the kT used
	# for each onesetofMetropolismoves()
	full_energy_list = []
	kT_byround_list = []

	########## ADD YOUR CODE HERE ############
	for i in range(number_of_kT_halvings+1):
		energy_list, rotated_pose = onesetofMetropolismoves(rotated_pose, starting_kT)
		full_energy_list.extend(energy_list)
		kT_byround_list.append(starting_kT)
		starting_kT /= 2
	########## END YOUR CODE HERE ############

	return full_energy_list, kT_byround_list, rotated_pose






def part_e_num_angles_different(pose_1, pose_2):
	"""
	- Takes in two structures (PyRosetta Poses), which must have the same number
		of residues
	- Returns: num_angles_different, the number of phi and psi angles
		that differ by more than 1 degree between the two
		(if there are x residues, a maximum of 2*x can be returned)
	- NOTE: you will want to consider corner cases (e.g. difference between
		12:59 and 1:00 is only one minute, even though it may appear as
		being ~12 hours different if you simply subtract the hours and
		minutes - what is the analogous case for phi and psi angles?) - 
		it may be helpful to write and call a separate helper function
		for this (can be written ouside of the indicated code bounds)
	"""

	num_residues_1 = pose_1.total_residue()
	num_residues_2 = pose_2.total_residue()
	# make sure the two structures have the same number of residues
	assert( num_residues_1 == num_residues_2 )
	
	num_angles_different = 0
	########## ADD YOUR CODE HERE ############
	for i in range(1, num_residues_1 + 1):
		anglediff = abs(pose_1.phi(i)-pose_2.phi(i))
		print 'phi diff: ',anglediff
		if (anglediff > 1 and not abs(360-anglediff) < 1):
			print 'differ by more than 1, phi'
			num_angles_different += 1
		anglediff = abs(pose_1.psi(i)-pose_2.psi(i))
		print 'psi diff: ',anglediff
		if (anglediff > 1 and not abs(360-anglediff) < 1):
			print 'differ by more than 1, psi'
			num_angles_different += 1
	########## END YOUR CODE HERE ############
	return num_angles_different






# If you run this script from the command line, the script will enter this loop
if __name__ == "__main__":
	
	### load the two Poses provided to you
	clean_pose = pose_from_pdb("1EK8.clean.pdb")
	rotated_pose = pose_from_pdb("1EK8.rotated.pdb")
	
	###### PART (A) ########
	
	### you must implement part_a_energy_minimization()
	energy_list, energy_min_pose = part_a_energy_minimization()
	
	### save the Pose as a .pdb file for future analysis
	#energy_min_pose.dump_pdb('1EK8.energy_min.pdb')
	### the function make_energy_plot is already implemented for you and
	### will automatically make a plot
	#make_energy_plot(energy_list, "energy_minimization_plot.pdf",\
		#title = "Energy Minimization at each iteration")
	
	##### END PART (A) #####



	
	###### PART (C) ########

	### you should assign kT_from_part_B
       	kT_from_part_B = 2.885
	
	### part_c() calls onesetofMetropolismoves(), which in turn
	### should call make_backbone_change_metropolis(); you need to 
	### implement onesetofMetropolismoves() and
	### make_backbone_change_metropolis()
       	#energy_list = part_c(kT_from_part_B)
	
	### Plots the energy at each iteration
	#metropolis_title = "Energy of Metropolis Algorithm, kT=", kT_from_part_B
	#make_energy_plot(energy_list, "energy_plot_metropolis.pdf",
		#title = metropolis_title, plot_dots = False)

	##### END PART (C) #####


	
	###### PART (D) ########
	
	### your calculation of the number of halvings such that the last iteration
	### is the first with kT less than 1
	number_of_kT_halvings = 2
	energy_list_sim_annealing, kT_byround_list, sim_annealing_pose =\
		part_d(kT_from_part_B, number_of_kT_halvings)	
	
	### write the final pose to a .pdb file for future analysis
	sim_annealing_pose.dump_pdb('1EK8.simulatedannealing.pdb')
	
	### The following code will plot the energy at each iteration, with kT
	### being labeled on the x-axis at the end of that set's moves; you do
	### not need to edit any of the code below here
	metropolis_title = "Energy of Metropolis Algorithm, initial kT={0}; {1} halvings of kT".format(\
		kT_from_part_B, number_of_kT_halvings)
	# the x-ticks (final move # for that kT)
	x_ticks = []
	rotated_pose = pose_from_pdb("1EK8.rotated.pdb")
	num_residues = rotated_pose.total_residue()
	for i in range(number_of_kT_halvings + 1):
		x_ticks.append((100*num_residues)*(i+1))
	# the x-tick labels  (the kT for the set of moves ending at that x_tick)
	x_ticks_labels = []
	for kT_that_round in kT_byround_list:
		x_ticks_labels.append(round(kT_that_round, 1))
	make_energy_plot(energy_list_sim_annealing, "energy_plot_metropolis_part_d.pdf",
		title = metropolis_title, plot_dots = False, xticks_list = x_ticks,
		xticks_labels_list = x_ticks_labels, x_label = "kT for previous round")
	
	###### END PART (D) ########



	###### PART (E) ########
	
	### you must implement part_e_num_angles_different()
	### number of angles different between ROTATED and clean poses
	num_angles_diff_in_rotated = part_e_num_angles_different(\
		rotated_pose, clean_pose)
	print num_angles_diff_in_rotated,\
		" angles different between Rotated and Clean poses"
	### number of angles different between ENERGY_MIN and clean poses
	energy_min_pose = pose_from_pdb("1EK8.energy_min.pdb")
	num_angles_diff_in_energy_min = part_e_num_angles_different(\
		energy_min_pose, clean_pose)
	print num_angles_diff_in_energy_min,\
		" angles different between Energy Min. and Clean poses"
	### number of angles different between SIMULATED ANNEALING and clean poses
	simulated_annealing_pose = pose_from_pdb("1EK8.simulatedannealing.pdb")
	num_angles_diff_in_sim_annealing_pose = part_e_num_angles_different(\
		sim_annealing_pose, clean_pose)
	print num_angles_diff_in_sim_annealing_pose,\
		" angles different between Simulated Annealing and Clean poses"

	### Get the RMSD between the three poses and the clean poses
	print "RMSD between Rotated and Clean poses: ",\
		CA_rmsd(rotated_pose, clean_pose) 
	print "RMSD between Energy Min. and Clean poses: ",\
		CA_rmsd(energy_min_pose, clean_pose) 
	print "RMSD between Simulated Annealing and Clean poses: ",\
		CA_rmsd(simulated_annealing_pose, clean_pose)
	
	###### END PART (E) ########


