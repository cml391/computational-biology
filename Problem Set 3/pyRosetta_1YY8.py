

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

# import all functions from PyRosetta
from rosetta import *
init()

def part_b():
	# load the PDB file into a "pose", Rosetta's data container object for a structure
	pdb_protein = pose_from_pdb("1YY8.clean.pdb")

	# print out the phi and psi angles of all the residues
	for i in range(1, pdb_protein.total_residue() + 1):
		print "Residue ", i, ": phi =", pdb_protein.phi(i), "psi =", pdb_protein.psi(i)

	# Get a scoring function, which by default includes Full-Atom (fa) energy terms
	scorefxn = get_fa_scorefxn()

	# print out the energy terms contributing to the pdb_protein
	print "\n\nPre-packing energy of protein:\n\n"
	scorefxn.show( pdb_protein )


def part_c():
	# load the PDB file into a "pose", Rosetta's data container object for a structure
	pdb_protein = pose_from_pdb("1YY8.clean.pdb")
	
	scorefxn = get_fa_scorefxn()

	# Implements Monte-Carlo Side-Chain packing. Note that this is packing
	# all residues, not just a particular residue as in the PyRosetta tutorial
	# - see http://graylab.jhu.edu/pyrosetta/downloads/documentation/Workshop6_PyRosetta_Packing_Design.pdf
	#		for a tutorial
	task_pack = standard_packer_task( pdb_protein )
	task_pack.restrict_to_repacking()
	task_pack.or_include_current(True) # include original side-chains
	pack_mover = PackRotamersMover(scorefxn, task_pack) 
	print "Applying the side-chain packing..."
	pack_mover.apply( pdb_protein ) # implement side-chain packing

	################ ADD CODE TO PRINT THE POST-PACKING SCORE ##################
	scorefxn = get_fa_scorefxn()
	scorefxn.show(pdb_protein)	
	
	################ END CODE TO PRINT THE POST-PACKING SCORE ##################

def part_d():
	print "Implementing side-chain packing..."	
	########### ADD CODE TO LOAD 1YY8.rotated.pdb and pack side chains ########
	# load the rotated PDB file into a "pose"
	pdb_protein = pose_from_pdb("1YY8.rotated.pdb")

	# print out the energy terms contributing to the protein before
	# side chain packing
	scorefxn = get_fa_scorefxn()
	scorefxn.show(pdb_protein)

	# Implement Monte-Carlo Side-Chain packing
	task_pack = standard_packer_task( pdb_protein )
        task_pack.restrict_to_repacking()
        task_pack.or_include_current(True) # include original side-chains
        pack_mover = PackRotamersMover(scorefxn, task_pack)
        print "Applying the side-chain packing..."
        pack_mover.apply( pdb_protein ) # implement side-chain packing

        # print out the energy terms contributing to the protein after
        # the side chain packing
        scorefxn = get_fa_scorefxn()
        scorefxn.show(pdb_protein)
	
	########### END CODE TO LOAD 1YY8.rotated.pdb and pack side chains ########

def part_e():
	pdb_protein = pose_from_pdb("1YY8.clean.pdb")
	# reload the rotated structure (without packed side-chains)
	rot_pdb_protein = pose_from_pdb("1YY8.rotated.pdb")
	################### ADD CODE TO DETERMINE WHICH ANGLE CHANGED #################
	# iterate through resides to determine which residue and angle (phi or psi) is changed
	# print out the phi and psi angles of all the residues
	difResidue = 0
        for i in range(1, pdb_protein.total_residue() + 1):
                print "Residue ", i, ": phi =", pdb_protein.phi(i), "psi =", pdb_protein.psi(i)
		print "Rotated Residue ", i, ": phi =", rot_pdb_protein.phi(i), "psi =", rot_pdb_protein.psi(i)
		if abs(pdb_protein.phi(i) - rot_pdb_protein.phi(i)) > 1 or abs(pdb_protein.psi(i) - rot_pdb_protein.psi(i)) > 1:
			difResidue = i
			print 'These residues are different!'
	print "Differing Angle Residue ", difResidue	
	################### END CODE TO DETERMINE WHICH ANGLE CHANGED #################

def part_f():
	rot_pdb_protein = pose_from_pdb("1YY8.rotated.pdb")
	scorefxn = get_fa_scorefxn()
	# for the residue and angle you determined from (e) was changed, calculate the energy
	# of the structure with each of the 360 possible integer angles. Add these
	# energies in the list energies_of_angle so it can be plotted.
	lowest_modified_energy_so_far = 1000000
	angle_with_lowest_energy_so_far = -181
	possible_angles = range(-180, 181)
	# you'll need to fill in this list so it's the same length as possible_angles;
	# then, a plot will be created below
	energies_of_angle = [0]*361
	### ADD CODE TO CALCULATE ENERGIES OF ALL POSSIBLE ANGLES FOR YOUR RESIDUE, ANGLE ###
	for i in range(361):
		rot_pdb_protein.set_phi(50, possible_angles[i])
		energies_of_angle[i] = scorefxn(rot_pdb_protein)
		if energies_of_angle[i] < lowest_modified_energy_so_far:
			lowest_modified_energy_so_far = energies_of_angle[i]
			angle_with_lowest_energy_so_far = possible_angles[i]
	### END CODE TO CALCULATE ENERGIES OF ALL POSSIBLE ANGLES FOR YOUR RESIDUE, ANGLE ###
	
	print "Lowest energy: ", lowest_modified_energy_so_far, " with angle: ", angle_with_lowest_energy_so_far
	# This will make a .pdf plot of the structure's energy vs. angle of the
	# structure once you've filled in the energies_of_angle list
	if (len( energies_of_angle ) == len( possible_angles )):
		
		max_energy = max( energies_of_angle )
		min_energy = min( energies_of_angle )
		# makes a figure, axis, and plots the points on it
		fig = plt.figure()
		ax = fig.add_subplot(111)
		plt.rc('font', family='serif')
		YellowOrangeRed = plt.cm.YlOrRd
		for i, y in enumerate(energies_of_angle):
			plt.scatter(possible_angles[i], y,
				color = YellowOrangeRed((max_energy-y)/(max_energy-min_energy)),
				s = 15, edgecolor = 'none')	
		# MODIFY THE TITLE to say which residue number and angle (phi or psi) this plot is showing
		ax.set_title("Energy of 1YY8 by 50 phi degree")	
		ax.set_ylabel("Energy",fontsize=14)
		ax.set_xlabel("Angle Degree",fontsize=14)
		# this will save the figure in the directory from which this script is run
		fig.savefig("energy_vs_angle.pdf")

	else:
		print "Length of energies_of_angle and possible_angles are not equal!", len(energies_of_angle)," != ", len(possible_angles)
		print "Write code to fill in energies_of_angle so it's as long as possible_angles"


# If you run this script from the command line, the script will enter this loop
if __name__ == "__main__":
	#part_b()
	#part_c()
	#part_d()
	#part_e()
	part_f()

