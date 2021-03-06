import os
import sys
import time
import pprint
import subprocess
import numpy as np
import pprint
import pickle
from math import log

amino_acids_list = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",\
		    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

def parse_alignment_file(aligned_fasta_file):
	"""
	- Goes through an aligned FASTA file, extracts all the sequences,
	determines which of the sequence positions have residues (i.e.,
	not gaps) in 50% or more of the sequences, henceforth referred to
	as "wel-represented positions"
	
	INPUT: a FASTA file of a multiple sequence alignment (MSA)
	RETURNS: a tuple (seqs_byposition_dict, human_sequence), where
	(1) seqs_byposition_dict is a dictionary
	that contains each of the well-represented positions and the
	residue at that position in each of the sequences:
	seqs_byposition_dict[pos][seq_name], and
	(2) human_sequence is a string of the aligned human sequence
	for later analysis
	"""
	sequences_dict = {}

	with open(aligned_fasta_file, "r") as f:
		# since sequences span multiple lines, we'll keep track of
		# one as we're going through in prev_line_seq, and then
		# deposit it in the dictionary when we hit the next line
		# starting with ">"
		prev_line_seq = ""
		for line in f:
			# if line begins with ">", it's the name of a sequence,
			# so we want to store the previous sequence in the
			# dictionary and start a new record for this sequence
			if (line.strip()[0] == ">"):
				if (len(prev_line_seq) > 0):
					# add the previous sequence to the dict
					sequences_dict[seq_name] = prev_line_seq
					prev_line_seq = ""
				seq_name = line.strip()[1:].split("/")[0]

			# if line doesn't begin with ">", it's sequence that
			# we want to keep track of
			else:
				prev_line_seq += line.strip()

	# add the final sequence to the dictionary
	sequences_dict[seq_name] = prev_line_seq

	# Now go through each of the positions in the alignment and determine
	# which positions have residues (i.e., non-gaps) in at least 50% of
	# sequences - henceforth referred to as "well-represented positions" -
	# we'll only consider these positions going forward
	seq_names = sequences_dict.keys()
	length_of_seqs = len(sequences_dict[seq_names[0]])
	seqs_byposition_dict = {}

	plotted_position = 0
	for pos in range(length_of_seqs):
		tot_positions_with_res = 0
		for seq_name in seq_names:
			# see if that sequence has an amino acid or gap at pos
			if sequences_dict[seq_name][pos] in amino_acids_list:
				tot_positions_with_res += 1
		if (float(tot_positions_with_res)/len(seq_names) > 0.5):
			seqs_byposition_dict[pos] = {}
			# go through and include each of the sequences for that
			# position
			for seq_name in seq_names:
				seqs_byposition_dict[pos][seq_name] =\
				    sequences_dict[seq_name][pos]

	# get the human_sequence ("CGL_HUMAN") for later analysis
	human_sequence = sequences_dict["CGL_HUMAN"]
	
	return (seqs_byposition_dict, human_sequence)






def get_mono_freqs_and_di_counts(seqs_byposition_dict):
	"""
	- INPUT: seqs_byposition_dict, a dictionary returned by 
		the parse_alignment_file() function above
	- RETURNS: a tuple of (mono_freqs_dict, di_counts_dict),
		two dictionaries that have
		(1) the (mono-)amino acid frequencies at each of the well-
			represented positions
		(2) the di-amino acid counts at each pair of well-represented
			positions
			- format is di_counts_dict[pos][sec_pos][aa_1 + aa_2],
				where pos is the left column, sec_pos is the
				right column, aa_1+aa_2 is a string of the di-a.a.
	"""
	
	# a list of the "well-represented positions"
	positions_list = seqs_byposition_dict.keys()
	positions_list.sort()
	seq_names_list = seqs_byposition_dict[positions_list[0]].keys()
	
	# the dictionaries that contain the frequences of mono- and di-
	# amino acids, which will be returned (the _counts_ dictionaries are
	# temporary dictionaries to store counts, which will then be normalized
	# in into _freqs_)
	mono_freqs_dict = {}
	di_counts_dict = {}
	for pos_num, pos in enumerate(positions_list):
		print "Calculating amino acid frequencies for original MSA position ",\
		    pos," of", positions_list[-1],\
		    " (well-represented position #", pos_num, ")"
		mono_counts_dict = {}
		mono_freqs_dict[pos] = {}
		di_counts_dict[pos] = {}
		
		# initialize the count and frequency to be 0 for each amino acid
		for aa in amino_acids_list:
			mono_counts_dict[aa] = 0
			mono_freqs_dict[pos][aa] = 0
		# go through each of the sequences and add the counts
		non_gaps_at_pos = 0
		for seq_name in seq_names_list:
			aa = seqs_byposition_dict[pos][seq_name]
			# make sure that it's an amino acid (not a gap)
			if aa in amino_acids_list:
				non_gaps_at_pos += 1
				mono_counts_dict[aa] += 1
		# normalize the counts into frequencies
		for aa in amino_acids_list:
			mono_freqs_dict[pos][aa] =\
			    mono_counts_dict[aa]/float(non_gaps_at_pos)
		
		# Get the di-a.a. counts by iterating through
		# all positions after pos to the end as the second position of
		# the di-a.a.
		for sec_pos in positions_list[(pos_num+1):]:
			di_counts_dict[pos][sec_pos] = {}
			for aa_1 in amino_acids_list:
				for aa_2 in amino_acids_list:
					di_counts_dict[pos][sec_pos][aa_1+aa_2]=0
			
			# go through each of the sequences and get the aa pair
			# in positions (pos, sec_pos)
			for seq_name in seq_names_list:
				aa_1 = seqs_byposition_dict[pos][seq_name]
				aa_2 = seqs_byposition_dict[sec_pos][seq_name]
				# make sure that both are amino acids (not gaps)
				if (aa_1 in amino_acids_list) and\
					    (aa_2 in amino_acids_list):
					di_counts_dict[pos][sec_pos][aa_1+aa_2]+=1
	
	return (mono_freqs_dict, di_counts_dict)





def part_A_get_information_content(mono_freqs_dict):
	"""
	- Calculates the information content at each position and returns
	a list, info_content_list, that contains tuples
	[(orig_pos_in_seq, info), (orig_pos_in_seq, info), ...],
	where orig_pos_in_seq is the index in the original MSA
	(i.e., the keys of mono_nts_freqs_dict) and info is a float of the
	information content at that position
	"""
	info_content_list = []

	########## ADD YOUR CODE HERE ############
	for pos in mono_freqs_dict:
		info = log(20, 2)
		pos_freqs = mono_freqs_dict[pos]
		for aa in pos_freqs:
			if pos_freqs[aa] != 0:
				info += pos_freqs[aa]*log(pos_freqs[aa], 2)
		info_content_list.append((pos, info))
	########## END YOUR CODE HERE ############

	return info_content_list
				


def plot_info_content_list(info_content_list):
	"""
	- This will create 2 plots of the information content at each position,
		one indexed by the original MSA position and one indexed by the
		well-represented position #
		- matplotlib must be installed; otherwise, use the output
			from the above function (info_content_list) to plot
			the information at each position in Excel or another
			package of your choice
	"""
	import matplotlib
	matplotlib.use('agg')
	import matplotlib.pyplot as plt

	print "Making 2 plots of the information content at each position..."

	max_info = -1
	min_info = 1000
	for tupl in info_content_list:
		info = tupl[1]
		if (info > max_info):
			max_info = info
		if (info < min_info):
			min_info = info

	# makes a figure, axis, and plots the points on it relative to the
	# well-positioned index
	fig = plt.figure()
	ax = fig.add_subplot(111)
	#plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	OrRd = plt.cm.OrRd
	for plot_pos, info_tupl in enumerate(info_content_list):
		info = info_tupl[1]
		plt.scatter(plot_pos, info,\
				color=OrRd((info-min_info)/(max_info-min_info)),\
				s = 10, edgecolor = 'none')

	sm = plt.cm.ScalarMappable(cmap=plt.cm.OrRd,\
		norm=plt.normalize(vmin=min_info, vmax=max_info))
	# include a bar that labels the colors with their information content
	sm._A = []
	ticks_list = []
	ticks_labels = []
	for i in range(11):
		ticks_list.append( min_info + i*(max_info - min_info)/10. )
		ticks_labels.append( "{0:.2f}".format(min_info+\
			i*(max_info - min_info)/10.))
	pp = plt.colorbar(sm, ticks = ticks_list)
	pp.ax.get_yaxis().set_ticklabels(ticks_labels)

	# set the title and labels
	ax.set_title("Information scatterplot")	
	ax.set_xlim((0, len(info_content_list) + 1))
	ax.set_ylim((0, 1.1*max_info))
	ax.set_xlabel("Well-represented Position",fontsize=14)
	ax.set_ylabel("Information",fontsize=14)
	# this will save the figure in the directory from which this script
	# is run
	fig.savefig("information_plot_well_represented_positions.pdf")
	print "\n\tPlot saved as information_plot_well_represented_positions.pdf"

	# makes a figure, axis, and plots the points on it relative to the
	# original MSA positions
	fig = plt.figure()
	ax = fig.add_subplot(111)
	#plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	OrRd = plt.cm.OrRd
	for plot_pos, info_tupl in enumerate(info_content_list):
		pos = info_tupl[0]
		info = info_tupl[1]
		plt.scatter(pos, info,\
				color=OrRd((info-min_info)/(max_info-min_info)),\
				s = 10, edgecolor = 'none')

	sm = plt.cm.ScalarMappable(cmap=plt.cm.OrRd,\
		norm=plt.normalize(vmin=min_info, vmax=max_info))
	# include a bar that labels the colors with their information content
	sm._A = []
	ticks_list = []
	ticks_labels = []
	for i in range(11):
		ticks_list.append( min_info + i*(max_info - min_info)/10. )
		ticks_labels.append( "{0:.2f}".format(min_info +\
			i*(max_info - min_info)/10.) )
	pp = plt.colorbar(sm, ticks = ticks_list)
	pp.ax.get_yaxis().set_ticklabels(ticks_labels)

	# set the title and labels
	ax.set_title("Information scatterplot")	
	ax.set_xlim((0, pos + 1))
	ax.set_ylim((0, 1.1*max_info))
	ax.set_xlabel("Original MSA Position",fontsize=14)
	ax.set_ylabel("Information",fontsize=14)
	# this will save the figure in the directory from which this script
	# is run
	fig.savefig("information_plot_original_MSA_positions.pdf")
	print "\n\tPlot saved as information_plot_original_MSA_positions.pdf"
	plt.clf()


def get_MI_at_pairs_of_positions(di_counts_dict):
	"""
	- Calculates the MI at each pair of positions (i, j) in the multiple
	sequence alignment, where i < j are indexes corresponding to positions
	in the original MSA (i.e., are keys of mono_freqs_dict)
	- Returns a dictionary, mutual_information_dict, where
	mutual_information_dict[i][j] = MI_at_ij_positions
	- Note: when calculating the mutual information for a pair of positions,
		you can scale the di-a.a. counts directly into frequencies to
		get the joint distribution f^(i,j)_x,y, and you should sum over
		this joint distribution to get the marginal distributions
		f^(i)_x and f^(j)_y
	"""

	mutual_information_dict = {}
	
	########## ADD YOUR CODE HERE ############
	di_freqs_dict = {}
	pos1_freqs_dict = {}
	pos2_freqs_dict = {}
	total_count = 0
	
	for pos1 in di_counts_dict:
		for pos2 in di_counts_dict[pos1]:
			pos_counts = di_counts_dict[pos1][pos2]
			for di_aa in pos_counts:
				total_count += pos_counts[di_aa]

	for pos1 in di_counts_dict:
		di_freqs_dict[pos1] = {}
		pos1_freqs_dict[pos1] = {}
		pos2_freqs_dict[pos1] = {}
		for pos2 in di_counts_dict[pos1]:
			di_freqs_dict[pos1][pos2] = {}
			pos_counts = di_counts_dict[pos1][pos2]
			for di_aa in pos_counts:
				di_freqs_dict[pos1][pos2][di_aa] = pos_counts[di_aa]/float(total_count)

	for pos1 in di_freqs_dict:
		for pos2 in di_freqs_dict[pos1]:
			for di_aa in di_freqs_dict[pos1][pos2]:
				if di_aa[0] in pos1_freqs_dict[pos1]:
					pos1_freqs_dict[pos1][di_aa[0]] += di_freqs_dict[pos1][pos2][di_aa]
				else:
					pos1_freqs_dict[pos1][di_aa[0]] = di_freqs_dict[pos1][pos2][di_aa]
				if di_aa[1] in pos2_freqs_dict[pos2]:
					pos2_freqs_dict[pos2][di_aa[1]] += di_freqs_dict[pos1][pos2][di_aa]
				else:
					pos2_freqs_dict[pos2][di_aa[1]] = di_freqs_dict[pos1][pos2][di_aa]

	max_mi = 0.
	max_mi_pos = (0,0)

	for pos1 in di_counts_dict:
		mutual_information_dict[pos1] = {}
		for pos2 in di_counts_dict[pos1]:
			mi_val = 0
			for aa1 in amino_acids_list:
				for aa2 in amino_acids_list:
					joint = di_freqs_dict[pos1][pos2][aa1+aa2]
					divisor = pos1_freqs_dict[pos1][aa1]*pos2_freqs_dict[pos2][aa2]
					if divisor != 0 and joint !=0:
						mi_val += joint*log((joint/divisor), 2)
			mutual_information_dict[pos1][pos2] = mi_val
			if (mi_val > max_mi):
				max_mi = mi_val
				max_mi_pos = (pos1, pos2)
			#print 'mi of pos1:',pos1,' pos2:',pos2,' is = ', mi_val

	print 'max MI: ', max_mi, ' at pos1:',max_mi_pos[0], ' and at pos2:',max_mi_pos[1]

	########## END YOUR CODE HERE ############
	
	return mutual_information_dict





def part_c_get_highest_MI_block_of_10(mutual_information_dict):
	"""
	- Finds the block of 10 consecutive position pairs (i.e., keys of
		mutual_information_dict) with the highest average MI:
		between (pos_1, pos_2)
			   (pos_2, pos_3)
			   ...
			   ...
			   (pos_9, pos_10)		
	- Your code should return the list
		[pos_1, pos_2, ..., pos_10], where each of the pos_i
		are consecutive keys of mutual_information_dict
	"""

	highest_average_positions_list = []

	########## ADD YOUR CODE HERE ############
	positions = mutual_information_dict.keys()
	max_group_MI = 0
	highest_average_positions_list = []
	for i in range(len(positions)-9):
		avgMI = 0
		for j in range(9):
			avgMI += mutual_information_dict[positions[i+j]][positions[i+j+1]]
		print avgMI
		if avgMI > max_group_MI:
			pos_list = []
			for j in range(10):
				pos_list.append(positions[i+j])
			highest_average_positions_list = pos_list


		
	########## END YOUR CODE HERE ############
	
	return highest_average_positions_list




def plot_mutual_information(mutual_information_dict):
	"""
	- This will create a plot of the mutual information
		at each pair of positions
		- matplotlib must be installed; otherwise, use the output
			from the above function (mutual_information_dict)
			to plot the information at each position pair in
			Excel or another package of your choice
	"""
	import matplotlib
	matplotlib.use('agg')
	import matplotlib.pyplot as plt
	from matplotlib.patches import Rectangle
	positions_list = mutual_information_dict.keys()
	positions_list.sort()
	min_pos = min(positions_list)
	max_pos = max(positions_list)
	num_positions = len(positions_list)

	array_to_plot = []
	for i in range(num_positions):
		temp_list = [0] * num_positions
		array_to_plot.append(temp_list)
	
	# get the maximum and minimum MIs in the dictionary
	min_mi = 100000
	max_mi = -1
	for pos_num, pos in enumerate(positions_list):
		for sec_pos_num, sec_pos in enumerate(positions_list):
			try:
				mi = mutual_information_dict[pos][sec_pos]
				array_to_plot[pos_num][sec_pos_num] = mi
				if (mi < min_mi):
					min_mi = mi
				if (mi > max_mi):
					max_mi = mi
			except KeyError:
				if (pos != sec_pos):
					mi = mutual_information_dict[sec_pos][pos]
					array_to_plot[pos_num][sec_pos_num] = mi
	# makes a figure, axis, and plots the mutual information between
	# residues i and j as a box at (i, j) with the color of the box
	# indicating the Mutual Information; here (i, j) are the #s
	# corresponding to what well-represented positions it is
	print "\nPlotting indexed relative to well-represented positions..."
	fig = plt.figure();
	ax = fig.add_subplot(111);
	#plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	this_cmap = plt.cm.Reds
	plt.imshow(np.array(array_to_plot), cmap=this_cmap, vmin=0., vmax=max_mi)
	ax.set_aspect("equal")
	sm = plt.cm.ScalarMappable(cmap=this_cmap,\
			norm=plt.normalize(vmin=0., vmax=max_mi));
	# include a bar indicating what colors correspond to what MI values
	sm._A = []
	ticks_list = []
	ticks_labels = []	
	for i in range(11):
		ticks_list.append( i*(max_mi-min_mi)/10. )
		ticks_labels.append( "{0:.2f}".format(i*(max_mi-min_mi)/10.) )
	pp = plt.colorbar(sm, ticks = ticks_list)
	pp.ax.get_yaxis().set_ticklabels(ticks_labels)

	ax.set_title("Mutual information heatmap")	
	ax.set_ylim((0, len(positions_list)+1))
	ax.set_xlim((0, len(positions_list)+1))
	ax.set_xlabel("Well-represented Position 1",fontsize=14)
	ax.set_ylabel("Well-represented Position 2",fontsize=14)
	# this will save the figure in the directory from which this script is run
	fig.savefig("mutual_information_plot_well_represented_positions.pdf")
	print "\n\tPlot saved as mutual_information_plot_well_represented_positions.pdf"
	plt.clf()

	## does the same thing, but plots according to the original position in
	## the multiple sequence alignment, not which "well-represented position"
	## it is
	array_to_plot = []
	for i in range(max_pos+1):
		temp_list = [0] * (max_pos+1)
		array_to_plot.append(temp_list)
	
	# get the maximum and minimum MIs in the dictionary
	min_mi = 100000
	max_mi = -1
	for pos_num, pos in enumerate(positions_list):
		for sec_pos_num, sec_pos in enumerate(positions_list):
			try:
				mi = mutual_information_dict[pos][sec_pos]
				array_to_plot[pos][sec_pos] = mi
				if (mi < min_mi):
					min_mi = mi
				if (mi > max_mi):
					max_mi = mi
			except KeyError:
				if (pos != sec_pos):
					mi = mutual_information_dict[sec_pos][pos]
					array_to_plot[pos][sec_pos] = mi
	# makes a figure, axis, and plots the mutual information between
	# residues i and j as a box at (i, j) with the color of the box
	# indicating the Mutual Information; here (i, j) are the #s
	# corresponding to what well-represented positions it is
	print "\nPlotting indexed relative to original MSA positions..."
	fig = plt.figure();
	ax = fig.add_subplot(111);
	#plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	this_cmap = plt.cm.Reds
	plt.imshow(np.array(array_to_plot), cmap=this_cmap, vmin=0., vmax=max_mi)
	ax.set_aspect("equal")
	sm = plt.cm.ScalarMappable(cmap=this_cmap,\
			norm=plt.normalize(vmin=0., vmax=max_mi));
	# include a bar indicating what colors correspond to what MI values
	sm._A = []
	ticks_list = []
	ticks_labels = []	
	for i in range(11):
		ticks_list.append( i*(max_mi-min_mi)/10. )
		ticks_labels.append( "{0:.2f}".format(i*(max_mi-min_mi)/10.) )
	pp = plt.colorbar(sm, ticks = ticks_list)
	pp.ax.get_yaxis().set_ticklabels(ticks_labels)

	ax.set_title("Mutual information heatmap")	
	ax.set_ylim((min_pos, max_pos + 1))
	ax.set_xlim((min_pos, max_pos + 1))
	ax.set_xlabel("Original MSA Position 1",fontsize=14)
	ax.set_ylabel("Original MSA Position 2",fontsize=14)
	# this will save the figure in the directory from which this script
	# is run
	fig.savefig("mutual_information_plot_original_MSA_positions.pdf")
	print "\n\tPlot saved as mutual_information_plot_original_MSA_positions.pdf"







if __name__ == "__main__":
	### the FASTA MSA file should be the command line argument
	aligned_fasta_file = sys.argv[1]

	### parse the FASTA MSA file to get back a dictionary which gives
	### the positions that have residues (non-gaps) in 50%
	### or more of the sequences and the residue at each of those positions
	### in all of the sequences
	###    - also returns the human sequence in the MSA for future analysis
	seqs_byposition_dict, human_sequence = \
	    parse_alignment_file(aligned_fasta_file)

	### get the frequency of a.a.s and counts of di-a.a.s - see the function
	### for the structure of how these frequencies are represented in the
	### dictionaries that are returned
	mono_freqs_dict, di_counts_dict =\
	get_mono_freqs_and_di_counts(seqs_byposition_dict)
	
	############ PART (A) ############

	### get the information content at each position
	###	- you must implement part_A_get_information_content()
	info_content_list = part_A_get_information_content(mono_freqs_dict)
	max_info = max(info_content_list,key=lambda item:item[1])
	print 'max info: ', max_info[1], ' at pos: ', max_info[0]

	### you can uncomment the line below to plot the
	### information content at each position if you have matplotlib
	### installed; otherwise, plot info_content_list with some other
	### tool of your choice
	plot_info_content_list(info_content_list)
	############# END PART (A) ############



	############# PART (B) ############
	### Calculate the mutual information between all pairs of
	### "well-represented positions"
	###	- you must implement get_MI_at_pairs_of_positions()
	mutual_information_dict = get_MI_at_pairs_of_positions(di_counts_dict)
	
	### This will "pickle" mutual_information_dict (this is Python-speak
	### for saving the dictionary in binary format - "wb" corresponds
	### to Write in Binary mode - as mutual_information_dict.pkl)
	### so you can later load it without having to
	### run the above commands to regenerate it
	with open("mutual_information_dict.pkl", "wb") as f:
		pickle.dump(mutual_information_dict, f)

	### This will make 2 heatmap plots of the mutual information
	### dictionary, one indexed by the
	### original positions in the MSA and a condensed one indexed by the
	### "well-represented positions"
	### - only uncomment and use this if you have matplotlib installed
	plot_mutual_information(mutual_information_dict)
	
	############ END PART (B) ############


	############ PART (C) ############
	### determine which residues in the human sequence correspond
	### to the 10 consecutive positions with the highest avg. MI

	### if you have previously pickled mutual_information_dict,
	### you can uncomment this to load it without having to run
	### the above functions!
	#with open("mutual_information_dict.pkl") as f:
	#	mutual_information_dict = pickle.load(f)

	### get the 10 consecutive positions with the highest avg. MI
	###	- you must implement part_c_get_highest_MI_block_of_10
	highest_MI_block_of_10 = part_c_get_highest_MI_block_of_10(\
		mutual_information_dict)

	### This will print out
	### the MSA entries of the human sequence corresponding to the range
	### implied by the returned highest_MI_block_of_10
	if (len(highest_MI_block_of_10) == 10):
		start_of_MI_block = highest_MI_block_of_10[0]
		end_of_MI_block = highest_MI_block_of_10[-1]
		print "\nResidues ",start_of_MI_block,"-",end_of_MI_block,\
			" of the human sequence in the MSA are:"
		print human_sequence[start_of_MI_block:(end_of_MI_block+1)]

		print "\nNow go figure out where these residues are in the ungapped"
		print "human protein!"
	else:
		print "Fill in part_c_get_highest_MI_block_of_10() !!"	
	############ END PART (C) ############


