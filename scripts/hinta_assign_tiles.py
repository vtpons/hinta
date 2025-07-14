'''
Haplotype-based Inference of Non-transmitted Alleles (HINTA)

--------------------------------
Year: 2021

Authors:
- Hanna van Loo
- Annique Claringbould
- Priscilla Kamphuis
- Victória Trindade Pons

--------------------------------
Description:
In this script the haplotypes from the child are compared to the haplotypes from the parent(s)
to determine which haplotype is from which parent haplotype.
This works with trios (child and both parents) and with child-parent pairs (so one parent is present).
The goal of this script is to determine which haplotypes from the parent(s) are transmitted and which ones
are non-transmitted.
Recombination is also implemented in this script. There is one overall best matching parental haplotype
per child haplotype and if present also a recombination parental haplotype. (if some tiles are matching better with
another parental haplotype).

Input: 1) .haps file
	2) .sample file
	3) Tile size (default 150)
	4) Output filename .haps transmitted alleles file (example: chr20_transmitted.haps)
	5) Output filename .sample transmitted alleles (example: chr20_transmitted.sample)
	6) Output filename .haps non-transmitted alleles file (example: chr20_non_transmitted.haps)
	7) Output filename .sample non-transmitted alleles (example: chr20_non_transmitted.sample)

Output: 1) The overall best matching parental haplotype printed to the terminal per child with the matching percentage
	2) 2 .haps files: transmitted alleles and non-transmitted alleles, with corresponding .sample files
        (the sample files are identical)

--------------------------------
'''

#IMPORTS
import numpy as np
import pandas as pd
import argparse
import re
import csv
import collections

def argparser():
	'''
	This function is to parse the commandline arguments.
	Input: 1) .haps file (string)
		2) .sample file (string)
		3) tile size (integer)
		4) output transmitted .haps file name (string)
		5) output transmitted .sample file name transmitted alleles (string)
		6) output non-transmitted .haps file name non-transmitted alleles (string)
        7) output non-transmitted .sample file name non-transmitted alleles (string)
	Output: the commandline arguments
	'''
	parser = argparse.ArgumentParser(description="Determine non-transmitted alleles")
	parser.add_argument("haplo_file", type=str, help="Path to the .haps file to use") # .haps file
	parser.add_argument("sample_file", type=str, help="Path to the .sample file to use") # .sample file
	parser.add_argument("tile_size", type=int, help="Tile size to use, number of bases in tile") # Tile size
	parser.add_argument("output_transmitted_haps", type=str, help="Path and filename for output .map file") # output 1, map file
	parser.add_argument("output_transmitted_sample", type=str, help="path and filename for output .ped transmitted file") # output 2, transmitted .ped file
	parser.add_argument("output_nontransmitted_haps", type=str, help="path and filename for output .ped non-transmitted file") # output 3, non-transmitted .ped file
	parser.add_argument("output_nontransmitted_sample", type=str, help="path and filename for output .ped non-transmitted file") # output 3, non-transmitted .ped file
	args = parser.parse_args()
	haplo_file = args.haplo_file
	sample_file = args.sample_file
	tile_size = args.tile_size
	output_transmitted_haps = args.output_transmitted_haps
	output_transmitted_sample = args.output_transmitted_sample
	output_nontransmitted_haps = args.output_nontransmitted_haps
	output_nontransmitted_sample = args.output_nontransmitted_sample
	return haplo_file, sample_file, tile_size, output_transmitted_haps, output_transmitted_sample, output_nontransmitted_haps, output_nontransmitted_sample

##EDITED: removed ped_file_info, collects info for sample instead (sample_file_info)
def process_sample_file(sample_file):
	'''
	In this function the .sample file is used to fill the dictionaries children and individuals and also
	collecting the information for the .ped file (without the genotype).
	Input: .sample file format
	Output: A dict with all individuals and their index, a dict with all children and info for output .sample file (child ID with their parents IDs) 
    
	'''
	index = 0 #counter
	individuals = {} # Make dict with all individuals
	children = collections.OrderedDict() # Make dict for only the children
	sample_file_info = []
	sample_file = open(sample_file, "r")
	sample = sample_file.readlines()[2:] # Read lines without header (there are 2 header lines in a .sample file)
	for line in sample: # In this for loop fill the individuals dict
		line = line.strip().split(' ')
		individuals[line[1]] = index
		index += 1
	for line in sample:
		line = line.strip().split(' ')
		father = line[3] # Column with fathers ID
		mother = line[4] # Column with mothers ID 
		if (father != "0" and father in individuals) or (mother != "0" and mother in individuals):
			children[line[1]] = line[3:5] # Add individual to the dict of children
			tmp_dict = {}
			tmp_dict["ID_1"] = line[0]
			tmp_dict["ID_2"] = line[1]
			tmp_dict["missing"] = line[2]
			tmp_dict["father"] = line[3]
			tmp_dict["mother"] = line[4]
			tmp_dict["sex"] = line[5]
			tmp_dict["phenotype"] = line[6]
			sample_file_info.append(tmp_dict)	           
	return individuals, children, sample_file_info

def compare(child_hap, parent_hap):
	'''
	This function is to check the similarity of a child and parent haplotype tile
	Input: child haplotype and parental haplotype to compare
	Output: the matching percentage of the two compared haplotypes as a float
	'''
	if child_hap != '0' and child_hap != '1' and not isinstance(child_hap, list):
		child_hap = ord(child_hap) # Make ords of input if it is not given with 0s and 1s, to also compare nucleotides like A, C etc.
	if parent_hap != '0' and parent_hap != '1' and not isinstance(parent_hap, list):
		parent_hap = ord(parent_hap)
	elements = np.array(child_hap) == np.array(parent_hap)
	match = np.sum(elements)
	if isinstance(elements, np.bool_):
		elements = np.array([elements])
	total = len(elements)
	percentage = round(np.true_divide(match, total) * 100.0, 2) # Round with 2 decimals
	return percentage # Float with 2 decimals as percentage


def get_best_match(compare, child_hap, parent1_hap1, parent1_hap2, parent2_hap1, parent2_hap2, parent1_index, parent2_index):
	'''
	This function determines the best match for the child haplotype, using all haplotypes from parent(s).

	Input: The function compare one child haplotype, four (if present) haplotypes from the parents,
	parent1_index and parent2_index (where the parents are in order, if not present the index is None)

	Output: The best matching parent haplotype and the percentage for each tile, so a list with for each tile a best matching
	parent and a matching percentage
	'''
	outcome1 = outcome2 = outcome3 = outcome4 = 0
	if parent1_index != None: # If parent 1 is present
		outcome1 = compare(child_hap, parent1_hap1)
		outcome2 = compare(child_hap, parent1_hap2)
	elif parent1_index == None:
		outcome1 = outcome2 = 0
	if parent2_index != None: # If parent 2 is present
		outcome3 = compare(child_hap, parent2_hap1)
		outcome4 = compare(child_hap, parent2_hap2)
	elif parent2_index == None:
		outcome3 = outcome4 = 0
	outcomes = [outcome1, outcome2, outcome3, outcome4]
	best_outcome = max(outcomes) # Check which percentage from the four outcomes is the highest
	parent_list = ["father_hap1", "father_hap2", "mother_hap1", "mother_hap2"]
	if best_outcome:
		highest_perc = max(outcomes) # Get the max outcome in the list
		best_match_indices = [i for i, j in enumerate(outcomes) if j == highest_perc] # Select all indices with the max outcome
		best_match_parents = [parent_list[i] for i in best_match_indices]# Select with indices all parents with max outcome
	else:
		best_match_parents = None
	outcomes.append(best_match_parents)
	return outcomes # This is a list with for each tile the best matching parent haplotype and the percentage


# Helper function to make sure lists are the same length
def verify_same_length(*args):
	"""
	This function checks if the given arguments are of the same length
	Input: a string
	Output: TRUE if arguments are of the same length, FALSE if not
	"""
	if all(len(args[0]) == len(_arg) for _arg in args[1:]):
		return True
	return False


def make_tiles(tile_size_two_thirds, tile_size, *args):
	'''
	Input: a list with each item a nucleotide or a string with the nucleotides
	This function makes from the input tiles from tile size tile_size with an overlap of one third.
	Output: a list with lists, each item containing one nucleotide
	'''
	# Input is a string
	if not verify_same_length(args):
		raise ValueError("Hap lists/strings not of same length")

	global_pos = 0
	results = [[] for i in range(len(args))]
	while global_pos < len(args[0])-tile_size:
		for i,a in enumerate(args):
			results[i].append(list(a[int(round(global_pos)):int(round(global_pos+tile_size))]))
		global_pos += tile_size_two_thirds

	if global_pos is not len(args[0])-tile_size:
		for i,a in enumerate(args):
			results[i].append(list(a[int(round(global_pos)):]))
	return results # Results is a list containing lists




def calculate_other_parent(child, child_hap_number, make_tiles, compare, get_best_match, tile_size_two_thirds, tile_size, parent1_index, parent2_index, parent_list, child_hap1_complete, child_hap2_complete, parent1_hap1_complete, parent1_hap2_complete, parent2_hap1_complete, parent2_hap2_complete):
	'''
	This function is called when both child haplotypes are matching with the same parent.
	One haplotype of the child is inherited of one parent and the other child haplotype is
	inherited from the other parent.
	Input: The information to calculate the new percentages for each tile and the tiles from child and parent(s)
	Output: The new best parental haplotype and the matching percentage for this tile
	'''
	print("Calculating again for the other parent!")
	children_results2 = {}
	parent_counter2 = {}
	start_tile_line_nr = 0
	children_results2[child] = {1:[], 2:[]}
	parent_counter2[child] = {1:{"father_hap1":0, "father_hap2":0, "mother_hap1":0, "mother_hap2":0}, 2:{"father_hap1":0, "father_hap2":0, "mother_hap1":0, "mother_hap2":0}}

	child_hap1, child_hap2, parent1_hap1, parent1_hap2, parent2_hap1, parent2_hap2 = make_tiles(tile_size_two_thirds, tile_size, child_hap1_complete, child_hap2_complete, parent1_hap1_complete, parent1_hap2_complete, parent2_hap1_complete, parent2_hap2_complete) # Save tiles

	for (child_hap1_tile, child_hap2_tile, parent1_hap1_tile, parent1_hap2_tile, parent2_hap1_tile, parent2_hap2_tile) in zip(child_hap1, child_hap2, parent1_hap1, parent1_hap2, parent2_hap1, parent2_hap2):
		if child_hap_number == 1:
			best_match_result = get_best_match(compare, child_hap1_tile, parent1_hap1_tile, parent1_hap2_tile, parent2_hap1_tile, parent2_hap2_tile, parent1_index, parent2_index) # Return best parents and percentage
			children_results2[child][child_hap_number].append(best_match_result)
			if best_match_result[4]:
				for parent in best_match_result[4]:
					parent_counter2[child][child_hap_number][parent] += 1
		elif child_hap_number == 2:
			best_match_result = get_best_match(compare, child_hap2_tile, parent1_hap1_tile, parent1_hap2_tile, parent2_hap1_tile, parent2_hap2_tile, parent1_index, parent2_index) # Return best parents and percentage
			children_results2[child][child_hap_number].append(best_match_result)
			if best_match_result[4]:
				for parent in best_match_result[4]:
					parent_counter2[child][child_hap_number][parent] += 1
		else:
			print("Wrong child haplotype number, give 1 or 2")


	best_parent = max(parent_counter2[child][child_hap_number], key = parent_counter2[child][child_hap_number].get)
	best_parent_index = parent_list.index(best_parent)
	par = 0
	for tile_result in children_results2[child][child_hap_number]:
		par += tile_result[best_parent_index]
	mean_perc_best_parent = round(np.true_divide(par, len(children_results2[child][child_hap_number])), 2)
	return best_parent, mean_perc_best_parent # return the best parental haplotype with the matching percentage


##EDITED: make dict for haplo snps + remove conversion from 0/1 to nucleotides
def process_haplo_file(individuals, children, haplo_file, tile_size, tile_size_two_thirds):
	'''
	This function determines the best parent and the matching percentage for each child in the dataset for both child haplotypes.
	If the best parent is the same for both child haplotypes, then the matching percentage for the other parent (if present)
	is calculated for the lowest percentage haplotype. It also collects the complete haplotypes for the child and parent(s)
	It prints the results to the commandline for each child.

	Input: The individuals dict (containing each individual in the dataset and their index in the .sample file),
	the children dict (containing all children from the dataset and their parent(s)), the .haps file,
	the given tile_size and the tile_size_two_thirds to let the thiles have overlap

	Output: The best parents in a list per child for haplotype 1 (best_parent1_list) and haplotype 2 (best_parent2_list),
	the complete haplotypes for the child and their parent(s) in a list.
	'''
	children_results = {} # Make dict for lists of results for all children
	child_nr = 1
	best_parent1_list = []
	best_parent2_list = []
	perc_hap1_list = []
	perc_hap2_list = []
	parent1_index_list = []
	parent2_index_list = []
	haplo_file_snps = []
    
	with open(haplo_file) as haps:
		haps_lines = [line.rstrip() for line in haps]
 #get haplo file info (SNPs)
		for line in haps_lines:
			snps = line.strip().split(' ')
			tmp_dict={}
			tmp_dict["CHR"] = snps[0]
			tmp_dict["SNP"] = snps[1]
			tmp_dict["LOC"] = snps[2]
			tmp_dict["A1"] = snps[3]
			tmp_dict["A2"] = snps[4]
			haplo_file_snps.append(tmp_dict)
        
 #resume script 
		parent1_hap1_list, parent1_hap2_list, parent2_hap1_list, parent2_hap2_list, child_hap1_list, child_hap2_list = [], [], [], [], [], []
		for child in children:
			print("\n--------\nchild " + str(child_nr) + "\n--------") # Print the number of the child
			parent_counter = {} # Make dict for counting which parent was most often the best in a tile
			children_results[child] = {1:[], 2:[]}
			parent_counter[child] = {1:{"father_hap1":0, "father_hap2":0, "mother_hap1":0, "mother_hap2":0}, 2:{"father_hap1":0, "father_hap2":0, "mother_hap1":0, "mother_hap2":0}}
			#plot_list_hap1 = []# Per tile a list with matching percentages
			#plot_list_hap2 = []

			########## Collect complete haplotypes for parents and child in this chunk #########
			parent1_hap1_complete, parent1_hap2_complete, parent2_hap1_complete, parent2_hap2_complete, child_hap1_complete, child_hap2_complete = [], [], [], [], [], []
			parent1_index_haplotype = None
			if children[child][0] != "0" and children[child][0] in individuals:
				parent1_index_haplotype = 2 * (individuals[children[child][0]])
			parent2_index_haplotype = None
			if children[child][1] != "0" and children[child][1] in individuals:
				parent2_index_haplotype = 2 * (individuals[children[child][1]])
			child_index_haplotype = 2 * (individuals[child])

##EDITED: instead of converting to nucleotides, keep 0s and 1s
			for line in haps_lines:
				line = line.strip().split(" ")
				line = line[3:]
				for index, item in enumerate(line):
					if item == "0": # Replace 0 with first variant
						line[index] = "0"
					elif item == "1": # Replace 1 with the second variant
						line[index] = "1"
				line = line[2:]

				if parent1_index_haplotype != None:
					parent1_hap1_complete.append(line[parent1_index_haplotype])
					parent1_hap2_complete.append(line[parent1_index_haplotype + 1])

				if parent2_index_haplotype != None:
					parent2_hap1_complete.append(line[parent2_index_haplotype])
					parent2_hap2_complete.append(line[parent2_index_haplotype + 1])

				child_hap1_complete.append(line[child_index_haplotype])
				child_hap2_complete.append(line[child_index_haplotype + 1])

##EDITED: instead of '0', 'N' for missing

			if parent1_index_haplotype == None:
				parent1_hap1_list.append("N")
				parent1_hap2_list.append("N")
			else:
				parent1_hap1_list.append(parent1_hap1_complete)
				parent1_hap2_list.append(parent1_hap2_complete)

			if parent2_index_haplotype == None:
				parent2_hap1_list.append("N")
				parent2_hap2_list.append("N")
			else:
				parent2_hap1_list.append(parent2_hap1_complete)
				parent2_hap2_list.append(parent2_hap2_complete)

			child_hap1_list.append(child_hap1_complete)
			child_hap2_list.append(child_hap2_complete)
			################## end ##############################################################
			start_tile_line_nr = 0
			parent1_index = None
			if children[child][0] != "0" and children[child][0] in individuals:
				parent1_index = 5+2*(individuals[children[child][0]])
			parent2_index = None
			if children[child][1] != "0" and children[child][1] in individuals:
				parent2_index = 5+2*(individuals[children[child][1]])
			child_index = 5+2*(individuals[child]) # Save child index

			if parent1_index == None:
				print("No father in the data for this child")
			if parent2_index == None:
				print("No mother in the data for this child")

			# Make the tiles list for each person:
			child_hap1, child_hap2, parent1_hap1, parent1_hap2, parent2_hap1, parent2_hap2 = make_tiles(tile_size_two_thirds, tile_size, child_hap1_complete, child_hap2_complete, parent1_hap1_complete, parent1_hap2_complete, parent2_hap1_complete, parent2_hap2_complete)
			# Loop through the tiles
			tile_nr_count=1
			for (child_hap1_tile, child_hap2_tile, parent1_hap1_tile, parent1_hap2_tile, parent2_hap1_tile, parent2_hap2_tile) in zip(child_hap1, child_hap2, parent1_hap1, parent1_hap2, parent2_hap1, parent2_hap2):
				#plot_list_per_tile_hap1 = []
				#plot_list_per_tile_hap2 = []
				# Child hap 1
				best_match_result = get_best_match(compare, child_hap1_tile, parent1_hap1_tile, parent1_hap2_tile, parent2_hap1_tile, parent2_hap2_tile, parent1_index, parent2_index)
				children_results[child][1].append(best_match_result) # Add to dict
				if best_match_result[4]:
					for parent in best_match_result[4]:
						parent_counter[child][1][parent] += 1

				#plot_list_per_tile_hap1.append(tile_nr_count)
				#plot_list_per_tile_hap1.append(best_match_result[0])
				#plot_list_per_tile_hap1.append(best_match_result[1])
				#plot_list_per_tile_hap1.append(best_match_result[2])
				#plot_list_per_tile_hap1.append(best_match_result[3])
				# Child hap 2
				best_match_result = get_best_match(compare, child_hap2_tile, parent1_hap1_tile, parent1_hap2_tile, parent2_hap1_tile, parent2_hap2_tile, parent1_index, parent2_index)
				children_results[child][2].append(best_match_result)
				if best_match_result[4]:
					for parent in best_match_result[4]:
						parent_counter[child][2][parent] += 1
				#plot_list_per_tile_hap2.append(tile_nr_count)
				#plot_list_per_tile_hap2.append(best_match_result[0])
				#plot_list_per_tile_hap2.append(best_match_result[1])
				#plot_list_per_tile_hap2.append(best_match_result[2])
				#plot_list_per_tile_hap2.append(best_match_result[3])

				#plot_list_hap1.append(plot_list_per_tile_hap1) # per child the information for the plots
				#plot_list_hap2.append(plot_list_per_tile_hap2) # per child the information for the plots
				tile_nr_count += 1

			# Write per child a file with the name {childID}{_hapnumber}.csv
			#filename_hap1 = "{}{}.csv".format(child,"_hap1")
			#filename_hap2 = "{}{}.csv".format(child, "_hap2")

			#write_list_to_csv(plot_list_hap1, filename_hap1)
			#write_list_to_csv(plot_list_hap2, filename_hap2)

			parent_list = ["father_hap1", "father_hap2", "mother_hap1", "mother_hap2"]
			best_parent_1 = max(parent_counter[child][1], key=parent_counter[child][1].get)
			best_parent_2 = max(parent_counter[child][2], key=parent_counter[child][2].get)
			best_parent_1_index = parent_list.index(best_parent_1)
			best_parent_2_index = parent_list.index(best_parent_2)
			par = 0
			for tile_result in children_results[child][1]:
					par += tile_result[best_parent_1_index]
			mean_perc_best_parent_1 = round(np.true_divide(par, len(children_results[child][1])), 2)

			par = 0
			for tile_result in children_results[child][2]:
				par += tile_result[best_parent_2_index]
			mean_perc_best_parent_2 = round(np.true_divide(par, len(children_results[child][2])), 2)

			############# In this code chunk check if the same parent is the best one for both child haplotypes, if yes calculate the percentage and best matching haplotype from the other parent (if present) ######
			if re.match("father_hap[1-2]", best_parent_1) and re.match("father_hap[1-2]", best_parent_2):
				print("Both haplotypes matches with father!")
				if mean_perc_best_parent_1 < mean_perc_best_parent_2:
					#Check if other parent is present, in this case parent 2: mom and for child hap 1
					if parent2_index != None:
						# Use function to calculate for other parent (in this case mother)
						best_parent_1, mean_perc_best_parent_1 = calculate_other_parent(child, 1, make_tiles, compare, get_best_match, tile_size_two_thirds, tile_size, None, parent2_index, parent_list, child_hap1_complete, child_hap2_complete, parent1_hap1_complete, parent1_hap2_complete, parent2_hap1_complete, parent2_hap2_complete)
					else:
						best_parent_1 = "Mother_hap"
						mean_perc_best_parent_1 = None

				elif mean_perc_best_parent_1 > mean_perc_best_parent_2:
					#Check if other parent is present, in this case parent 2: mom and this is for child hap 2
					if parent2_index != None:
						best_parent_2, mean_perc_best_parent_2 = calculate_other_parent(child, 2, make_tiles, compare, get_best_match, tile_size_two_thirds, tile_size, None, parent2_index, parent_list, child_hap1_complete, child_hap2_complete, parent1_hap1_complete, parent1_hap2_complete, parent2_hap1_complete, parent2_hap2_complete)
					else:
						best_parent_2 = "Mother_hap"
						mean_perc_best_parent_2 = None
				else:
					print("Percentages are the same")

			if re.match("mother_hap[1-2]", best_parent_1) and re.match("mother_hap[1-2]", best_parent_2):
				print("Both haplotypes matches with mom!")
				if mean_perc_best_parent_1 < mean_perc_best_parent_2:
					if parent1_index != None:
						best_parent_1, mean_perc_best_parent_1 = calculate_other_parent(child, 1, make_tiles, compare, get_best_match, tile_size_two_thirds, tile_size, parent1_index, None, parent_list, child_hap1_complete, child_hap2_complete, parent1_hap1_complete, parent1_hap2_complete, parent2_hap1_complete, parent2_hap2_complete)
					else:
						best_parent_1 = "Father_hap"
						mean_perc_best_parent_1 = None

				elif mean_perc_best_parent_1 > mean_perc_best_parent_2:
					if parent1_index != None:
						best_parent_2, mean_perc_best_parent_2 = calculate_other_parent(child, 2, make_tiles, compare, get_best_match, tile_size_two_thirds, tile_size, parent1_index, None, parent_list, child_hap1_complete, child_hap2_complete, parent1_hap1_complete, parent1_hap2_complete, parent2_hap1_complete, parent2_hap2_complete)
					else:
						best_parent_2 = "Father_hap"
						mean_perc_best_parent_2 = None

				else:
					print("Percentages are the same")
					########### End of code chunk to check if same parent is the best parent for both child haplotypes ################################################
			print("The best matching parent for haplotype 1 is: " + str(best_parent_1) + " with a matching percentage of " + str(mean_perc_best_parent_1) + "%")
			print("The best matching parent for haplotype 2 is: " + str(best_parent_2) + " with a matching percentage of " + str(mean_perc_best_parent_2) + "%")
			best_parent1_list.append(best_parent_1)
			best_parent2_list.append(best_parent_2)
			perc_hap1_list.append(mean_perc_best_parent_1)
			perc_hap2_list.append(mean_perc_best_parent_2)
			parent1_index_list.append(parent1_index)
			parent2_index_list.append(parent2_index)
			child_nr += 1

	parent1_hap1_list = list(map("".join, parent1_hap1_list))
	parent1_hap2_list = list(map("".join, parent1_hap2_list))
	parent2_hap1_list = list(map("".join, parent2_hap1_list))
	parent2_hap2_list = list(map("".join, parent2_hap2_list))
	child_hap1_list = list(map("".join, child_hap1_list))
	child_hap2_list = list(map("".join, child_hap2_list))
	return best_parent1_list, best_parent2_list, parent1_hap1_list, parent1_hap2_list, parent2_hap1_list, parent2_hap2_list, child_hap1_list, child_hap2_list, perc_hap1_list, perc_hap2_list, parent1_index_list, parent2_index_list, haplo_file_snps


def replacer(s, newstring, index, nofail=False):
	"""
	This function is to replace part of a string with another string
	Input: s as original string, a new string (the replacer string) and the index
	on which the string need to be replaced.
	Output: The new string with the replaced part
	"""
	# raise an error if index is outside of the string
	if not nofail and index not in range(len(s)):
		print("Length",len(s))
		print("Index",index)
		raise ValueError("index outside given string")

	# if not erroring, but the index is still not in the correct range..
	if index < 0:  # add it to the beginning
		return newstring + s
	if index > len(s):  # add it to the end
		return s + newstring

	# insert the new string between "slices" of the original
	return s[:index] + newstring + s[index + len(newstring):]

def convert_to_list(string_to_convert):
	"""This function is to convert a string to a list
	Input: string to convert
	Output: list from the string
	"""
	list1=[]
	list1[:0]=string_to_convert
	return list1

def get_most_frequent(list_to_search):
	"""
	This function is to find the most frequent item in a list
	Input: the list
	Output: Most frequent item in the list
	"""
	if not list_to_search:
		return None
	elif None in list_to_search:
		raise ValueError("List contains None value that should not be there")
	mylist = np.array(list_to_search)
	(unique, counts) = np.unique(mylist, return_counts=True)
	return unique[np.where(counts == np.amax(counts))][0]


def implement_recombination(children, tile_size, make_tiles, child_hap1_list, child_hap2_list, perc_hap1_list, perc_hap2_list, get_best_match, best_parent1_list, best_parent2_list, replacer, tile_size_two_thirds, parent1_hap1_list, parent1_hap2_list, parent2_hap1_list, parent2_hap2_list, parent1_index_list, parent2_index_list):
	"""
	This function is to find the recombination spots/tiles and recombination parental haplotype.
	Then implement recombination for these tiles with the recombination parental haplotype.

	Input: The dictionary with all children, tile size, the child haplotypes,
	the overall best matching parental haplotypes (best_parent_list) and their percentages
	(perc_hap_list) and the replacer, make_tiles and get_best_match function

	Output: The transmitted alleles for child hap 1, child hap 2 and the non-transmitted alleles
	(with implemented recombination) for child hap 1 and 2.
	"""
	transmitted1 = [] # Transmitted haplotypes for child haplotype 1
	transmitted2 = [] # Transmitted haplotypes for child haplotype 2
	non_transmitted_1 = []
	non_transmitted_2 = []
	tile_size_recombination = tile_size # Use same tile size for recombination
	tile_size_two_thirds_recombination = (float(2)/float(3)) * tile_size_recombination

	for index, childID in enumerate(children):
		# Make tiles for this child:
		child_hap1, child_hap2, parent1_hap1, parent1_hap2, parent2_hap1, parent2_hap2 = make_tiles(tile_size_two_thirds, tile_size, child_hap1_list[index], child_hap2_list[index], parent1_hap1_list[index], parent1_hap2_list[index], parent2_hap1_list[index], parent2_hap2_list[index])
		print("\nchildID", index +1)
		parent_counter = {} # for each child the parent counter 1 and 2
		parent_counter = {1:{"father_hap1":0, "father_hap2":0, "mother_hap1":0, "mother_hap2":0}, 2:{"father_hap1":0, "father_hap2":0, "mother_hap1":0, "mother_hap2":0}}
		if child_hap1_list[index] != None: # Check if allele 1 is present, normally it is
			transmitted1.append(child_hap1_list[index])
		else:
			transmitted1.append("0")

		if child_hap2_list[index] != None: # Check if allele 2 is present, normally it is
			transmitted2.append(child_hap2_list[index])
		else:
			transmitted2.append("0")
		# Compare original child haps with parental haps with tiles to see if some parts are matching better
		# with another parent: recombination. Determine with which parent hap the recombination happens and the locations.
		run = 0 # To count the number of tiles
		if perc_hap1_list[index] != None:
			if child_hap1_list[index] != None and perc_hap1_list[index] < 99.9: # First check if allele 1 for child is present and if there is recombination
				tile_numbers_recombination1 = []
				for (child_hap1_tile, child_hap2_tile, parent1_hap1_tile, parent1_hap2_tile, parent2_hap1_tile, parent2_hap2_tile) in zip(child_hap1, child_hap2, parent1_hap1, parent1_hap2, parent2_hap1, parent2_hap2):
					# Compare child hap 1 with all parents to get best_parent1_list[index] and the recombination parent if present
					best_parent_recombination_output = get_best_match(compare, child_hap1_tile, parent1_hap1_tile, parent1_hap2_tile, parent2_hap1_tile, parent2_hap2_tile, parent1_index_list[index], parent2_index_list[index])
					if best_parent_recombination_output[4]: # check if list is not empty
						for parent in best_parent_recombination_output[4]: # for each chosen best parent in this tile
							if parent != best_parent1_list[index]: # check if parent is not the best parent overall for this haps
								parent_counter[1][parent] += 1
								tile_numbers_recombination1.append(run)
					run += 1
			elif perc_hap1_list[index] >= 99.9:
				print("no recombination in this hap")
		else:
			print("The best parent for this child hap is not present")

		if sum(parent_counter.get(1).values()) != 0:
			recombination_parent1 = max(parent_counter[1], key=parent_counter[1].get)
			# get only unique values from tile_numbers_recombination1, because if there are more than one parents in the list, the location will appear multiple times
			tile_numbers_recombination1 = list(set(tile_numbers_recombination1))
		else:
			recombination_parent1 = None
			tile_numbers_recombination1 = None

		haps2 = list(zip(child_hap2_list[index], parent1_hap1_list[index], parent1_hap2_list[index], parent2_hap1_list[index], parent2_hap2_list[index]))
		run2 = 0 # To count the number of tiles
		if perc_hap2_list[index] != None:
			if child_hap2_list[index] != None and perc_hap2_list[index] < 99.9:
				print(perc_hap2_list[index])
				tile_numbers_recombination2 = []
				for (child_hap1_tile, child_hap2_tile, parent1_hap1_tile, parent1_hap2_tile, parent2_hap1_tile, parent2_hap2_tile) in zip(child_hap1, child_hap2, parent1_hap1, parent1_hap2, parent2_hap1, parent2_hap2):
					# Compare child hap 1 with all parents to get best_parent1_list[index] and the recombination parent if present
					best_parent_recombination_output = get_best_match(compare, child_hap2_tile, parent1_hap1_tile, parent1_hap2_tile, parent2_hap1_tile, parent2_hap2_tile, parent1_index_list[index], parent2_index_list[index])
					if best_parent_recombination_output[4]: # Check if list is not empty
						for parent in best_parent_recombination_output[4]:
							if parent != best_parent2_list[index]:
								parent_counter[2][parent] += 1
								tile_numbers_recombination2.append(run2)
					run2 += 1
			elif perc_hap2_list[index] >= 99.9:
				print("no recombination in this hap")
		else:
			print("The best parent for this child hap is not present")

		if sum(parent_counter.get(2).values()) != 0:
			recombination_parent2 = max(parent_counter[2], key=parent_counter[2].get)
			# Get unique values from tile_numbers_recombination2
			tile_numbers_recombination2 = list(set(tile_numbers_recombination2))
		else:
			recombination_parent2 = None
			tile_numbers_recombination2 = None

		non_transmitted_haps1_tmp = ""
		if best_parent1_list[index] != None: # Check if best parent for hap 1 is present
			if best_parent1_list[index] == "father_hap1": # This is the best parent hap for transmitted, so use for non-transmitted the other hap from this parent
				if not recombination_parent1: # No recombination parents
					non_transmitted_haps1_tmp = parent1_hap2_list[index]
				else:
					non_transmitted_haps1_tmp = parent1_hap2_list[index]
					for index_recombination in tile_numbers_recombination1:
						index_recombination_start = int(round(index_recombination * tile_size_two_thirds_recombination))
						index_recombination_end = int(round((index_recombination * tile_size_two_thirds_recombination) + tile_size_recombination))
						# These recombination parents are based on transmitted, so the other one need to be used
						# Replace the opposite of best_parent in this case because of this.
						if recombination_parent1 == "father_hap2":
							non_transmitted_haps1_tmp = replacer(non_transmitted_haps1_tmp, parent1_hap1_list[index][index_recombination_start:index_recombination_end], index_recombination_start)
						elif recombination_parent1 == "mother_hap1":
							non_transmitted_haps1_tmp = replacer(non_transmitted_haps1_tmp, parent2_hap2_list[index][index_recombination_start:index_recombination_end], index_recombination_start)
						elif recombination_parent1 == "mother_hap2":
							non_transmitted_haps1_tmp = replacer(non_transmitted_haps1_tmp, parent2_hap1_list[index][index_recombination_start:index_recombination_end], index_recombination_start)

			elif best_parent1_list[index] == "father_hap2": # This is the best parent hap for transmitted, so use for non-transmitted the other hap from this parent
				if not recombination_parent1: # No recombination parents
					non_transmitted_haps1_tmp = parent1_hap1_list[index]
				else:
					non_transmitted_haps1_tmp = parent1_hap1_list[index]
					for index_recombination in tile_numbers_recombination1:
						index_recombination_start = int(round(index_recombination * tile_size_two_thirds_recombination))
						index_recombination_end = int(round((index_recombination * tile_size_two_thirds_recombination) + tile_size_recombination))
						if recombination_parent1 == "father_hap1":
							non_transmitted_haps1_tmp = replacer(non_transmitted_haps1_tmp, parent1_hap2_list[index][index_recombination_start:index_recombination_end], index_recombination_start)
						elif recombination_parent1 == "mother_hap1":
							non_transmitted_haps1_tmp = replacer(non_transmitted_haps1_tmp, parent2_hap2_list[index][index_recombination_start:index_recombination_end], index_recombination_start)
						elif recombination_parent1 == "mother_hap2":
							non_transmitted_haps1_tmp = replacer(non_transmitted_haps1_tmp, parent2_hap1_list[index][index_recombination_start:index_recombination_end], index_recombination_start)

			elif best_parent1_list[index] == "mother_hap1": # This is the best parent hap for transmitted, so use for non-transmitted the other hap from this parent
				if not recombination_parent1: # No recombination parents
					non_transmitted_haps1_tmp = parent2_hap2_list[index]
				else:
					non_transmitted_haps1_tmp = parent2_hap2_list[index]
					for index_recombination in tile_numbers_recombination1:
						index_recombination_start = int(round(index_recombination * tile_size_two_thirds_recombination))
						index_recombination_end = int(round((index_recombination * tile_size_two_thirds_recombination) + tile_size_recombination))
						if recombination_parent1 == "father_hap1":
							non_transmitted_haps1_tmp = replacer(non_transmitted_haps1_tmp, parent1_hap2_list[index][index_recombination_start:index_recombination_end], index_recombination_start)
						elif recombination_parent1 == "father_hap2":
							non_transmitted_haps1_tmp = replacer(non_transmitted_haps1_tmp, parent1_hap1_list[index][index_recombination_start:index_recombination_end], index_recombination_start)
						elif recombination_parent1 == "mother_hap2":
							non_transmitted_haps1_tmp = replacer(non_transmitted_haps1_tmp, parent2_hap1_list[index][index_recombination_start:index_recombination_end], index_recombination_start)

			elif best_parent1_list[index] == "mother_hap2": # This is the best parent hap for transmitted, so use for non-transmitted the other hap from this parent
				if not recombination_parent1: # No recombination parents
					non_transmitted_haps1_tmp = parent2_hap1_list[index]
					print("no recombination")
				else:
					non_transmitted_haps1_tmp = parent2_hap1_list[index]
					for index_recombination in tile_numbers_recombination1:
						index_recombination_start = int(round(index_recombination * tile_size_two_thirds_recombination))
						index_recombination_end = int(round((index_recombination * tile_size_two_thirds_recombination) + tile_size_recombination))
						if recombination_parent1 == "father_hap1":
							non_transmitted_haps1_tmp = replacer(non_transmitted_haps1_tmp, parent1_hap2_list[index][index_recombination_start:index_recombination_end], index_recombination_start)
						elif recombination_parent1 == "father_hap2":
							non_transmitted_haps1_tmp = replacer(non_transmitted_haps1_tmp, parent1_hap1_list[index][index_recombination_start:index_recombination_end], index_recombination_start)
						elif recombination_parent1 == "mother_hap1":
							non_transmitted_haps1_tmp = replacer(non_transmitted_haps1_tmp, parent2_hap2_list[index][index_recombination_start:index_recombination_end], index_recombination_start)
			else:
				non_transmitted_haps1_tmp = "0"
		else:
			non_transmitted_haps1_tmp = "0" # Add string with 0 of length parent hap that is present if child hap 1 is not present
		non_transmitted_1.append(non_transmitted_haps1_tmp) # Add for this child the non-transmitted hap to the list

		non_transmitted_haps2_tmp = ""
		#Fill non-transmitted hap 2
		if best_parent1_list[index] != None:
			if best_parent2_list[index] == "father_hap1":
				if not recombination_parent2:
					non_transmitted_haps2_tmp = parent1_hap2_list[index]
				else:
					non_transmitted_haps2_tmp = parent1_hap2_list[index]
					for index_recombination in tile_numbers_recombination2:
						index_recombination_start = int(round(index_recombination * tile_size_two_thirds_recombination))
						index_recombination_end = int(round((index_recombination * tile_size_two_thirds_recombination) + tile_size_recombination))
						if recombination_parent2 == "father_hap2":
							non_transmitted_haps2_tmp = replacer(non_transmitted_haps2_tmp, parent1_hap1_list[index][index_recombination_start:index_recombination_end], index_recombination_start)
						elif recombination_parent2 == "mother_hap1":
							non_transmitted_haps2_tmp = replacer(non_transmitted_haps2_tmp, parent2_hap2_list[index][index_recombination_start:index_recombination_end], index_recombination_start)
						elif recombination_parent2 == "mother_hap2":
							non_transmitted_haps2_tmp = replacer(non_transmitted_haps2_tmp, parent2_hap1_list[index][index_recombination_start:index_recombination_end], index_recombination_start)

			elif best_parent2_list[index] == "father_hap2":
				if not recombination_parent2:
					non_transmitted_haps2_tmp = parent1_hap1_list[index]
				else:
					non_transmitted_haps2_tmp = parent1_hap1_list[index]
					for index_recombination in tile_numbers_recombination2:
						index_recombination_start = int(round(index_recombination * tile_size_two_thirds_recombination))
						index_recombination_end = int(round((index_recombination * tile_size_two_thirds_recombination) + tile_size_recombination))
						if recombination_parent2 == "father_hap1":
							non_transmitted_haps2_tmp = replacer(non_transmitted_haps2_tmp, parent1_hap2_list[index][index_recombination_start:index_recombination_end], index_recombination_start)
						elif recombination_parent2 == "mother_hap1":
							non_transmitted_haps2_tmp = replacer(non_transmitted_haps2_tmp, parent2_hap2_list[index][index_recombination_start:index_recombination_end], index_recombination_start)
						elif recombination_parent2 == "mother_hap2":
							non_transmitted_haps2_tmp = replacer(non_transmitted_haps2_tmp, parent2_hap1_list[index][index_recombination_start:index_recombination_end], index_recombination_start)

			elif best_parent2_list[index] == "mother_hap1":
				if not recombination_parent2:
					non_transmitted_haps2_tmp = parent2_hap2_list[index]
				else:
					non_transmitted_haps2_tmp = parent2_hap2_list[index]
					for index_recombination in tile_numbers_recombination2:
						index_recombination_start = int(round(index_recombination * tile_size_two_thirds_recombination))
						index_recombination_end = int(round((index_recombination * tile_size_two_thirds_recombination) + tile_size_recombination))
						if recombination_parent2 == "father_hap1":
							non_transmitted_haps2_tmp = replacer(non_transmitted_haps2_tmp, parent1_hap2_list[index][index_recombination_start:index_recombination_end], index_recombination_start)
						elif recombination_parent2 == "father_hap2":
							non_transmitted_haps2_tmp = replacer(non_transmitted_haps2_tmp, parent1_hap1_list[index][index_recombination_start:index_recombination_end], index_recombination_start)
						elif recombination_parent2 == "mother_hap2":
							non_transmitted_haps2_tmp = replacer(non_transmitted_haps2_tmp, parent2_hap1_list[index][index_recombination_start:index_recombination_end], index_recombination_start)

			elif best_parent2_list[index] == "mother_hap2":
				if not recombination_parent2:
					non_transmitted_haps2_tmp = parent2_hap1_list[index]
				else:
					non_transmitted_haps2_tmp = parent2_hap1_list[index]
					for index_recombination in tile_numbers_recombination2:
						index_recombination_start = int(round(index_recombination * tile_size_two_thirds_recombination))
						index_recombination_end = int(round((index_recombination * tile_size_two_thirds_recombination) + tile_size_recombination))
						if recombination_parent2 == "father_hap1":
							non_transmitted_haps2_tmp = replacer(non_transmitted_haps2_tmp, parent1_hap2_list[index][index_recombination_start:index_recombination_end], index_recombination_start)
						elif recombination_parent2 == "father_hap2":
							non_transmitted_haps2_tmp = replacer(non_transmitted_haps2_tmp, parent1_hap1_list[index][index_recombination_start:index_recombination_end], index_recombination_start)
						elif recombination_parent2 == "mother_hap1":
							non_transmitted_haps2_tmp = replacer(non_transmitted_haps2_tmp, parent2_hap2_list[index][index_recombination_start:index_recombination_end], index_recombination_start)
			else:
				non_transmitted_haps2_tmp = "0"
		else:
			non_transmitted_haps2_tmp = "0" # Add string with 0 of length parent hap that is present if child hap 2 is not present
		non_transmitted_2.append(non_transmitted_haps2_tmp) # Add for this child the non-transmitted hap to the list
	return transmitted1, transmitted2, non_transmitted_1, non_transmitted_2

#new output functions: write haps files    

def write_transmitted_haps(transmitted1, transmitted2,haplo_file_snps,output_transmitted_haps):
    print("Writing output transmitted .haps files..")
    
    transmitted1_sep = [[x for x in y] for y in transmitted1]
    transmitted2_sep = [[x for x in y] for y in transmitted2]
    
    res2 = []
    for idx in range(len(transmitted1_sep)):
        res2.append(transmitted1_sep[idx])
        res2.append(transmitted2_sep[idx])
            
    
    transmitted_haplos = pd.DataFrame(res2).T
    snps = pd.DataFrame.from_dict(haplo_file_snps) 
    concat = pd.concat([snps,transmitted_haplos],axis=1)
    concat.to_csv(output_transmitted_haps,sep=" ",index=False,header=False)
    

def write_nontransmitted_haps(non_transmitted_1, non_transmitted_2,haplo_file_snps,output_nontransmitted_haps):
    print("Writing output non-transmitted .haps files..")
    
    #if the list in non_transmitted2_sep is = ['0'], replace it with
    #list with the length of the other lists with "N"s or whatever else
    #for NA.
    
    for i in range(len(non_transmitted_1)):
        if non_transmitted_1[i] == "0":
            non_transmitted_1[i] = "N" * len(haplo_file_snps)
        if non_transmitted_2[i] == "0":
            non_transmitted_2[i] = "N" * len(haplo_file_snps)
    
    non_transmitted1_sep = [[x for x in y] for y in non_transmitted_1]
    non_transmitted2_sep = [[x for x in y] for y in non_transmitted_2]
    

    res2 = []
    for idx in range(len(non_transmitted1_sep)):
        res2.append(non_transmitted1_sep[idx])
        res2.append(non_transmitted2_sep[idx])
            
    non_transmitted_haplos = pd.DataFrame(res2).T
    snps = pd.DataFrame.from_dict(haplo_file_snps) 
    concat = pd.concat([snps,non_transmitted_haplos],axis=1)
    concat.to_csv(output_nontransmitted_haps,sep=" ",index=False,header=False)
     
    
#write sample files

def write_sample_files(sample_file_info,output_transmitted_sample,output_nontransmitted_sample):
	'''
	This function is used to write two identical .sample files 
    with only the children information 
    (this is for convenience of having two files with names corresponding to the .haps) 
	Input: sample_file_info, which is a list with dicts (for each line a dict)
	Output: Two .sample files with information of the children
	'''
	print("Writing output sample files..")
    
	sampleinfo = pd.DataFrame.from_dict(sample_file_info)
    
    #add "2nd header"
	header2 = pd.DataFrame(data=['0', '0', '0', 'D', 'D', 'D', 'B']).T
	header2.columns = sampleinfo.columns
	samplefile = pd.concat([header2,sampleinfo])
    
	samplefile.to_csv(output_transmitted_sample,sep=" ",index=False)
	samplefile.to_csv(output_nontransmitted_sample,sep=" ",index=False)   
    

def main():
	'''
	This is the main function
	'''
	haplo_file, sample_file, tile_size, output_transmitted_haps, output_transmitted_sample, output_nontransmitted_haps, output_nontransmitted_sample = argparser()
	individuals, children, sample_file_info = process_sample_file(sample_file)
	tile_size_two_thirds = (float(2)/float(3)) * tile_size # output: 2/3 from tile_size, overlap of 1/3 in tiles
	best_parent1_list, best_parent2_list, parent1_hap1_list, parent1_hap2_list, parent2_hap1_list, parent2_hap2_list, child_hap1_list, child_hap2_list, perc_hap1_list, perc_hap2_list, parent1_index_list, parent2_index_list,haplo_file_snps = process_haplo_file(individuals, children, haplo_file, tile_size, tile_size_two_thirds)

	transmitted1, transmitted2, non_transmitted_1, non_transmitted_2 = implement_recombination(children, tile_size, make_tiles, child_hap1_list, child_hap2_list, perc_hap1_list, perc_hap2_list, get_best_match, best_parent1_list, best_parent2_list, replacer, tile_size_two_thirds, parent1_hap1_list, parent1_hap2_list, parent2_hap1_list, parent2_hap2_list, parent1_index_list, parent2_index_list)

	write_transmitted_haps(transmitted1,transmitted2,haplo_file_snps,output_transmitted_haps)
	write_nontransmitted_haps(non_transmitted_1,non_transmitted_2,haplo_file_snps,output_nontransmitted_haps)
	write_sample_files(sample_file_info,output_transmitted_sample,output_nontransmitted_sample)

if __name__ == "__main__":
	main()


