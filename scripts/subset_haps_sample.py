"""
This script filters phased data (haps/sample) in chunks of N=1000

Input: Uphased .haps and .sample files
Output: .haps and .sample datasets with around N=1000 offspring + parents

Arguments: chr + path to original .haps file + path to original sample file + path to output directory

"""

# imports
import argparse
import pandas as pd
import numpy as np
import csv
import os

def argparser():
	"""
	Function to parse the command line arguments: chr + .haps file and .sample file + output directory
	"""
	parser = argparse.ArgumentParser(description="Chunking")
	parser.add_argument("chromosome", type=str, help="Chromosome")
	parser.add_argument("haps_file", type=str, help="Path to .haps file")
	parser.add_argument("sample_file", type=str, help="Path to .sample file")
	parser.add_argument("data_dir", type=str, help="Path to folder where files should be stored")
	parser.add_argument("--chunk-size", type=int, default=1000, help="Number of offspring per chunk (default: 1000)")
	args = parser.parse_args()

    return args.chromosome, args.haps_file, args.sample_file, args.data_dir, args.chunk_size

def sample_chunker(sample_file, data_dir):
    """
	This function filters and creates the trios .sample files

    In: original .sample file
    Out: one .sample file with trios
         list with indices of the IIDs so that .haps file can be filtered
	"""

    # open sample file as df
    sample_UGLI = pd.read_csv(sample_file, sep = " ")

    # save and drop second header
    header = (sample_UGLI.iloc[0]).to_frame().T
    sample_UGLI = sample_UGLI.drop([0])
    
    # whole df to str and reset idx
    sample_UGLI = sample_UGLI.astype(str).reset_index(drop=True)

        # condition: father OR mother available
    offspring  = sample_UGLI[((((sample_UGLI['father'] != 0) & (sample_UGLI['father'].isin(sample_UGLI['ID_2'])))
             | ((sample_UGLI['mother'] != 0) & (sample_UGLI['mother'].isin(sample_UGLI['ID_2'])))))]

        # split offspring into chunks of 1000
    chunks = list()
    chunk_size = 1000
    num_chunks = len(offspring) // chunk_size + 1
    for i in range(num_chunks):
        chunks.append(offspring[i*chunk_size:(i+1)*chunk_size])
        # add parents in each dataset and remove duplicates
    for i, dataset in enumerate(chunks):
        chunks[i] = pd.concat([dataset, sample_UGLI[sample_UGLI['ID_2'].isin(np.append(dataset['father'], dataset['mother']))]]).drop_duplicates()

        # save list of indices and create .sample files
    list_idx = []
    for i, dataset in enumerate(chunks):
        list_idx.append(list(dataset.index)) 
        chunks[i] = pd.concat([header,dataset]) # include second header
        filename= "offspring.{}.sample".format(i)
        filepath = os.path.join(data_dir, filename)
        chunks[i].to_csv(filepath,sep=" ",index=False)

    return list_idx

def filter_haps(haps_file, idx, haps_file_out):
    """
    This function filters the .haps file by including only information of
    the individuals present in the .sample file
    """

    haps_file_info = []
    with open(haps_file) as haps:
            for line in haps:
                haps1, haps2 = [], []
                tmp_dict = {}
                line = line.strip().split(" ")
                tmp_dict["chr_nr"] = line[0] # Fill the temporary dictionary with the right information per line
                tmp_dict["variant_ID"] = line[1]
                tmp_dict["position"] = line[2]
                tmp_dict["first_option"] = line[3]
                tmp_dict["second_option"] = line[4]
                line = line[5:] # From 5th position the genotypes
                for index in idx:
                    haps1.append(line[2 * index]) # First haplotype
                    haps2.append(line[(2 * index) + 1]) # Second haplotype
                combined_haps = []
                for (item1, item2) in zip(haps1, haps2): # Combine haps1 and haps2
                    combined_haps.append("".join(map("".join, zip(item1, item2))))
                combined_haps_string = "".join(list(combined_haps))
                tmp_dict["genotypes"] = combined_haps_string # Add genotypes to the temporary dict
                haps_file_info.append(tmp_dict)

    haps_file_out = open(haps_file_out, "w")
    cols = ["Column" + str(x) for x in range(len(haps_file_info[0]["genotypes"].replace(" ", "")))] # Make unique column name for all genotype columns
    w = csv.DictWriter(haps_file_out, ["chr_nr", "variant_ID", "position", "first_option", "second_option"] + cols, delimiter = " ")
    for snp in haps_file_info:
            trans = dict(zip(cols, snp["genotypes"].replace(" ", ""))) # Remove spaces between nucleotides
            snp_dict = snp.copy()
            snp_dict.update(trans)
            del snp_dict["genotypes"] # Remove original genotypes and use new ones
            w.writerow(snp_dict)
    haps_file_out.close()


def main():
	'''
	This is the main function
	'''
	chromosome, haps_file, sample_file, data_dir = argparser()
	list_idx = sample_chunker(sample_file, data_dir)
	for i, idx in enumerate(list_idx):
		haps_file_out = os.path.join(data_dir,"chr{}.offspring.{}.haps".format(chromosome,i))
		filter_haps(haps_file, idx, haps_file_out)

if __name__ == "__main__":
	main()
