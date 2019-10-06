#!/usr/bin/env python3
# -*- coding: utf-8 -*-


__author__ = "Shashwat Deepali Nagar"
__version__ = "1.0"
__email__ = "shashwat.d.nagar@gmail.com"


'''
Script for calculating population-specific allele frequencies from a VCF file.
'''


# Python imports
from sys import exit
from os.path import exists
from collections import defaultdict
from argparse import ArgumentParser
from gzip import open as gzip_open


def file_checks(input_file, grouping_file):
    '''
    Input:  input_file (str)
            Name of input file specified by the user
            
            grouping_file (str)
            Name of grouping file specified by the user

    Output: None
            Exits the program if files not found.
    '''
    if not exists(input_file):
        exit(f"Terminating.  {input_file} does not exist.  \nPlease check the specified input.")

    if not exists(grouping_file):
        exit(f"Terminating.  {grouping_file} does not exist.  \nPlease check the specified input.") 


def get_groupings(grouping_input):
    '''
    Input:  grouping_input (str)
            File that contains grouping information in a tab-separated format.  It
            can contain more than one grouping - i.e. there can be more than one 
            grouping system.
    Output: grouping_dict (defaultdict)
            A dictionary of dictionaries that contains grouping information for all
            individuals in a grouping file.
    '''
    
    # Storing different groupings based on headers in the grouping file
    grouping_dict = defaultdict(dict)
    
    # Keeping track of grouping systems (headers)
    grouping_headers = defaultdict(set)
    
    with open(grouping_input, 'r') as fh:
        
        # Getting header information
        header_names = {}
        header = fh.readline()
        header_fields = header.rstrip().split()
        for group_index in range(len(header_fields[1:])):
            header_names.update({group_index : header_fields[group_index + 1]})

        # Populating grouping list
        for line in fh:
            fields = line.rstrip().split()
            for group_index in range(len(fields[1:])):
                grouping_dict[fields[0]].update({header_names[group_index]: fields[group_index + 1]})
                grouping_headers[header_names[group_index]].add(fields[group_index + 1])
        
    return grouping_dict, grouping_headers
    


def get_frequencies(genotype_list, group_index_list):
    '''
    Input:  genotype_list (list)
            List of genotype calls from VCF file
            
            group_index_list (defaultdict)
            A dictionary of dictionaries that contains grouping information for all
            individuals in a grouping file.
    Output: frqs (dict)
            Dictionary of frequencies for different groups
    '''
    # Initializing a dictionary that will hold frequeny values for different groups
    sums = {}
    counts = {}
    frqs = {}
    
    # Initializing dictionaries with group names
    for group_name in group_index_list:
        sums[group_name] = 0
        counts[group_name] = 0
        frqs[group_name] = 0
    
    for individual_index in range(len(genotype_list)):
        for group_name in group_index_list:
            if individual_index in group_index_list[group_name]:
                sums[group_name] += genotype_list[individual_index].count('1')
                counts[group_name] += sum(character.isdigit() for character in genotype_list[individual_index])
        
    for group_name in group_index_list:
        frqs[group_name] = sums[group_name]/counts[group_name]
    
    return frqs, counts


def main():
    # Getting user input
    parser = ArgumentParser(description = "A population-specific SNP frequency calculator that let's you specify different population groups.")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument('--vcf', required=True, help = "Input VCF file - may be gzipped.")
    required.add_argument('--groups', required=True, help = "Grouping file - tab-separated file that defines group membership.  Each row is one individual followed by group(s).  See the Github () for a sample file.")
    required.add_argument('--output', required=True)

    optional.add_argument('--use-fid', action='store_true', help = "Use Family ID (the second to last word after splitting individual identifiers on '_'.  Defaults to using the individual ID (last word after split).  This is the ID that is used for group mapping - choose wisely.")
    parser._action_groups.append(optional)
    
    args = parser.parse_args()


    # Uncompressed input VCF file
    input_file = args.vcf

    # Can specify more than one grouping - REQUIRES a header
    grouping_file = args.groups

    # Ouput file location
    output_file = args.output

    # Use IID?
    if args.use_fid:
        use_iid = 2
    else:
        use_iid = 1

    # Check input files
    file_checks(input_file, grouping_file)

    grouping_info, grouping_systems = get_groupings(grouping_file)

    # Determine whether input file is Gzipped.  Only uses the name of the file.  
    # Can use magic bit, but that would require installing an additional dependency.
    if input_file.endswith("gz"):
        open_mode = gzip_open(input_file, 'rt')
    else:
        open_mode = open(input_file, 'r')

    with open_mode as fh:
        with open(output_file, 'w') as ofh:
            for line in fh:
                if not line.startswith('#'):
                    # Reading in data and calculating frequencies
                    fields = line.rstrip().split("\t")
                    
                    ofh.write("\t".join(fields[:5]))
                    freqs, count = get_frequencies(fields[9:], group_sets)
                    for group_name in freqs:
                        ofh.write(f'\t{freqs[group_name]}\t{count[group_name]}')
                    ofh.write('\n')
                    
                elif line.startswith("#CHROM"):
                    # Reading in header to get individual names
                    sample_names = [sample_name.split("_")[-use_iid] for sample_name in line.rstrip().split("\t")[9:]]
                    
                    # Creating group-specific sets of indices
                    group_sets = defaultdict(set)
                    
                    for individual_index in range(len(sample_names)):
                        for group_type in grouping_info[sample_names[individual_index]].values():
                            group_sets[group_type].add(individual_index)
                    
                    # Printing header to output file
                    ofh.write('Chr\tPosition\trsID\tRef\tAlt')
                    for group_name in group_sets:
                        ofh.write(f'\t{group_name}_freq\t{group_name}_count')
                    ofh.write('\n')



if __name__ == '__main__':
    main()