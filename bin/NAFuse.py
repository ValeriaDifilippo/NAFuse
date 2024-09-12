#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NAFuse v2.0.0
Created on Tue Sep 10
@author: Valeria Difilippo
"""

import argparse
import glob
import os

# Set up the argument parser
parser = argparse.ArgumentParser(description="NAFuse allows for the detection of breakpoints in both partner genes/regions on both the RNA and DNA levels without having to manually search output files generated during downstream processing of BAM-files.")
parser.add_argument('-dir', required=True, help="Directory containing the DNA and RNA files")
parser.add_argument('-o', required=True, help="Output directory for gene fusions results")

# Add --version flag to display the script's version
parser.add_argument('--version', action='version', version='NAFuse v2.0.0')

# Parse the arguments
args = parser.parse_args()

# Define the input directory and output directory based on arguments
input_dir = args.dir
output_dir = args.o

# Ensure the output directory exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Find all DNA and RNA file pairs
dna_files = glob.glob(os.path.join(input_dir, '*_DNA.bed'))
rna_files = glob.glob(os.path.join(input_dir, '*_RNA.txt'))

# Create a dictionary to map file prefixes to RNA files
rna_file_map = {os.path.basename(f).replace('_RNA.txt', ''): f for f in rna_files}

# Process each DNA file and find its matching RNA file
for dna_file in dna_files:
    prefix = os.path.basename(dna_file).replace('_DNA.bed', '')

    # Check if there is a matching RNA file
    if prefix in rna_file_map:
        rna_file = rna_file_map[prefix]
        output_file = os.path.join(output_dir, f'{prefix}_supporting_reads.txt')
        gene_pairs_file = os.path.join(output_dir, f'{prefix}_validated_gene_fusion.txt')  # New file for unique gene pairs

        print(f'Processing DNA: {dna_file} with RNA: {rna_file}, output will be {output_file} and {gene_pairs_file}')

        # Read RNA file and store the gene pairs in a set for faster lookup
        gene_pairs = set()
        with open(rna_file, 'r') as rna_f:
            for line in rna_f:
                genes = line.strip().split()
                if len(genes) == 2:  # Ensuring there are exactly two genes
                    gene_pairs.add(tuple(genes))

        # Set to store unique gene pairs for the second output file
        unique_gene_pairs = set()

        # Open the output file for writing
        with open(output_file, 'w') as output, open(gene_pairs_file, 'w') as gene_output:
            # Read DNA file and check if the gene pairs exist in the RNA file
            with open(dna_file, 'r') as dna_f:
                line_count = 0  # Count the lines processed

                for line in dna_f:
                    line_count += 1
                    columns = line.strip().split()

                    # Ensure the line has at least 9 columns (since we access column[4] and column[8])
                    if len(columns) >= 9:
                        gene1, gene2 = columns[4], columns[8]  # Extract gene pair from DNA file

                        # Check if the pair or the reverse of the pair exists in RNA file
                        if (gene1, gene2) in gene_pairs or (gene2, gene1) in gene_pairs:
                            output.write(line)  # Write the matching line to the output file
                            unique_gene_pairs.add((gene1, gene2))  # Add the unique pair to the set
                            print(f"Match found: {gene1}, {gene2}")  # Optional: print the match
                    else:
                        # Debugging: Not enough columns in this line
                        print(f"Skipping line {line_count}: Insufficient columns ({len(columns)} found), requires at least 9.")

                # Write the unique gene pairs to the second output file
                for gene1, gene2 in unique_gene_pairs:
                    gene_output.write(f"{gene1}\t{gene2}\n")

                # Debugging: Print a message when analysis is complete
                print(f"Finished processing {line_count} lines from {dna_file}.")
                print(f"Outputs are written to {gene_pairs_file} and {output_file}.")

    else:
        print(f"No matching RNA file for {dna_file}")
