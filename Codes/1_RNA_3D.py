#!/usr/bin/env python
# coding: utf-8



import math
import numpy as np
import pandas as pd





# Function to parse a PDB file and extract relevant atom coordinates
def parse_pdb_file(file_path):
    # Initialize a list to store extracted atom coordinates
    coordinates = []
    # Open the PDB file in read mode
    with open(file_path, 'r') as pdb_file:
        # Iterate over each line in the file
        for line in pdb_file:
            # Check if the line starts with 'ATOM' and contains information about 'C3\'' atom
            if line.startswith('ATOM') and 'C3\'' in line: #get lanes that contains information about C3' atom
                # Split the line into columns
                columns = line.split()
                # Get the nucleotide information from the columns
                atome = columns[3]
                #Get coordinates of each C3' atom of the nucleotide
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                # Append the nucleotide and its coordinates to the list
                coordinates.append((atome,(x, y, z)))
    # Return the list of extracted atom coordinates
    return coordinates

# Function to compute Euclidean distance between two atoms
def compute_distance(atom1, atom2):
    # Calculate the Euclidean distance using the formula: sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
    return math.sqrt((atom1[0] - atom2[0])**2 + (atom1[1] - atom2[1])**2 + (atom1[2] - atom2[2])**2)

# Function to compute interatomic distances for a list of atom coordinates
def compute_interatomic_distances(coordinates):
    # Initialize a list to store computed distances
    distances = []
    # Iterate over each pair of coordinates in the input list
    for i in range(len(coordinates)):
        # Only consider residues separated by at least 3 positions on the sequence
        for j in range(i + 4, len(coordinates), 1):
             # Calculate the Euclidean distance between the atoms at positions i and j
            distance = compute_distance(coordinates[i][1], coordinates[j][1])
             # Append the combination of nucleotides and the computed distance to the distances list
            distances.append((coordinates[i][0]+coordinates[j][0],distance))
    # Return the list of computed distances
    return distances

#def compute_frequencies(distances, bin_width=1, max_distance=20):
#    bins = np.arange(0, max_distance + bin_width, bin_width)
#    frequencies, _ = np.histogram(distances, bins=bins)
#    return frequencies / len(distances)

# Function to process distances by rounding them to the nearest integer
def process_distances(distances):
    # Initialize a list to store processed distances
    new_distances = []
    # Iterate over each pair of distances in the input list
    for pair in distances:
        # Unpack the pair into code and original distance value
        code, value = pair
        # Round the original distance value to the nearest integer using ceil (ceiling) function
        new_value = np.ceil(value)
         # Append the code and the rounded distance value to the new_distances list
        new_distances.append((code, new_value))
    # Return the list of processed distances
    return new_distances

# Function to calculate observed frequencies for a specific code (base pair)
def calculate_obs_frequencies(code, processed_distances):
    # Convert the code to a set for flexible matching
    code_set = set(code)

    # Filter tuples based on the condition
    filtered_tuples = [(code, pair[1]) for pair in processed_distances if set(pair[0]) == code_set]

    # Calculate frequencies of values between 1 and 20
    frequencies = {i: 0 for i in range(1, 21)}

    for _, value in filtered_tuples:
        if 1<=value<=20:
            frequencies[int(value)] += 1
    
    
    for frequency in frequencies.keys():
        frequencies[frequency]=frequencies[frequency]/len(filtered_tuples)
    
    return (code,frequencies)


# Function to compute reference frequencies for all processed distances
def compute_reference_frequencies(processed_distances):
    # Initialize a dictionary to store frequencies for values between 1 and 20
    frequencies = {i: 0 for i in range(1, 21)}
     # Iterate over processed distances
    for _, value in processed_distances:
        # Check if the value is within the range of 1 to 20
        if 1<=value<=20:
            # Increment the corresponding frequency in the dictionary
            frequencies[int(value)] += 1
    # Compute frequencies by dividing by the total number of processed distances
    for frequency in frequencies.keys():
        frequencies[frequency]=frequencies[frequency]/len(processed_distances)
    # Return the computed frequencies
    return frequencies

# Function to compute log ratios between observed and reference frequencies
def compute_log_ratio(observed_frequencies, reference_frequencies):
    # Add a small constant to avoid division by zero
    epsilon = 1e-10
    observed_frequencies = observed_frequencies + epsilon
    reference_frequencies = reference_frequencies + epsilon

    # Calculate log ratios, handle situations where observed_frequencies are zero
    log_ratios = -np.log(observed_frequencies / reference_frequencies)
    log_ratios[np.isinf(log_ratios)] = 10  # Set infinity to an arbitrary maximum value (e.g., 10)
    # Return the computed log ratios
    return log_ratios

# Function to write computed scores to a file
def write_scores_to_file(scores, base_pair):
    # Generate a file name based on the current base pair
    file_name = f'{base_pair}_scores.txt'
    # Open the file in write mode
    with open(file_name, 'w') as output_file:
        # Iterate over each score in the list
        for score in scores:
            # Write the score to the file with 6 decimal places
            output_file.write(f'{score:.6f}\n')

# Main function to execute the analysis
def main():
    pdb_files = ['1c2x.pdb']  # replace with your PDB files
    base_pairs = ['AA', 'AU', 'AC', 'AG', 'UU', 'UC', 'UG', 'CC', 'CG', 'GG']
# Iterate over each PDB file in the list
    for pdb_file in pdb_files:
        coordinates = parse_pdb_file(pdb_file)# Parse the PDB file to extract atom coordinates
        distances=compute_interatomic_distances(coordinates)# Compute interatomic distances based on the extracted coordinates
        processed_distances=process_distances(distances)# Process the distances (rounding to the nearest integer)
        
        reference_frequencies=compute_reference_frequencies(processed_distances)# Compute reference frequencies for all processed distances
        reference_frequencies_array = np.array([reference_frequencies[i] for i in range(1, 21)])# Convert reference frequencies to a numpy array for efficient computation
         # Iterate over each base pair
        for base_pair in base_pairs:
            # Calculate observed frequencies for the current base pair
            _,observed_frequencies = calculate_obs_frequencies(base_pair, processed_distances)
            # Convert observed frequencies to a numpy array for efficient computation
            observed_frequencies_array = np.array([observed_frequencies[i] for i in range(1, 21)])
            # Compute log ratios between observed and reference frequencies
            scores = compute_log_ratio(observed_frequencies_array, reference_frequencies_array)
 
            
            # Cap log ratios at an arbitrary maximum value (e.g., 10)
            scores = np.minimum(scores, 10)
            # Write the computed scores to a file for the current base pair
            write_scores_to_file(scores, base_pair)
        
if __name__ == "__main__":
    main()






