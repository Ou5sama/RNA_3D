#!/usr/bin/env python
# coding: utf-8

import math
import numpy as np
import pandas as pd
import os

base_pairs = ['AA', 'AU', 'AC', 'AG', 'UU', 'UC', 'UG', 'CC', 'CG', 'GG']
x_values = [i + 0.5 for i in range(20)]

def parse_pdb_file(file_path):
    coordinates = []
    with open(file_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('ATOM') and 'C3\'' in line:
                columns = line.split()
                atome = columns[3]

                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coordinates.append((atome,(x, y, z)))
    return coordinates

def compute_distance(atom1, atom2):
    return math.sqrt((atom1[0] - atom2[0])**2 + (atom1[1] - atom2[1])**2 + (atom1[2] - atom2[2])**2)

def compute_interatomic_distances(coordinates,base_pairs):
    distances = []
    for i in range(len(coordinates)):
        for j in range(i + 4, len(coordinates), 1):
            distance = compute_distance(coordinates[i][1], coordinates[j][1])
            if coordinates[i][0]+coordinates[j][0] in base_pairs:
                code=coordinates[i][0]+coordinates[j][0]
            else:
                code=coordinates[j][0]+coordinates[i][0]
            distances.append((code,distance))
    return distances

# Function to read values from a file and return them as a list of floats
def read_values_from_file(file_path):
    # Open the file and read each line, converting each line to a float and stripping whitespace
    with open(file_path, 'r') as file:
        return [float(line.strip()) for line in file]

# Function to create a dictionary of values for each base pair from a list of files
def create_values_dic(files_list, base_pairs):
     # Initialize an empty dictionary to store values for each base pair
    values_dic={base_pair:None for base_pair in base_pairs}
     # Iterate over each file in the list
    for file_name in files_list:
        # Extract the code (assumed format: "XY_scores.txt")
        code = file_name.split('_')[0]
        # Read Y-axis values from the file using the read_values_from_file function
        y_values = read_values_from_file(file_name)
         # Update the dictionary with the code as the key and corresponding Y-axis values
        values_dic[code]=y_values
    # Return the populated values dictionary
    return values_dic


# Function to interpolate scores based on distances using numpy's interp function
def interpolate(distances,values_dic,x_values):
    # Initialize dictionaries to store scores and distances for each key
    scores={}
    distances_dict = {}
    # Iterate over the distances and values
    for key, value in distances:
    # Check if the key is already in the dictionary
        if key in distances_dict:
            distances_dict[key].append(value)
        else:
        # If the key is not in the dictionary, create a new entry
            distances_dict[key] = [value]
    # Iterate over the keys and values in the distances dictionary        
    for key,value in distances_dict.items():
        try:
            # Use numpy's interp function to interpolate scores based on distances
            scores[key]=np.interp(value,x_values,values_dic[key])
        except KeyError:
            # Handle the case where a key is not found in the values dictionary
            print("this is a key")
    # Return the interpolated scores
    return scores
    



