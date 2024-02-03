#!/usr/bin/env python
# coding: utf-8



import matplotlib.pyplot as plt
import os
import sys

# Get the directory of the current script

# Function to read values from a file and return them as a list of floats
def read_values_from_file(file_path):
    # Open the file and read each line, converting each line to a float and stripping whitespace
    with open(file_path, 'r') as file:
        return [float(line.strip()) for line in file]
# Function to plot XY scores from a file
def plot_xy_scores(file_name):

    
    x_values = [i + 0.5 for i in range(20)]  # X-axis values from 0.5 to 19.5 with a step of 1
    

    # Extract the code (assumed format: "XY_scores.txt")    
    code = file_name.split('_')[0]
    y_values = read_values_from_file(file_name)# Read Y-axis values from the file
    plt.plot(x_values, y_values)  # Use code value as a label

    # Set labels and title
    plt.xlabel('distance')
    plt.ylabel('score')
    plt.title(f'{code} scores')

    # Save the plot with a file name containing the XY value
    
    plot_filename = f'{code}_Scores_Plot.png'
    plt.savefig(plot_filename)

    # Show the plot
    plt.show()

    #print(f"Plot saved as: {plot_filename}")



# Get a list of files in the specified directory ending with "_scores.txt"
files_list=[file for file in os.listdir('PATH_TO_FILES') if file.endswith("_scores.txt")]



# Iterate over each file and plot its XY scores
for file in files_list:
    plot_xy_scores(file)





