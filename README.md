# M2_GENIOMHE_RNA
The 1_RNA_3D script takes as an input the pdb file and calculate the distacences then the frequences (obs and ref) to out put in the end 10 files for each of the base pairs : 'AA', 'AU', 'AC', 'AG', 'UU', 'UC', 'UG', 'CC', 'CG', 'GG'.
The seconde script, 2_profile_plot, takes the files of produced by the 1_RNA_3D to plot the scores.
The third script, 3_scoring, interpolate scores based on the computed interatomic distances using numpy's interp function from a pdb file.