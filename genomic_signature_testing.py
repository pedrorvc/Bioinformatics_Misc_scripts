#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 17:20:28 2019

@author: pcerqueira
"""

import os
import random

import numpy as np
import pandas as pd
import seaborn as sns

from Bio import SeqIO

from itertools import islice
from collections import Counter

from sklearn.manifold import TSNE
from matplotlib import pyplot as plt


def k_mers(sequence, k):
    """ Divide string into substrings of length k
    
        Args:
            sequence (str): String to be divided
            k (int): Length of the substring to be created
        Returns:
            generator containing the resulting kmers
        
    """
    it = iter(sequence)
    result = tuple(islice(it, k))
    if len(result) == k:
        yield "".join(result)
    for elem in it:
        result = result[1:] + (elem,)
        yield "".join(result)
       
        
def normalized(a, axis=-1, order=2):
    """ Normalize arrays
    
        Args:
            a (np.array): Numpy arrar to be normalized
        Returns:
            normalized numpy array
    """
    l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
    l2[l2==0] = 1
    return a / np.expand_dims(l2, axis)


def is_fasta(filename: str):
    """ Checks if a file is a FASTA file.

        Args:
            filename (str): the full path to the FASTA file

        Returns:
            True if FASTA file,
            False otherwise
    
    """
    
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        
        # returns True if FASTA file, False otherwise
        return any(fasta)
    

def check_if_list_or_folder(folder_or_list):
    """ Checks if the input is a file or a directory.

        Args: 
            folder_or_list (str): the full path to the file or directory

        Returns:
            list_files (str) if folder_or_list is a path to a file,
            list_files (list) if folder_or_list is a path to a directory,
            Raises Exception otherwise
    """
    
    # check if input argument is a file or a directory
    if os.path.isfile(folder_or_list):
        list_files = folder_or_list
    
    elif os.path.isdir(folder_or_list):
        
        fasta_files = []
            
        for genome in os.listdir(folder_or_list):
                
            genepath = os.path.join(folder_or_list, genome)
            
            # do not continue if genepath is a dir
            if os.path.isdir(genepath):
                continue
            
            # check if file is a FASTA file
            if is_fasta(genepath):
                fasta_files.append(os.path.abspath(genepath))
        
        # if there are FASTA files
        if fasta_files:
            # store full paths to FASTA files
            with open("listGenes.txt", "w") as f:
                for genome in fasta_files:
                    f.write(genome + "\n")
        else:
            raise Exception("There were no FASTA files in the given directory. Please provide a directory \
                            with FASTA files or a file with the list of full paths to the FASTA files.")

        list_files = "listGenes.txt"
    
    else:
        raise Exception("Input argument is not a valid directory or file with a list of paths. \
                        Please provide a valid input, either a folder with FASTA files or a file with \
                        the list of full paths to FASTA files (one per line).")

    return list_files


# Get data dirs
tutorial_gbs = "/home/pcerqueira/DATA/chewbbaca/chewBBACA_tutorial/complete_genomes"
legio_test = "/home/pcerqueira/github_repos/nothing_else_meta_testing/legio_complete"

# Get input files
fiels = check_if_list_or_folder(tutorial_gbs)

listGenes = []
with open(fiels, "r") as gf:
    for gene in gf:
        gene = gene.rstrip('\n')
        listGenes.append(gene)
listGenes.sort()

try:
    os.remove("listGenes.txt")
except:
    pass

# Obtain list containing the kmer frequencies (highly inefficient)
# Subests are created to ensure that the resulting arrays have the same length
# should change kmer function
hey = []
for i in listGenes:
    for rec in SeqIO.parse(i, "fasta"):
        # get seqeuence
        seq = str(rec.seq)
        # Divide sequence in kmers
        kmer = k_mers(seq, 6)
        # Calculate kmer frequency
        kmer_freq = Counter(kmer)
        # Obtain list of kmer freqeuncies
        val_list = list(kmer_freq.values())
        # Create a subset of the kmer frequencies
        # Ensures that the arrays have the same length
        subset_val_list = random.sample(val_list, 50)
        hey.append(subset_val_list)

# Labels of the genoems
ya = ["GBS"]*len(hey)


# rinse and repeat above for second species
fiels_legio = check_if_list_or_folder(legio_test)


listGenes2 = []
with open(fiels_legio, "r") as gf:
    for gene in gf:
        gene = gene.rstrip('\n')
        listGenes2.append(gene)
listGenes2.sort()

try:
    os.remove("listGenes.txt")
except:
    pass

hey_legio = []

for i in listGenes2:
    for rec in SeqIO.parse(i, "fasta"):
        seq = str(rec.seq)
        kmer = k_mers(seq, 6)
        kmer_freq = Counter(kmer)
        val_list = list(kmer_freq.values())
        subset_val_list = random.sample(val_list, 50)
        hey_legio.append(subset_val_list)
        
ya2 = ["Legio"]*len(hey_legio)


# Create arrays 
length_gbs = max(map(len, hey))
arrr_gbs = np.array([xi+[None]*(length_gbs-len(xi)) for xi in hey])
norm_arrr_gbs = normalized(arrr_gbs)

length_legio = max(map(len, hey_legio))
arrr_legio = np.array([xi+[None]*(length_legio-len(xi)) for xi in hey_legio])
norm_arrr_legio = normalized(arrr_legio)


#stack arrays and labels
hey_hey = np.row_stack((arrr_gbs, arrr_legio)) 

joinedList = [*ya, *ya2]

# Perform dimensionality reduction using tSNE
tsne = TSNE(n_components=2, random_state=0)

X_2d = tsne.fit_transform(hey_hey)

tsne_df = pd.DataFrame({'X': X_2d[:,0],
                       'Y': X_2d[:,1],
                       'species': joinedList})

tsne_df.head()

# Visualize scatterplot
plt.figure(figsize=(8,8))
sns.scatterplot(x="X", y="Y",
              data=tsne_df,
              hue="species",
              palette=["red", "green"],
              legend="full")
plt.xlim(-10, 10)
plt.ylim(-10, 10)

