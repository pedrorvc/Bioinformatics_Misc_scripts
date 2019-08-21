#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 13:11:20 2019

@author: Pedro Cerqueira
"""

import argparse
import os

import pandas as pd
import numpy as np

from Bio.SeqUtils import GC
from Bio import SeqIO


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

def get_species_dict():
    """
    """
    
    # Read file
    species_scheme = pd.read_csv("scheme_species_map.tab", 
                              delimiter = '\t',
                              encoding = 'utf-8',
                              header = 0)
    
    # Replace Numpy NaN with blank string
    species_scheme = species_scheme.replace(np.nan, "", regex = True)
    
    # Concatenate the genus column with the species column
    species_scheme_concat = pd.DataFrame(species_scheme["GENUS"].str.cat(species_scheme["SPECIES"], sep = " "))
    
    # Rename concatenated column
    species_scheme_concat.rename(columns={"GENUS": "Species"}, inplace=True)
    
    # Add the scheme to search for in the report
    species_scheme_concat["Scheme"] = species_scheme["#SCHEME"]
    
    # Convert the dataframe to dictionary
    species_dict = species_scheme_concat.set_index('Species')["Scheme"].to_dict()
    
    return species_dict

def analyse_report(report, mlst, total_bps, nr_contigs, gc_content):
    
    """ Classify samples using report results
    
        Args:
            report (dict): Contains the assembly and annotation reports
            mlst (str): Species to compare with the annotation results
            total_bps (int): Total assembled base pairs
            nr_contigs (int): number of generated contigs
            gc_content (float): 
                
        Returns:
            
            result (dict) : Contains the results of the analysis
    """
    

    result = {}
    
    for record in report:
        
        result[record["Sample"]] = []
        
        if record["Total assembly length"] > total_bps:
            result[record["Sample"]].append(["FAIL", "Total_BPs"])
            
        if record["Number of contigs"] > nr_contigs:
            result[record["Sample"]].append(["FAIL", "Nr_contigs"])
            
        if record["GC content"] > gc_content:
            result[record["Sample"]].append(["FAIL", "GC_content"])
            
        if record["mlst"] != mlst:
            result[record["Sample"]].append(["FAIL", "Species"])
        
        if result[record["Sample"]] == []:
            result[record["Sample"]].append("PASS")
    
    return result


def pilon_report_analysis(species, assemblies, output, pilon_report_path, 
                          mlst_report_path, total_bps, nr_contigs, gc_content):
    """
    """
    
    # Read report files
    pilon_report = pd.read_csv(pilon_report_path,
                               encoding = 'utf-8')
    
    
    mlst_report = pd.read_csv(mlst_report_path, 
                              delimiter = '\t',
                              encoding = 'utf-8',
                              header = None,
                              usecols = [0,1])
    
    
    # Merge reports
    final_report = pilon_report.copy()
    
    final_report["mlst"] = mlst_report[1]
    
    #for idx, row in final_report.iterrows():
    #    print(idx, row["mlst"])
    
    # Convert merged report to dict
    final_report_dict = final_report.to_dict(orient = 'records')
    
    # Analyse report 
    result = analyse_report(final_report_dict, species, total_bps, nr_contigs, gc_content)
    
    
    # Build the final report
    final_report2 = final_report.copy()
    
    s = pd.Series(result, name='Result')
    
    s.index.name = 'Sample'
    
    s = s.reset_index()
    
    final_report2["Result"] = s["Result"]
    
    # Write file
    final_report2.to_csv(os.path.join(output, "final_report.tsv"), 
                         sep = '\t', 
                         encoding = 'utf-8', 
                         index = False)
    
    print("File written in {}".format(output))


def generate_report(assemblies):
    """
    """
    
    results = {}
    
    for path in assemblies:
        records = SeqIO.parse(path, "fasta")
        
        for rec in records:
            name = rec.name
            seq = rec.seq
            len_seq = len(seq)
            

def main(species, assemblies, output, pilon_report_path, mlst_report_path,
         total_bps, nr_contigs, gc_content):
    """
    """
    
    # get species mlst name 
    species_dict = get_species_dict()
    #print(species)
    try:
        species_mlst = species_dict[species]
    
    except KeyError:
        print("Invalid or unknown species. Please provide a valid species name.")
        return
    
    if assemblies:
        assemblies_list = check_if_list_or_folder(assemblies)
        
        generate_report(assemblies_list)
    
    if pilon_report_path and mlst_report_path:
        
        pilon_report_analysis(species_mlst, assemblies, output, pilon_report_path, mlst_report_path,
                              total_bps, nr_contigs, gc_content)
        
    print("Execution Finished")
    return
    



def parse_arguments():
    
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
        
    parser.add_argument('-s', '--species', type=str, required=True, dest='species', default=False,
                    help='Path to the mlst report file.')
    
    parser.add_argument('-i', '--assemblies', type=str, required=False, dest='assembly_path', default=False,
                        help='Path to the directory containing the assemblies.')
    
    parser.add_argument('-o', '--output', type=str, required=True, dest='output_path', default=False,
                        help='Path to the output directory.')
    
    parser.add_argument('--pilon_report', type=str, required=False, dest='pilon_report_path', default=False,
                        help='Path to the pilon report file.')
    
    parser.add_argument('--mlst_report', type=str, required=False, dest='mlst_report_path', default=False,
                    help='Path to the mlst report file.')
    
    parser.add_argument('--total_bps', type=str, required=False, dest='total_bps', default=2300000,
                help='Number of total assembled base pairs to test against.')
    
    parser.add_argument('--nr_contigs', type=str, required=False, dest='nr_contigs', default=350,
                    help='Number of contigs allowed for each assembly.')

    parser.add_argument('--gc_content', type=str, required=False, dest='gc_content', default=0.45,
                    help='Percentage of GC content allowed.')
    
    args = parser.parse_args()
    
    return [args.species, args.assembly_path, args.output_path, args.pilon_report_path,
            args.mlst_report_path, args.total_bps, args.nr_contigs, args.gc_content]
    
if __name__ == '__main__':
    
    args = parse_arguments()
    main(args[0], args[1], args[2], args[3],
         args[4], args[5], args[6], args[7])

    




    











