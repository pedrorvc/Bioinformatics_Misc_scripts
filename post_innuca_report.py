#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Pedro Cerqueira
    github: @pedrorvc
    
    
"""

import argparse
import os
from multiprocessing import Pool, Process, cpu_count
import shutil
import subprocess

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
    """ Creates a dictionary mapping species to mlst scheme names
        
        Example: 
            {"Streptococcus agalactiae" : "sagalactiae"}
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
            gc_content (float): Percentage of GC content allowed
                
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


def write_final_report(report, analysed_report, output):
    """ Writes a tab-delimited file of analysis results.
    
        Args:
            report (pandas.DataFrame): DataFrame containing the assembly and annotation results
            analysed_report (dict): Contains the results of the report analysis
            output (str): Path to the output directory
            
        Returns:
            Tab-delimited file in the specified output director
    """
    
    # Build the final report (DataFrame)
    final_report = report.copy()
    
    # Convert analysed report (dict) to pandas.Series
    s = pd.Series(analysed_report, name='Result')
    
    # Change the name of the index
    s.index.name = 'Sample'
    
    # Reset index
    s = s.reset_index()
    
    # Add the result Series to the final_report DataFrame
    final_report["Result"] = s["Result"]
    
    # Write file
    final_report.to_csv(os.path.join(output, "final_report.tsv"), 
                         sep = '\t', 
                         encoding = 'utf-8', 
                         index = False)
    
    print("File written to {}".format(output))


def pilon_report_analysis(species, assemblies, output, pilon_report_path, 
                          mlst_report_path, total_bps, nr_contigs, gc_content):
    
    """ Analyses the report generated by the Innuca pipeline
    
        Args:
            species (str): MLST scheme name of the species
            assemblies (list): Contains the assembly files to be analysed
            output (str): Path to the output directory
            pilon_report_path (str): Path to the pilon report file generated by Innuca
            mlst_report_path (str): Path to the mlst report file generated by Innuca
            total_bps (int): Total size in base pairs allowed for the analysis
            nr_contigs (int): Max number of contigs allowed for the analysis
            gc_content (float): Max percentage of GC content allowed for the analysis
            
        Returns:
            Report in tab-delimited format with the results of the pipeline and analysis.
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
    
    # Convert merged report to dict
    final_report_dict = final_report.to_dict(orient = 'records')
    
    # Analyse report 
    result = analyse_report(final_report_dict, species, total_bps, nr_contigs, gc_content)
    
    
    # Write the final report
    write_final_report(final_report, result, output)
      

def calc_n50(sizes):
    """ Calculates the N50 of an assembly.
    
        Args:
            sizes (list): Contains the sizes of the contigs from an assembly file
            
        Returns:
            n50 (int)
    """
    
    # Calculate the cumulative sum of the sizes
    csum = np.cumsum(sizes)
    
    # Get half of the sum of sizes
    n2 = int(sum(sizes) / 2)
    
    # Get index for csum >= N/2
    csumn2 = min(csum[csum >= n2])
    idx = np.where(csum == csumn2)
    
    # Get the N50
    n50 = sizes[int(idx[0])]
    
    return n50
    


def generate_report(assemblies, output, species, total_bps, nr_contigs, gc_content):
    """ Creates a report for assembly files
    
        Args:
            assemblies (list): Contains the assembly files to be analysed
            output (str): Path to the output directory
            species (str): MLST scheme name of the species
            total_bps (int): Total size in base pairs allowed for the analysis
            nr_contigs (int): Max number of contigs allowed for the analysis
            gc_content (float): Max percentage of GC content allowed for the analysis
            
        Returns:
            Report in tab-delimited format with assembly statistics 
            and subsequent analysis results
    """
    
    results = []
    
    for path in assemblies:
        
        # Get the sample name from the file
        sample = os.path.basename(path).split("_")[0]
        
        # Get the records of the assembly file
        records = list(SeqIO.parse(path, "fasta"))
        
        # Calculate the GC content
        all_gc_content = [GC(seq.seq) for seq in records]
        gc_content = np.mean(all_gc_content) / 100
        
        # Get the total number of contigs in the assembly file
        nr_contigs = len(records)
        
        # Get the contig sizes
        sizes = [len(seq) for seq in records]
        
        # Calculate the average contig size
        avg_size = np.average(sizes)
        
        # Calculate the total assembly length
        total_length = sum(sizes)
        
        # Calculate the N50
        n50 = calc_n50(sizes)
        
        # Determine missing data
        missing_data = sum([rec.seq.count("N") for rec in records])
         
        # Determine the allelic profile and save the identified species
        p = subprocess.run(["mlst", "--threads", str(cpu_count() - 2), path], capture_output=True)     
        mlst = p.stdout.decode("utf-8").split("\t")[1]

        # Save the results in a dictionary        
        results.append({"Sample" : sample, "Number of contigs" : nr_contigs,
                        "Average contig size" : avg_size, "N50" : n50,
                        "Total assembly length" : total_length, "GC content" : gc_content,
                        "Missing Data" : missing_data, "mlst" : mlst}) 
        
    # Analyse results
    analysed_report = analyse_report(results, species, total_bps, nr_contigs, gc_content)
    
    # Convert dictionary into pandas DataFrame
    report = pd.DataFrame(results)
    
    # Write the final report
    write_final_report(report, analysed_report, output)


    

def main(species, assemblies, output, pilon_report_path, mlst_report_path,
         total_bps, nr_contigs, gc_content):

    
    # get species mlst scheme name 
    species_dict = get_species_dict()

    try:
        species_mlst = species_dict[species]
    
    except KeyError:
        print("Invalid or unknown species. Please provide a valid species name.")
        return
    
    if assemblies:
        
        # Check if mlst is in PATH
        if shutil.which("mlst") == None:
            print("mlst is not on your PATH. Stopping execution...")
            return
        
        assemblies_file = check_if_list_or_folder(assemblies)
        
        listGenes = []
        with open(assemblies_file, "r") as gf:
            for gene in gf:
                gene = gene.rstrip('\n')
                listGenes.append(gene)
        listGenes.sort()

        try:
            os.remove("listGenes.txt")
        except:
            pass
        
#        proc = ""
#        
#        p = Pool(processes = proc)
#        r = p.map_async(generate_report, listGenes)
        
        generate_report(listGenes, output, species_mlst, total_bps, nr_contigs, gc_content)
    
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








