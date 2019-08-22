# Bioinformatics_Misc_scripts
Collection of Bioinformatics Python scripts for various analyses and tasks

## post_innuca_report.py

This script serves to analyse assemblies generated from a pipeline like Innuca and create a report based on that analysis. The sequence type(ST) of the assemblies is also determined during the analysis.
    
It also accepts a pilon and mlst reports generated by the Innuca pipeline and analyses it in order to filter false positives of the quality control checks.

### Usage

##### Report generation
To generate a report simply give it the path to the directory where teh assembly files are located. The files must be in `.fasta` format. The species needs to be provided beacuse the script will also perform multilocus sequence typing using [mlst](https://github.com/tseemann/mlst).

    % post_innuca_report.py -s "Streptococcus agalactiae" -i path/to/assemblies/dir \ 
                            -o path/to/output/dir --cpu 6 --total_bps 2300000 \
                            --nr_contigs 350 --gc_content 0.45
                            

##### Report analysis
The analyse a report, give it the path to the pilon and mlst report files.  

    % post_innuca_report.py -s "Streptococcus agalactiae" --pilon_report path/to/pilon/report \
                            --mlst_report path/to/mlst/report -o path/to/output/dir \
                            --cpu 6 --total_bps 2300000 \
                            --nr_contigs 350 --gc_content 0.45

### Dependencies
* [mlst](https://github.com/tseemann/mlst)
* [pandas](https://pandas.pydata.org/)
* [numpy](https://www.numpy.org/)
* [Biopython](https://biopython.org/)
