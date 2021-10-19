# E_coli_sORF_discovery
Readme
10/19/21

There are 4 scripts, each of which comes with its own set of input files. Details of each script are given below. Where possible, input files are also provided. Note that some input files are too large to include.


***Script #1***
***Ribo-RET_Ribo-Api-Pmn_ORF_caller_10_18_21.py***

This Python script takes data from Ribo-RET and Ribo-ApiPmn experiments and infers the position of ORFs. It also creates output files, including a graph for a publication-ready figure. Required input files are listed below, with example files included (these files are available on GitHub). 

The script runs in Python v3 and requires that the following modules are installed:
matplotlib.pyplot 
numpy as np
pandas
seaborn
datetime
pathlib
subprocess
os
math
random
scipy

The script also requires a local installation of the ViennaRNA Package that includes RNAfold.exe

The user is prompted to provide file paths for input and output files, and filenames for input files. Output files are written to a folder with the date and time included in the folder name.


------------------------------------------------------------------------------------------------------

List of required input files, with names of example files

1. .gff listing sequence read coverage genome-wide on both strands for Ribo-RET data. Must be at least 5 columns, with column 3 [counting starts at column 0] listing the genome coordinate and column 5 listing the RPM coverage. Negative values in column 5 indicate positions on the minus strand.
Example file: Ribo-seq.gff

2. .gff listing sequence read coverage genome-wide on both strands for Ribo-ApiPmn data. Must be at least 5 columns, with column 3 [counting starts at column 0] listing the genome coordinate and column 5 listing the RPM coverage. Negative values in column 5 indicate positions on the minus strand.
Example file: Ribo-seq.gff

3. .fna with one line of header and a single line of genome sequence (no gaps or line breaks)
Example file: CP001509_3.fna

4. .gff listing all ORFs from the existing genome annotation. Must be at least 6 columns, with column 3 [counting starts at column 0] listing the left-end coordinate, column 4 listing the right-end coordinate, and column 6 listing the strand.
Example file: BL21_genes.gff

------------------------------------------------------------------------------------------------------
List of output files (filenames contain date or time, shown here as #date# and #time#, respectively):

******************************************************************************************************

Filenames

ret_codon_frequency_#date#_#time#.txt

api_codon_frequency_#date#_#time#.txt

Description

Tab-delimited .txt showing the frequency of all trinucleotide sequences at positions from -50 to +10 relative to start codons for ret data and stop codons for Api-Pmn data
******************************************************************************************************

Filenames

RET all peaks_#date#_#time#.txt

Api all peaks_#date#_#time#.txt

Description

Tab-delimited .txt showing peak coordinate, strand and peak height

******************************************************************************************************

Filename

Ret-full report_#date#_#time#.txt

Description

Tab-delimited .txt showing Ret peak height, start codon, gap, strand, start codon coordinate, stop codon coordinate, nt sequence, aa sequence

******************************************************************************************************

Filename

Api-full report_#date#_#time#.txt

Description

Tab-delimited .txt showing ApiPmn peak coordinate, peak height, stop codon, gap, strand, stop codon coordinate

******************************************************************************************************

Filenames

Annotated_RBS_sequences_Annotated_Isoform_Novel_MFE_start_codon_scores_#date#.txt

Isoform_RBS_sequences_Annotated_Isoform_Novel_MFE_start_codon_scores_#date#.txt

Novel_RBS_sequences_Annotated_Isoform_Novel_MFE_start_codon_scores_#date#.txt

random_sequences_Annotated_Isoform_Novel_MFE_start_codon_scores_#date#.txt

Description

.fasta format files listing sequences used to predicted minimum free energy RNA structure in the regions around annotated/novel/isoform ORF start codons (-40 to +20 nt) and random sequences

******************************************************************************************************

Filenames

Annotated_MFE_deltaG_scores_Annotated_Isoform_Novel_MFE_start_codon_scores_#date#.txt

Isoform_MFE_deltaG_scores_Annotated_Isoform_Novel_MFE_start_codon_scores_#date#.txt

Novel_MFE_deltaG_scores_Annotated_Isoform_Novel_MFE_start_codon_scores_#date#.txt

Random_MFE_deltaG_scores_Annotated_Isoform_Novel_MFE_start_codon_scores.txt

Description

Single column text files listing delta G scores for predicted minimum free energy RNA structure in the regions around annotated/novel/isoform ORF start codons (-40 to +20 nt) and random sequences

******************************************************************************************************

Filename

Annotated_Isoform_Novel_MFE_start_codon_scores#date#.png

Description

Strip plot of delta G scores for predicted minimum free energy RNA structures around start codons of annotated/novel/isoform ORF start codons (-40 to +20 nt) and random sequences

******************************************************************************************************

Filename

MFE_statistical_summary_starts_#date#_#time#.txt

Description

Statistical comparisons of delta G scores for minimum free energy RNA structure predictions for regions around start codons of annotated/novel/isoform ORFs

******************************************************************************************************

Filename

ORF_list_#date#_#time#.txt

Description

Tab-delimited lists of ORFs showing: Ribo-RET peak coordinate, Ribo-Ret peak height, start codon, start codon to Ribo-RET peak distance (gap), strand, start codon coordinate, stop codon, stop codon coordinate, nt sequence, aa sequence, Ribo-ApiPmn peak coordinate, Ribo-ApiPmn peak height, stop codon to ApiPmn peak distance (gap), ORF category (Annotated, Isoform, Novel)

******************************************************************************************************

Filename

Summary_#date#.txt

Summary of the program output listing: Threshold, Window, Span, Fold over Span mean threshold, Acceptable start codon sequences and gaps, Acceptable stop codon sequences and gaps 

For each cycle, tab-delimited .txt showing number of Ribo-RET peaks, number of Ribo-RET peaks associated with a start codon, number of Ribo-ApiPmn peaks, number of Ribo-ApiPmn peaks associated with a stop codon after filtering for overlap with Ribo-RET peaks, number of unique Annotated ORFs, number of unique Isoform ORFs, number of unique Novel ORFs


------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------

***Script #2***
***adaptor_trimming.py***

This Python script takes a .fastq file from a Ribo-seq experiment where the RNA was ligated during library preparation. It trims the reads at the first instance of 'CTGTAGGCACC' between positions 19 and 44 of the read. If there is no 'CTGTAGGCACC' in that range, the read is discarded.

The script runs in Python v3

The user is prompted to provide file paths for input and output .fastq files.


------------------------------------------------------------------------------------------------------



Required input file

1. .fastq file from a Ribo-seq experiment where the RNA was ligated during library preparation



------------------------------------------------------------------------------------------------------


Output file:


******************************************************************************************************
Filename
Provided by user

Description
.fastq file with trimmed reads, but the same IDs as the original .fastq


------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------

***Script #3***
***write_gff_from_two_SAMs.py***

This Python script takes two .sam files, one where reads were mapped to the plus strand and one where reads were mapped to a reverse-complemented (i.e. minus strand) genome sequence. The script generates a .gff with coverage values at every genome position on both strands, where coverage is only calculated for sequence read 3' ends, i.e. a single genome position per read.

The script runs in Python v3

The user is prompted to provide file paths for input (x2) .sam files and output .gff file.

------------------------------------------------------------------------------------------------------



Required input files

1. .sam file from a Ribo-seq experiment .fastq that was mapped to the plus strand of the genome
2. .sam file from a Ribo-seq experiment .fastq that was mapped to the minus strand of the genome


------------------------------------------------------------------------------------------------------


Output file:


******************************************************************************************************
Filename
Provided by user

Description
.gff file with sequence read coverage at each genome position with non-zero coverage
Column 0: 'NA'
Column 1: ‘Agilent’
Column 2: ‘gff_label’
Column 3: Genome coordinate
Column 4: Genome coordinate
Column 5: Coverage. Negative values indicate coverage on the minus strand of the genome
Column 6: '.'
Column 7: '.'
Column 8: '.'

------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------

***Script #4***
***mask and normalize gff 9_13_21.py***

This Python script takes a .gff file, removes coverage at positions specified in a "mask" file, and normalizes the remaining values as RPM.

The script runs in Python v3

The user is prompted to provide file paths for input and output .gff files, and the mask .txt file


------------------------------------------------------------------------------------------------------



Required input files

1. .gff file, expected format:
Column 0: Can be anything
Column 1: Can be anything
Column 2: Can be anything
Column 3: Genome coordinate
Column 4: Genome coordinate
Column 5: Coverage. Negative values indicate coverage on the minus strand of the genome
Column 6: '.'
Column 7: '.'
Column 8: '.'

2. Tab-deliminted .txt file with strand/left coordinate/right coordinate for a set of regions to be masked in the .gff file
Example file: regions_to_mask.txt


------------------------------------------------------------------------------------------------------


Output file:


******************************************************************************************************
Filename
Provided by user

Description
.gff file with sequence read coverage at each genome position with non-zero coverage
Column 0: Same as input .gff file
Column 1: Same as input .gff file
Column 2: 'NA3'
Column 3: Genome coordinate
Column 4: Genome coordinate
Column 5: Coverage (RPM). Negative values indicate coverage on the minus strand of the genome
Column 6: Same as input .gff file
Column 7: Same as input .gff file
Column 8: Same as input .gff file


