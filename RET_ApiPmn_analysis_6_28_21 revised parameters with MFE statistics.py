
#6/28/21 - include all IERF/TERF positions in the output files, even if they don't match a start/stop codon
#switch colours on the strip plots

    
#############
# Functions #
#############




###########
# revcomp #
###########

#takes a nucleotide sequence as input and returns the reverse complement
def revcomp(sequence):
    revseq = sequence[::-1]

    revseq = revseq.replace('A','k')
    revseq = revseq.replace('C','l')
    revseq = revseq.replace('G','m')
    revseq = revseq.replace('T','n')

    revseq = revseq.replace('a','k')
    revseq = revseq.replace('c','l')
    revseq = revseq.replace('g','m')
    revseq = revseq.replace('t','n')
    
    revseq = revseq.replace('k','T')
    revseq = revseq.replace('l','G')
    revseq = revseq.replace('m','C')
    revseq = revseq.replace('n','A')

    return revseq




#################
# ORFtranslator #
#################

#takes a nucleotide sequence as input and returns a list of three protein sequences, one for each frame
def ORFtranslator(sequence):

    #lists of amino acids and codons
    aa_codes =  ['*'  ,'*'  ,'*'  ,'I'  ,'I'  ,'I'  ,'L'  ,'L'  ,'L'  ,'L'  ,'L'  ,'L'  ,'V'  ,'V'  ,'V'  ,'V'  ,'F'  ,'F'  ,'M'  ,'C'  ,'C'  ,'A'  ,'A'  ,'A'  ,'A'  ,'G'  ,'G'  ,'G'  ,'G'  ,'P'  ,'P'  ,'P'  ,'P'  ,'T'  ,'T'  ,'T'  ,'T'  ,'S'  ,'S'  ,'S'  ,'S'  ,'S'  ,'S'  ,'Y'  ,'Y'  ,'W'  ,'Q'  ,'Q'  ,'N'  ,'N'  ,'H'  ,'H'  ,'E'  ,'E'  ,'D'  ,'D'  ,'K'  ,'K'  ,'R'  ,'R'  ,'R'  ,'R'  ,'R'  ,'R'  ]
    aa_codons = ['TAA','TAG','TGA','ATT','ATC','ATA','CTT','CTC','CTA','CTG','TTA','TTG','GTT','GTC','GTA','GTG','TTT','TTC','ATG','TGT','TGC','GCT','GCC','GCA','GCG','GGT','GGC','GGA','GGG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','TCT','TCC','TCA','TCG','AGT','AGC','TAT','TAC','TGG','CAA','CAG','AAT','AAC','CAT','CAC','GAA','GAG','GAT','GAC','AAA','AAG','CGT','CGC','CGA','CGG','AGA','AGG']

    #codon --> amino acid dictionary
    amino_acids_dict = {}

    #populate the dictionary
    for i in range(0,len(aa_codes)):  
        amino_acids_dict[aa_codons[i]] = aa_codes[i]

    #this will be a list of amino acids sequences from the three possible reading frames
    open_reading_frames = []


    #go through each of the 3 possible frames
    for frame in range(0,3):

        #sequence of the protein for the corresponding nucleotide sequence (input for the function)        
        protein_sequence = ''

        #go through the nucleotide sequence 3 nt at a time
        for position in range(frame, len(sequence), 3):

            #break if there isn't a full codon at the end of the sequence
            if len(sequence[position : position + 3]) < 3:
                break

            #translate the current codon
            protein_sequence = protein_sequence + amino_acids_dict[sequence[position : position + 3]]
        open_reading_frames.append(protein_sequence)

    return(open_reading_frames)





######################
# coord_cov_selected #
######################

#this function takes Ribo-RET or Ribo-ApiPmn sequence read coverage data for a genome
#it is written to use a gff
#it then identifies positions that have coverage above threshold (sequence coverage threshold defined within the program)
#these positions must also be the highest or joint highest in a region of sequence defined in length by "window"
#and they must have "minfold" times as much coverage as the average of all positions within a region defined in length by "span"
def coord_cov_selected(len_of_genome, threshold, window, span, min_fold, gff_filename):#GFF filename
              
    #dictionary linking genome position to normalized sequence read coverage, plus strand only
    plus_strand_coverage_by_position = {}

    #dictionary linking genome position to normalized sequence read coverage, minus strand only
    minus_strand_coverage_by_position = {}

    #plus strand positions that meet all the criteria for being called as a "peak"
    plus_coord_selected = []

    #minus strand positions that meet all the criteria for being called as a "peak"
    minus_coord_selected = []

    #opens the gff to populate the dictionaries linking genome position to normalized sequence read coverage
    f=open(gff_filename,'r')

    
    for x in f:

        #position [5] in the .split line from the gff is the normalized sequence read coverage
        #the function assumes that negative coverage values indicate minus strand reads
        #requiring that x.split()[5] be > 0 selects for plus strand data
        if float(x.split()[5]) > 0:

            #position [4] in the .split line from the gff indicates genome position
            #this line creates a dictionary pair with position and coverage
            plus_strand_coverage_by_position[int(x.split()[4])] = float(x.split()[5])

            #wrap the first 1000 nt to the end of the genome since the genome is circular
            if int(x.split()[4]) < 1000:
                plus_strand_coverage_by_position[int(x.split()[4]) + len_of_genome] = float(x.split()[5])

        #repeat for the minus strand
        elif float(x.split()[5]) < 0:

            minus_strand_coverage_by_position[1 + len_of_genome - int(x.split()[4])] = -float(x.split()[5])
            if 1 + len_of_genome - int(x.split()[4]) < 1000:
                minus_strand_coverage_by_position[1 + len_of_genome - int(x.split()[4]) + len_of_genome] = -float(x.split()[5])

    f.close()

    #looks at every dictionary entry (position/coverage) for the plus strand
    for coordinate_coverage_dictionary_plus in plus_strand_coverage_by_position:

        #only consider positions within the genome, because the dictionary extends beyond the end of the genome
        if coordinate_coverage_dictionary_plus <= len_of_genome:

            #coverage must match or exceed threshold for a position to be called as a peak
            if plus_strand_coverage_by_position[coordinate_coverage_dictionary_plus] >= threshold:

                #list of coverage scores for every position +/- span from coordinate
                span_values = []

                #span_length is the length of the span tested; almost always 100, except in cases where the coordinate is at positions <100
                #span_length starts as -1 to account for the position being tested, which is not included in the span calculation
                span_length = -1

                #defines a range around the position of interest, starting span nt upstream (or position 0 if 0 is higher) and ending span nt downstream
                for span_around_coordinate in range(max(coordinate_coverage_dictionary_plus - span, 0), coordinate_coverage_dictionary_plus + span + 1):

                    #only consider positions with coverage, i.e. those represented in the dictionary
                    if span_around_coordinate in plus_strand_coverage_by_position:

                        #determine the normalized sequence read coverage at every position in the span range
                        span_values.append(plus_strand_coverage_by_position[span_around_coordinate])

                    #keep track of the size of the span range, which will almost always be 101 minus 1 for the position being tested
                    span_length += 1

                #calculate the mean coverage value in the span range
                span_value_mean = (float(sum(span_values)) - plus_strand_coverage_by_position[coordinate_coverage_dictionary_plus])/span_length #average all the coverage values in span excluding coverage at the position being tested 

                #list of coverage scores for every position +/- window from coordinate
                window_values = []

                #defines a range around the position of interest, starting window nt upstream (or position 0 if 0 is higher) and ending window nt downstream
                for window_around_coordinate in range(max(coordinate_coverage_dictionary_plus - window, 0), coordinate_coverage_dictionary_plus + window + 1):

                    #only consider positions with coverage, i.e. those represented in the dictionary
                    if window_around_coordinate in plus_strand_coverage_by_position:

                        #determine the normalized sequence read coverage at every position in the window range
                        window_values.append(plus_strand_coverage_by_position[window_around_coordinate])

                #add the position being tested to plus_coord_selected if
                #it is at least min_fold times higher than the average value in the span range
                #and it is the highest or joint-highest in the window range
                if plus_strand_coverage_by_position[coordinate_coverage_dictionary_plus] >= min_fold * span_value_mean:
                    if plus_strand_coverage_by_position[coordinate_coverage_dictionary_plus] == max(window_values):
                        plus_coord_selected.append(coordinate_coverage_dictionary_plus)


    
    #repeat with data from the minus strand
    for coordinate_coverage_dictionary_minus in minus_strand_coverage_by_position:
        if coordinate_coverage_dictionary_plus <= len_of_genome:
            if minus_strand_coverage_by_position[coordinate_coverage_dictionary_minus] >= threshold:
                span_values = []
                span_length = -1
                for span_around_coordinate in range(max(coordinate_coverage_dictionary_minus - span, 0), coordinate_coverage_dictionary_minus + span + 1):
                    if span_around_coordinate in minus_strand_coverage_by_position:
                        span_values.append(minus_strand_coverage_by_position[span_around_coordinate])
                    span_length += 1
                span_value_mean = (float(sum(span_values)) - minus_strand_coverage_by_position[coordinate_coverage_dictionary_minus])/span_length #average all the coverage values in span excluding coverage at the position being tested 
                window_values = []
                for window_around_coordinate in range(max(coordinate_coverage_dictionary_minus - window, 0), coordinate_coverage_dictionary_minus + window + 1):
                    if window_around_coordinate in minus_strand_coverage_by_position:
                        window_values.append(minus_strand_coverage_by_position[window_around_coordinate])
                if minus_strand_coverage_by_position[coordinate_coverage_dictionary_minus] >= min_fold * span_value_mean:
                    if minus_strand_coverage_by_position[coordinate_coverage_dictionary_minus] == max(window_values):
                        minus_coord_selected.append(coordinate_coverage_dictionary_minus)
       

    return(plus_strand_coverage_by_position, plus_coord_selected, minus_strand_coverage_by_position, minus_coord_selected)




###################
# codon_frequency #
###################

#takes as input:
#(i) genome sequence
#(ii) revcomp genome sequence
#(iii) genome length
#(iv) dictionary of position/coverage for plus strand
#(v) selected positions on the plus strand
#(vi) dictionary of position/coverage for minus strand
#(vii) selected positions on the minus strand
#(viii) filename for writing the output file
#(ix) count, i.e. the cycle the program is on
#codon_frequency is only run on the first cycle, i.e. count = 0

def codon_frequency(gen, revcompgen, len_of_genome, plus_coord_selected, minus_coord_selected, codonfilename):

    #a list of all sequences from -50 to +10 relative to selected positions in the input (plus_coord_selected and minus_coord_selected)
    list_of_sequences_around_peaks = []


    #consider each position in plus_coord_selected
    for position in plus_coord_selected:

        #don't include positions that are within 50 of the start or 10 of the end of the genome
        if len(gen[position - 50 : position + 10]) == 60:

            #make a list that includes position, strand, and sequence from -50 to +10 relative to the position
            list_of_sequences_around_peaks.append([position, '+', gen[position - 50 : position + 10]])

    #repeat for the minus strand
    #positions on the minus strand are counted from the right-hand side of the genome
    for position in minus_coord_selected:
        if len(revcompgen[position - 50 : position + 10]) == 60:
            list_of_sequences_around_peaks.append([len_of_genome + 1 - position, '-', revcompgen[position - 50 : position + 10]])

    
    #a quick way to make a list of all possible codons
    nucleotides = ['A','C','G','T']
    codons = []
    for pos1 in nucleotides:
        for pos2 in nucleotides:
            for pos3 in nucleotides:
                codons.append(pos1 + pos2 + pos3)


    #codon counts for the sum of all sequences in list_of_sequences_around_peaks
    list_of_codon_counts = []

    #loop through a range that is as long as each sequence, minus 2, since trinucleotides will be considered
    #each trinucleotide will be considered in turn, looking at all sequences before moving to the next position downstream
    for position in range(len(list_of_sequences_around_peaks[0][2]) - 2):

        #start out with a list of zeroes corresponding to each of the 64 trinucleotide sequences
        codon_count = [0 for i in range(len(codons))]

        #consider each sequence in turn
        for sequence in list_of_sequences_around_peaks:

            #extract the trinucleotide sequence at the defined position for the defined sequence
            trinucleotide = sequence[2][position : position + 3]
            
            #add 1 to the position in codon_count that corresponds to the sequence of trinucleotide
            codon_count[codons.index(trinucleotide)] += 1

        #each "row" of list_of_codon_counts is a codon_count list for a trinucleotide position
        #e.g. the first row has codon_count values for the trinucleotide starting at position -50
        list_of_codon_counts.append(codon_count)


    #write list_of_codon_counts to a file
    with open(codonfilename + '_' + date + time + '.txt', 'w') as grid:

        #for each trinucleotide sequence
        for codon in range(len(codons)):

            #start out by writing a row with each trinucleotide sequence
            grid.write(str(codons[codon]))

            #consider each row in codon_counts, where each row corresponds to a different trinucleotide position
            for codon_count in range(len(list_of_codon_counts)):

                #write the codon count for that codon and position
                grid.write('\t' + str(list_of_codon_counts[codon_count][codon]))
            grid.write('\n')

    return(list_of_sequences_around_peaks, list_of_codon_counts)





##################
# MFE statistics #
##################

#takes final lists of annotated/novel/isoform ORFs
#position_column specifies the column with the coordinate of interest, e.g. start codon
#strand_column specifies the column with strand
#sequence is extracted from -40 to +20 relative to the start codon
#MFE numbers are calculated by RNAfold
#these numbers are fed to the 'stripplot' function

def MFE_prediction_for_stripplot_input(master_list_of_datasets, position_column, strand_column, upstream_distance, downstream_distance, labels, y_label, log_or_linear, output_filename):

    #write a fasta file with IDs and sequences for each list in master_list_of_datasets

    #to keep track of the fasta files
    list_of_filenames = []
    
    for x in range(len(master_list_of_datasets)):
        f = open(filepath + labels[x] + '_RBS_sequences_' + output_filename + '_' + date + '.txt', 'w')
        list_of_filenames.append('"' + filepath + labels[x] + '_RBS_sequences_' + output_filename + '_' + date + '.txt' + '"')
        for y in master_list_of_datasets[x]:
            if y[strand_column] == '+':
                f.write('>' + str(y[position_column]) + '+\n') #unique ID with coordinate and strand
                f.write(gen[y[position_column] - upstream_distance :y[position_column] + downstream_distance] + '\n')
            elif y[strand_column] == '-':
                f.write('>' + str(y[position_column]) + '-\n') #unique ID with coordinate and strand
                #add 1 to make the range the same as on the plus strand
                f.write(revcomp(gen[y[position_column] - downstream_distance +1:y[position_column] + upstream_distance + 1]) + '\n')
        f.close()


    #add a set of 500 random sequences to the master_list_of_datasets
    f = open(filepath + 'random_sequences_' + output_filename + '_' + date + '.txt', 'w')
    for x in range(250):
        random_coordinate = random.randint(50, len(gen) - 50)
        f.write('>' + str(random_coordinate) + '+\n') #unique ID
        f.write(gen[random_coordinate - upstream_distance: random_coordinate + downstream_distance] + '\n')

        random_coordinate = random.randint(50, len(gen) - 50)
        f.write('>' + str(random_coordinate) + '-\n') #unique ID
        f.write(revcomp(gen[random_coordinate - downstream_distance + 1: random_coordinate + upstream_distance + 1]) + '\n')
    f.close()

    list_of_filenames.append('"' + filepath + 'random_sequences_' + output_filename + '_' + date + '.txt' + '"')
    labels.append('Random')
    
    
    #location of the RNAfold program
    RNAfold_path = '"C:\\Users\\jtw03\\AppData\\Local\\ViennaRNA Package\\RNAfold.exe"'


    #a list of lists where each sublist is a list of dG values
    list_of_list_of_dGs = []
    
    for x in range(len(list_of_filenames)):
        seq_file = list_of_filenames[x]
        out_file_cmd = '--outfile='
        out_file_name = 'RNAfold_out.txt' #temporary, will be deleted at the end
        RNAfold_output_file = filepath + labels[x] + '_MFE_deltaG_scores_' + output_filename + '.txt'
        
        #run cmd line for RNAfold then to get default user directory where RNAfold deposits out file
        cmd_list = RNAfold_path + ' -i ' + seq_file + ' ' + out_file_cmd + out_file_name #uses default parameters
        subprocess.run(cmd_list, stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True).stdout.decode('utf8')
        direct = subprocess.run(['dir'], stdout=subprocess.PIPE, shell=True).stdout.decode('utf8').split('\n')
        path = PureWindowsPath(direct[3].rstrip()[14:])
        
        dG_list = []
        
        #get dG from RNAfold out file
        with open(str(Path(path)) + '\\' + out_file_name, 'r') as r, open(RNAfold_output_file, 'w') as w:
            for line in r:
                ID = line.rstrip()[1:]
                seq = next(r).rstrip()
                fold_line = next(r).rstrip()
                dG = float(fold_line.split(' (')[1][0:-1].lstrip(' '))
                w.write(str(dG) + '\n')
                dG_list.append(dG)
        
        #delete postscript and temporary RNAout files
        user_files = os.listdir(str(Path(path)))
        ps_files = [f for f in user_files if f.endswith("_ss.ps")]
        for f in ps_files:
            path_to_file = os.path.join(str(Path(path)), f)
            os.remove(path_to_file)
        os.remove(str(Path(path)) + '\\' + out_file_name)
        
        list_of_list_of_dGs.append(dG_list)


    #Mann-Whitney U test descriptions and p-values
    statistical_comparisons = []
    for x in range(len(list_of_list_of_dGs)-1):
        for y in range(x + 1, len(list_of_list_of_dGs)):
            statistical_comparisons.append(['MFE ' + labels[x] + ' vs ' + labels[y], scipy.stats.mannwhitneyu(list_of_list_of_dGs[x], list_of_list_of_dGs[y], True, 'greater')])
            statistical_comparisons.append(['MFE ' + labels[x] + ' vs ' + labels[y], scipy.stats.mannwhitneyu(list_of_list_of_dGs[x], list_of_list_of_dGs[y], True, 'less')])

    
    x = stripplot(list_of_list_of_dGs, labels, y_label, log_or_linear, output_filename)

    return statistical_comparisons






##############
# Strip plot #
##############

#makes a stripplot
#list_of_datasets is a list of lists where each sublist is a set of numbers for one column of the stripplot
#labels is a list of labels that is the same length as list_of_datasets; each entry is the label for a dataset
#y_label is the y-axis label
#log_or_linear sets the y-axis scale to be log or linear

def stripplot(list_of_datasets, labels, y_label, log_or_linear, output_filename):


    combined_data = []
    combined_labels = []
    for x in range(len(list_of_datasets)):
        combined_data += list_of_datasets[x]
        for y in list_of_datasets[x]:
            combined_labels.append(labels[x])

    df = pandas.DataFrame()
    df[y_label] = combined_data
    df[' '] = combined_labels
    
    plt.clf()
    custom_stripplot = seaborn.stripplot(data = df, x = ' ', y = y_label, jitter = 0.4, size = 10) #"size" changes datapoint size
    plt.yscale(log_or_linear)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    
    plt.xlabel(xlabel = '', fontsize = 75, labelpad = 40)
    plt.ylabel(ylabel = y_label, fontsize = 75, labelpad = 40)
    plt.gcf().set_size_inches(24,36)
    plt.gcf().subplots_adjust(bottom = 0.25)
    plt.gcf().subplots_adjust(left = 0.25)

    #set the axis label numbers font size
    plt.tick_params(axis='both', which='major', labelsize=75)
    plt.xticks(rotation = 45)
    
    custom_stripplot.figure.savefig(filepath + output_filename + date + '.png')
    
    return None


######################################################################
######################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas
import seaborn
from pathlib import Path, PureWindowsPath
import subprocess
import os
import math
import random
import scipy
#files are saved with date and time information for version control
import datetime
date = str(datetime.datetime.now().month) + "_" + str(datetime.datetime.now().day) + "_" + str(datetime.datetime.now().year)
time = "__" + str(datetime.datetime.now().hour) + "_" + str(datetime.datetime.now().minute)

#set the file path
filepath = '/Users/jtw03/Documents/Projects/Ecoli ribo-seq/April 2021/'

#read in the genome sequence from a .fna file
f = open(filepath + 'CP001509_3.fna', 'r')

#the first line is a header and should be skipped
qqq = f.readline()

#add an 'x' to the start of the genome sequence so that genome position 1 is [1]
gen = 'x' + f.readline().split()[0].upper()
f.close()

#length of the genome; important for some calculations/list manipulations
len_of_genome = 4558952
    
#the program uses Ribo-RET and Ribo-ApiPmn data to identify ORFs
#there are multiple cycles, with the genome sequence and annotation being rotated 1000 bp each cycle
#the first cycle gets the actual numbers of ORF
#later cycles allow an estimate of FDR
#the more cycles, the better the FDR estimate
number_of_cycles = 11

#"peaks" in Ribo-RET or Ribo-ApiPmn data must have higher sequence coverage than the average in a surrounding range
#the size of the range is defined by "span" (number of nt upstream and downstream, e.g. span of 50 = 50 nt upstream to 50 nt downstream)
#the min_fold is the factor by which the peak coverage must be higher than the local average coverage
#note that the peak position itself is excluded from the local average coverage calculation
min_fold = 10
span = 50

#this is an RPM threshold that Ribo-RET or Ribo-ApiPmn data must exceed for positions to be called as peaks
threshold = 0.5

#peaks must have the highest (or joint highest) coverage in a range +/- "window" nt relative to the peak
window = 10



#make lists of + strand and - strand start and stop codon coordinates for the BL21 annotation
#these will be compared with ORFs identified from Ribo-RET/Ribo-ApiPmn data for ORF classification
with open(filepath + 'BL21_genes.gff', 'r') as BL21gff:
    BL21_start_coord_plus = []
    BL21_start_coord_minus = []
    BL21_stop_coord_plus = []
    BL21_stop_coord_minus = []
    for gene in BL21gff:
        if gene.split()[6] == '+':
            BL21_start_coord_plus.append(int(gene.split()[3]))
            BL21_stop_coord_plus.append(int(gene.split()[4]))
        if gene.split()[6] == '-':
            BL21_start_coord_minus.append(int(gene.split()[4]))
            BL21_stop_coord_minus.append(int(gene.split()[3]))

#this list will be written to a file at the end of the program
#it records the following:
#number of Ribo-RET peaks
#number of Ribo-RET peaks associated with a start codon
#number of Ribo-ApiPmn peaks
#number of Ribo-ApiPmn peaks associated with a stop codon after filtering for overlap with Ribo-RET peaks
#number of unique Annotated ORFs
#number of unique Isoform ORFs
#number of unique Novel ORFs
important_values_to_write_to_summary_file = []




    
    ############################################
    # Start of the loop for ORF identification #
    ############################################


#loop through the code as many times as "number_of_cycles"
for count in range(0, number_of_cycles):

    #take the reverse complement of the genome and add 'x' to the beginning so position 1 is [1]
    revcompgen = 'x' + revcomp(gen)


    #run the peak-calling function on Ribo-RET data
    plus_dictionary_ret, plus_coord_selected_ret, minus_dictionary_ret, minus_coord_selected_ret = coord_cov_selected(len_of_genome, threshold, window, span, min_fold, filepath + '4-23-21Ret-RPMgff.txt')




    if count == 0:
        with open(filepath + 'RET all peaks_' + date + time + '.txt', 'w') as ret_all_peaks:
            ret_all_peaks.write('RET Peak Coord\tStrand\tPeak Height\n')

            for plus_coord_entry in plus_coord_selected_ret:
                ret_all_peaks.write(str(plus_coord_entry) + '\t')
                ret_all_peaks.write('+\t')
                ret_all_peaks.write(str(plus_dictionary_ret[plus_coord_entry]) + '\n')

            for minus_coord_entry in minus_coord_selected_ret:
                ret_all_peaks.write(str(len_of_genome + 1 - minus_coord_entry) + '\t')
                ret_all_peaks.write('-\t')
                ret_all_peaks.write(str(minus_dictionary_ret[minus_coord_entry]) + '\n')






    #run the codon frequency function on Ribo-RET data
    #only do this on the first cycle, before the genome sequence has been rotated
    if count == 0:
        list_of_sequences_around_peaks_ret, lists_of_lists_ret = codon_frequency(gen, revcompgen, len_of_genome, plus_coord_selected_ret, minus_coord_selected_ret, filepath + 'ret_codon_frequency')




    
    ##############################################################
    # Look for Ribo-RET peaks with nearby potential start codons #
    ##############################################################


    #allowed start codons
    allowed_start_codons = ['ATG', 'GTG', 'TTG']

    #allowed start codon positions; note that not every combination of start codon and position is allowed; additional filtering occurs later
    allowed_start_gaps = [-18, -17, -16, -15, -14]


    #Ribo-RET peaks associated with suitably positioned start codons will be added to these lists
    start_codon_coord_plus_ret = []
    start_codon_coord_minus_ret = []

    #start_codon_coord_plus_ret is a list of Ribo-RET peaks with appropriate spacing to a potential start codon (to be further refined later)
    #with associated metadata
    #format: [Ribo-RET peak coordinate, coverage, start codon sequence, gap, strand, start codon coordinate]
    for ret_peak_coordinate in plus_coord_selected_ret:
        for start_codon in allowed_start_codons:
            for gap in allowed_start_gaps:
                if gen[ret_peak_coordinate + gap:ret_peak_coordinate + gap + 3] == start_codon:
                    start_codon_coord_plus_ret.append([ret_peak_coordinate, plus_dictionary_ret[ret_peak_coordinate], start_codon, gap, '+', ret_peak_coordinate + gap])

    #as above but for the minus strand
    for ret_peak_coordinate in minus_coord_selected_ret:
          for start_codon in allowed_start_codons:
            for gap in allowed_start_gaps:
                if revcompgen[ret_peak_coordinate + gap:ret_peak_coordinate + gap + 3] == start_codon:
                    #start_codon_coord_minus_ret.append(str((len_of_genome + 1) - m) + '\t' + str(minus_strand_coverage_ret[ret_peak_coordinate]) + '\t' + start_codon + '\t' + str(i) + '\t-\t' + str(int(len_of_genome + 1) - int(ret_peak_coordinate)-i) + '\n')
                    start_codon_coord_minus_ret.append([len_of_genome + 1 - ret_peak_coordinate, minus_dictionary_ret[ret_peak_coordinate], start_codon, gap, '-',  len_of_genome + 1 - ret_peak_coordinate - gap])
	




    #######################################################
    # Associate stop codons with start codons in the list #
    #######################################################


    #allowed stop codons
    allowed_stop_codons = ['TAA', 'TAG', 'TGA']

    #reverse complement of allowed stop codons
    revcomp_stop_codon = ['TTA', 'CTA', 'TCA']

    #reverse complement of allowed start codons
    revcomp_start_codon = ['CAT', 'CAC', 'CAA']

    #coordinates of stop codon linked to start codons associated with Ribo-RET peaks
    stop_codon_coord_ret = []

    #nucleotide sequences of ORFs linked to start codons associated with Ribo-RET peaks
    ORF_sequences_ret = []

    #reverse complement of nucleotide sequences of ORFs linked to start codons associated with Ribo-RET peaks
    ORF_sequencesrevcomp_ret = []




    #for every start codon, translate out to the end and add details of the stop codon and ORF/protein to the list
    for start_codon_list_entry in start_codon_coord_plus_ret:

        #start codon position; easier than writing it out every time
        start = start_codon_list_entry[5]
        
        #move across the genome in steps of 3 nt
        for genome_position in range(start, len_of_genome+1000,3):

            #prevents running off the end of the genome
            try:
                trinucleotide = gen[genome_position : genome_position + 3]
            except ValueError:
                print("error")
                continue

            #check to see if you have reached a stop codon
            if trinucleotide in allowed_stop_codons:

                #append [stop codon sequence, stop codon coordinate, ORF nucleotide sequence, protein sequence]
                #frames is the output of the ORFtranslator function, which returns protein sequences from all 3 frames of a nucleotide sequence input
                start_codon_list_entry.extend((trinucleotide, genome_position + 2, gen[start : genome_position + 3], ORFtranslator(gen[start : genome_position + 3])[0]))

                #break when you hit a stop codon
                break   




    #repeat on the minus strand

    for start_codon_list_entry in start_codon_coord_minus_ret:
        
        start = start_codon_list_entry[5]

        #to account for the minus strant
        revcompstart = (len_of_genome + 1) - start + 2
        
        for genome_position in range(revcompstart - 2, len_of_genome + 1000, 3):
            
            try:
                trinucleotide = revcompgen[genome_position : genome_position + 3]
            except IndexError:
                print("errorminus")
                continue
            
            if trinucleotide in allowed_stop_codons:
                start_codon_list_entry.extend((trinucleotide, len_of_genome + 1 - genome_position - 2, revcompgen[revcompstart - 2: genome_position + 3], ORFtranslator(revcompgen[revcompstart - 2 : genome_position + 3])[0]))

                break





    #############################################################################
    # Additional filtering of the ORF list for start codon position and spacing #
    #############################################################################


    combined_start_codon_list = start_codon_coord_plus_ret + start_codon_coord_minus_ret

    final_start_codon_list = []
    discarded_start_codons = []

    acceptable_start_codon_sequences_and_gaps = [['-14','ATG'],['-15','ATG'],['-16','ATG'],['-17','ATG'],['-18','ATG'],['-14','GTG'],['-15','GTG'],['-16','GTG'],['-14','TTG'],['-15','TTG']]


    for potential_ORF in combined_start_codon_list:
        for codon_gap_combination in acceptable_start_codon_sequences_and_gaps:
            if str(potential_ORF[3]) == codon_gap_combination[0] and potential_ORF[2] == codon_gap_combination[1]:
                final_start_codon_list.append(potential_ORF)

    for potential_ORF in combined_start_codon_list:
        if potential_ORF not in final_start_codon_list:
            discarded_start_codons.append(potential_ORF)








    if count == 0:
        with open(filepath + 'Ret-full report_' + date + time + '.txt', 'w') as ret_full:
            ret_full.write('Ret Peak Coord\tRet Peak Height\tStart Codon\tGap\tStrand\tStart Codon Coord\tStop Codon\tStop Codon Coord\tnt Seq\taa Seq\n')

            for entry in final_start_codon_list:
                for item in entry:
                    ret_full.write(str(item) + '\t')

                ret_full.write('\n')






    #####################################
    # Peak calling on Ribo-Api/Pmn data #
    #####################################


    #runs the peak-calling function on Ribo-Api/Pmn data
    plus_dictionary_api, plus_coord_selected_api_before_filtering, minus_dictionary_api, minus_coord_selected_api_before_filtering = coord_cov_selected(len_of_genome,threshold, window, span, min_fold, filepath + 'ApipmnRPMgff4_23_21_renamed.gff')


    if count == 0:
        with open(filepath + 'Api all peaks_' + date + time + '.txt', 'w') as api_all_peaks:
            api_all_peaks.write('Api/Pmn Peak Coord\tStrand\tPeak Height\n')

            for plus_coord_entry in plus_coord_selected_api_before_filtering:
                api_all_peaks.write(str(plus_coord_entry) + '\t')
                api_all_peaks.write('+\t')
                api_all_peaks.write(str(plus_dictionary_api[plus_coord_entry]) + '\n')

            for minus_coord_entry in minus_coord_selected_api_before_filtering:
                api_all_peaks.write(str(len_of_genome + 1 - minus_coord_entry) + '\t')
                api_all_peaks.write('-\t')
                api_all_peaks.write(str(minus_dictionary_api[minus_coord_entry]) + '\n')


    #Codon Frequency Analysis for Ribo-Api/Pmn peaks
    if count == 0:
        list_of_sequences_around_peaks_api, lists_of_lists_api = codon_frequency(gen, revcompgen, len_of_genome, plus_coord_selected_api_before_filtering, minus_coord_selected_api_before_filtering, filepath + 'api_codon_frequency')





    ######################################################
    # Remove Ribo-ApiPmn peaks that match Ribo-RET peaks #
    ######################################################


    #these lists will have removed Ribo-Api/Pmn peaks that are within 3 nt of a Ribo-RET peak
    plus_coord_selected_api = []
    minus_coord_selected_api = []

    #compare Ribo-ApiPmn peaks to Ribo-RET peaks and discard Ribo-Api/Pmn peaks within 3 nt of any Ribo-RET peak
    for Ribo_Api_peak_coordinate in plus_coord_selected_api_before_filtering:
        Toggle = True
        for Ribo_RET_peak_coordinate in plus_coord_selected_ret:
            if -3 <= Ribo_Api_peak_coordinate - Ribo_RET_peak_coordinate <= 3:
                Toggle = False
                break
        if Toggle == True:
            plus_coord_selected_api.append(Ribo_Api_peak_coordinate)

    #repeat for minus strand
    for Ribo_Api_peak_coordinate in minus_coord_selected_api_before_filtering:
        Toggle = True
        for Ribo_RET_peak_coordinate in minus_coord_selected_ret:
            if -3 <= Ribo_Api_peak_coordinate - Ribo_RET_peak_coordinate <= 3:
                Toggle = False
                break
        if Toggle == True:
            minus_coord_selected_api.append(Ribo_Api_peak_coordinate)


    
    ##############################################################
    # Look for Ribo-RET peaks with nearby potential start codons #
    ##############################################################



    allowed_stop_codons = ['TAA', 'TAG', 'TGA']
    
    allowed_stop_gaps = [-14, -13, -12]
    
    #Ribo-ApiPmn peaks associated with suitably positioned start codons will be added to these lists
    stop_codon_coord_plus_api = []
    stop_codon_coord_minus_api = []


    #stop_codon_coord_plus_api is a list of Ribo-ApiPmn peaks with appropriate spacing to a potential stop codon (to be further refined later)
    #with associated metadata
    #format: [Ribo-ApiPmn peak coordinate, coverage, stop codon sequence, gap, strand, stop codon coordinate]
    for api_peak_coordinate in plus_coord_selected_api:
        for stop_codon in allowed_stop_codons:
            for gap in allowed_stop_gaps:
                if gen[api_peak_coordinate + gap : api_peak_coordinate + gap + 3] == stop_codon:
                    stop_codon_coord_plus_api.append([api_peak_coordinate, plus_dictionary_api[api_peak_coordinate], stop_codon, gap, '+', api_peak_coordinate + gap + 2])

    #as above but for the minus strand
    for api_peak_coordinate in minus_coord_selected_api:
          for stop_codon in allowed_stop_codons:
            for gap in allowed_stop_gaps:
                if revcompgen[api_peak_coordinate + gap : api_peak_coordinate + gap + 3] == stop_codon:
                    stop_codon_coord_minus_api.append([len_of_genome + 1 - api_peak_coordinate, minus_dictionary_api[api_peak_coordinate], stop_codon, gap, '-', len_of_genome + 1 - api_peak_coordinate - gap - 2])


    #stop_codon_coord_api = stop_codon_coord_plus_api + stop_codon_coord_minus_api
    combined_stop_codon_list = stop_codon_coord_plus_api + stop_codon_coord_minus_api






    ####################################################################################
    # Additional filtering of the stop codon list for start codon position and spacing #
    ####################################################################################


    final_stop_codon_list = []
    discarded_stop_codons = []

    acceptable_stop_codon_sequences_and_gaps = [['-14','TAA'],['-13','TAA'],['-12','TAA'],['-14','TAG'],['-13','TAG'],['-12','TAG'],['-14','TGA'],['-13','TGA'],['-12','TGA']]


    for potential_stop in combined_stop_codon_list:
        for codon_gap_combination in acceptable_stop_codon_sequences_and_gaps:
            if str(potential_stop[3]) == codon_gap_combination[0] and potential_stop[2] == codon_gap_combination[1]:
                final_stop_codon_list.append(potential_stop)

    for potential_stop in combined_stop_codon_list:
        if potential_stop not in final_stop_codon_list:
            discarded_stop_codons.append(potential_stop)


    if count == 0:
        with open(filepath + 'Api-full report_' + date + time + '.txt', 'w') as api_full:
            api_full.write('ApiPmn Peak Coord\tPeak Height\tStop Codon\tGap\tStrand\tStop Codon Coord\n')

            for entry in final_stop_codon_list:
                for item in entry:
                    api_full.write(str(item) + '\t')

                api_full.write('\n')





    ####################################################
    # Select ORFs with Ribo-RET and Ribo-ApiPmn signal #
    ####################################################


    #list of ORFs with evidence from both Ribo-RET and Ribo-ApiPmn
    ORFs_found_with_Ret_and_Api = []
    
    for ORF in final_start_codon_list:
        for stop_codon in final_stop_codon_list:
            
            #Checks that the length of the ORF entry is 10, which it wouldn't be in very rare cases where the ORF extended beyond the genome end
            if len(ORF) == 10:

                #Stop codon coordinate is position 7 and 5, for the start/stop codon lists, respectively
                #Position 4 is the strand for both lists
                if ORF[7] == stop_codon[5] and ORF[4] == stop_codon[4]:
                    
                    #Add the Ribo-ApiPmn details
                    #Create a new list (extended_ORF) to account for cases where there are multiple Ribo-ApiPmn matches to the same start codon
                    extended_ORF = ORF + [stop_codon[0], stop_codon[1], stop_codon[3]]
                    ORFs_found_with_Ret_and_Api.append(extended_ORF)





    ###################################################
    # Classify ORFs based on comparison to annotation #
    ###################################################


    #Annotated ORFs are perfect matches to the existing annotation
    Annotated_ORFs = []

    #Isoform ORFs have a stop codon that matches a stop codon in the existing annotation, but start codons do not match
    Isoform_ORFs = []

    #Novel ORFs have start and stop codons that do not match the existing annotation
    Novel_ORFs = []


    #Go through the list of ORFs that have Ribo-RET and Ribo-ApiPmn peaks
    for ORF in ORFs_found_with_Ret_and_Api:

        #if Toggle = True then the ORF has been classified, and should not be assessed further
        #test here whether the ORF is an Annotated ORF by comparing its start codon coordinate to the BL21 annotation
        #this must be done separately for + and - strands in case there are matching coordinates for both strands
        Toggle = False
        if ORF[4] == '+':
            for BL21_start in BL21_start_coord_plus:
                if ORF[5] == BL21_start:
                    Annotated_ORFs.append(ORF)

                    #if an ORF is classified as Annotated then there is no need to assess whether it is an Isoform or Novel ORF
                    Toggle = True


            #test here whether the ORF is an Isoform ORF by comparing its stop codon coordinate to the BL21 annotation
            #Toggle is only False if the ORF is not an Annotated ORF, meaning it could be an Isoform or a Novel ORF
            if Toggle == False:
                for BL21_stop in BL21_stop_coord_plus:
                    if ORF[7] == BL21_stop:
                        Isoform_ORFs.append(ORF)

                        #if an ORF is classified as Isoform then there is no need to assess whether it is a Novel ORF
                        Toggle = True


            #Toggle is only False if the ORF is neither Annotated nor Isoform, meaning it must be Novel
            if Toggle == False:
                Novel_ORFs.append(ORF)


    #repeat for the minus strand
    for ORF in ORFs_found_with_Ret_and_Api:
        Toggle = False
        if ORF[4] == '-':
            for BL21_start in BL21_start_coord_minus:
                if ORF[5] == BL21_start:
                    Annotated_ORFs.append(ORF)
                    Toggle = True

            if Toggle == False:
                for BL21_stop in BL21_stop_coord_minus:
                   if ORF[7] == BL21_stop:
                        Isoform_ORFs.append(ORF)
                        Toggle = True

            if Toggle == False:
                Novel_ORFs.append(ORF)


    ####################################################
    # Write the ORFs and associated metadata to a file #
    ####################################################

    #only write the file for the first cycle
    if count == 0:
        
        ORF_list = open(filepath + 'ORF_list_' + date + time + '.txt', 'w')

        #add column descriptions
        ORF_list.write('Ribo-RET peak coord\tRibo-RET peak height\tStart codon\tStart codon to Ribo-RET peak distance\tStrand\t')
        ORF_list.write('Start codon coord\tStop codon\tStop codon coord\tnt seq\taa seq\t')
        ORF_list.write('Ribo-ApiPmn peak coord\tRibo-ApiPmn peak height\tStart codon to Ribo-ApiPmn peak distance\tCategory\n')
        
        for entry in Annotated_ORFs:
            for item in entry:
                ORF_list.write(str(item) + '\t')
                
            #add the classification as a final column
            ORF_list.write('Annotated\n')

        for entry in Isoform_ORFs:
            for item in entry:
                ORF_list.write(str(item) + '\t')
            ORF_list.write('Isoform\n')

        for entry in Novel_ORFs:
            for item in entry:
                ORF_list.write(str(item) + '\t')
            ORF_list.write('Novel\n')

        ORF_list.close()
    
    

    ORFsum = 0







    ##############################################
    # Analyze RNA structure and make strip plots #
    ##############################################

    if count == 0:
        MFE_statistical_comaprison_starts = MFE_prediction_for_stripplot_input([Annotated_ORFs, Isoform_ORFs, Novel_ORFs], 5, 4, 25, 16, ['Annotated', 'Isoform', 'Novel'], 'Delta G', 'linear', 'Annotated_Isoform_Novel_MFE_start_codon_scores')

        MFE_starts = open(filepath + 'MFE_statistical_summary_starts_' + date + time + '.txt', 'w')
        for entry in MFE_statistical_comaprison_starts:
            for item in entry:
                MFE_starts.write(str(item) + '\t')
            MFE_starts.write('\n')
        MFE_starts.close()


    
    ###################################################################
    # Count the number of unique ORFs in each of the three categories #
    ###################################################################


    #lists of ORFs from each category that just include the strand and start coordinate
    #these will be used to make sets (i.e. unique entries)
    #the lengths of the sets will indicate the number of unique ORFs in each category
    #this accounts for Ribo-RET and/or Ribo-ApiPmn peaks that are close together and lead to the repeated identification of the same ORF
    Annotated_ORFs_strand_and_start = []
    Isoform_ORFs_strand_and_start = []
    Novel_ORFs_strand_and_start = []
    
    for ORF in Annotated_ORFs:
        #make a string that consists of ORF strand and start codon coordinate
        #a list of strings can be converted to a set (a list of lists cannot)
        Annotated_ORFs_strand_and_start.append(ORF[4] + str(ORF[5]))

    for ORF in Isoform_ORFs:
        Isoform_ORFs_strand_and_start.append(ORF[4] + str(ORF[5]))

    for ORF in Novel_ORFs:
        Novel_ORFs_strand_and_start.append(ORF[4] + str(ORF[5]))

    #calculate the number of unique entries in each of the lists of start codon strands/coordinates
    number_of_unique_Annotated_ORFs = len(set(Annotated_ORFs_strand_and_start))
    number_of_unique_Isoform_ORFs = len(set(Isoform_ORFs_strand_and_start))
    number_of_unique_Novel_ORFs = len(set(Novel_ORFs_strand_and_start))
    


    
    ######################################################################################################################
    # Add information about peak counts and ORF counts to a list that will be written to a file at the end of the program#
    ######################################################################################################################


    #this could all be added to the master list in one go, but it's a lot of text, so easier to read if we build it up piece by piece first
    list_of_values_for_this_cycle = []

    #number of Ribo-RET peaks
    list_of_values_for_this_cycle.append(len(plus_coord_selected_ret) + len(minus_coord_selected_ret))

    #number of Ribo-RET peaks associated with a start codon
    list_of_values_for_this_cycle.append(len(final_start_codon_list))

    #number of Ribo-ApiPmn peaks
    list_of_values_for_this_cycle.append(len(plus_coord_selected_api) + len(minus_coord_selected_api))

    #number of Ribo-ApiPmn peaks associated with a stop codon after filtering for overlap with Ribo-RET peaks
    list_of_values_for_this_cycle.append(len(final_stop_codon_list))

    #number of unique Annotated ORFs
    list_of_values_for_this_cycle.append(number_of_unique_Annotated_ORFs)

    #number of unique Isoform ORFs
    list_of_values_for_this_cycle.append(number_of_unique_Isoform_ORFs)

    #number of unique Novel ORFs
    list_of_values_for_this_cycle.append(number_of_unique_Novel_ORFs)


    important_values_to_write_to_summary_file.append(list_of_values_for_this_cycle)



    
    ###########################################################################
    # Rotate the genome and the BL21 annotation by 1000 bp for the next cycle #
    ###########################################################################


    gen = 'x' + gen[-1000:] + gen[1:-1000]


    for start_codon_coord in range(len(BL21_start_coord_plus)):
        if BL21_start_coord_plus[start_codon_coord] < len_of_genome - 1000:
            BL21_start_coord_plus[start_codon_coord] += 1000
        else:
            BL21_start_coord_plus[start_codon_coord] += 1000 - len_of_genome

    for start_codon_coord in range(len(BL21_start_coord_minus)):
        if BL21_start_coord_minus[start_codon_coord] < len_of_genome - 1000:
            BL21_start_coord_minus[start_codon_coord] += 1000
        else:
            BL21_start_coord_minus[start_codon_coord] += 1000 - len_of_genome

    for stop_codon_coord in range(len(BL21_stop_coord_plus)):
        if BL21_stop_coord_plus[stop_codon_coord] < len_of_genome - 1000:
            BL21_stop_coord_plus[stop_codon_coord] += 1000
        else:
            BL21_stop_coord_plus[stop_codon_coord] += 1000 - len_of_genome

    for stop_codon_coord in range(len(BL21_stop_coord_minus)):
        if BL21_stop_coord_minus[stop_codon_coord] < len_of_genome - 1000:
            BL21_stop_coord_minus[stop_codon_coord] += 1000
        else:
            BL21_stop_coord_minus[stop_codon_coord] += 1000 - len_of_genome


    #########################################################################
    # End of the loop; repeat this cycle multiple times to estimate the FDR #
    #########################################################################

    print(count)       



##########################
# Write the summary file #
##########################


with open(filepath + 'Summary_' + date + '.txt', 'w') as summary:

    summary.write(str(datetime.datetime.now() ) + '\n')

    summary.write('Threshold = ' + str(threshold) + '\n')

    summary.write('Window (must be local maximum this distance upstream/downstream) = ' + str(window) + '\n')

    summary.write('Span (must be x-fold above the mean value at positions this distance upstrea/downstream) = ' + str(span) + '\n')

    summary.write('Fold over Span mean threshold = ' + str(min_fold) + '\n')

    summary.write('Acceptable start codon sequences and gaps \t')

    for combination in acceptable_start_codon_sequences_and_gaps:
        summary.write(str(combination) + '\t')

    summary.write('\n')
    
    summary.write('Acceptable stop codon sequences and gaps \t')

    for combination in acceptable_stop_codon_sequences_and_gaps:
        summary.write(str(combination) + '\t')

    summary.write('\n\n')

    summary.write('number of Ribo-RET peaks\tnumber of Ribo-RET peaks associated with a start codon\tnumber of Ribo-ApiPmn peaks\tnumber of Ribo-ApiPmn peaks associated with a stop codon after filtering for overlap with Ribo-RET peaks\tnumber of unique Annotated ORFs\tnumber of unique Isoform ORFs\tnumber of unique Novel ORFs\n')

    for entry in important_values_to_write_to_summary_file:
        for item in entry:
            summary.write(str(item) + '\t')
        summary.write('\n')




    
