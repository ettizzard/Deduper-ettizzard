#!/usr/bin/env python

#Importing necessary modules.
import argparse
import re

#Initializing argparse.
def get_args():
    parser = argparse.ArgumentParser(prog = "Deduper", description = "A PCR duplicate removal tool for single-end sequencing data which requires both a sorted SAM file as well as the list of Unique Molecular Identifiers used for sequencing as input. If duplicates are encountered, only the first copy in the input SAM file is saved in the output file.", epilog = "To achieve desired deduplication, please ensure that input files are formatted as instructed. Please note that BAM files are not acceptable, and they must be sorted and converted to SAM before using this tool.")
    parser.add_argument("-f", "--file", help = "Input SAM file sorted by chromosome (default 'samtools sort' method)", type = str)
    parser.add_argument("-o", "--outfile", help = "Name of output deduplicated SAM file", type = str)
    parser.add_argument("-u", "--umi",  help = "Input text file containing UMIs - formatted as 1 unique UMI sequence per line", type = str)
    return parser.parse_args()
args = get_args()


#Creating set of UMIs to ensure each alignment is properly identified
set_of_UMIs = set()
with open(args.umi, "r") as UMI_file:
    for line in UMI_file:
        line = line.strip()
        set_of_UMIs.add(line)
#print(set_of_UMIs)


#Function to determine which strand the alignment is on given the bitwise flag
def determine_strand(bitwise_flag: int) -> str:
    '''Receives land interprets the bitwise flag, then outputs strand orientation character "+" or "-".'''
    
    #Slicing alignment line by tab to isolate bitwise flag in 2nd column
    #bitwise_flag = int(alignment.split('\t')[1])
    
    if ((bitwise_flag & 16) == 16):
        strand = "-"
    
    else:
        strand = "+"

    return strand 

#Setting regular expression for CIGAR string as a variable to contain each operator and its corresponding count within a tuple, all
#tuples contained within a single list. e.g. [("2", "S"), ("3", "M"), ...]
cigar_string_regex_tuple = re.compile(r'(\d+)(\w)')



#Function to calculate the "true" starting position of the read when soft-clipping occurs
def determine_true_starting_position(cigar_string, given_position, strand: str):
    '''Receives alignment line and returns true starting position of read'''
    
    sum_softclipped_bases_neg_strand = 0
    
    #Positive strand logic
    if strand == "+":
        
        #If strand is + and soft-clipping occurs, the true starting position is the position given by the SAM file minus the number of clipped bases
        #on the left side of the read.
        if "S" in cigar_string:
            
            num_left_softclipped_bases_pos_strand = cigar_string.split('S')[0]
            
            #If soft-clipping occurs on both sides of read, we only need to count left-side clipping for the + strand. So, the .isdigit() if statement
            #selects only left-side clipping, because if right-side clipping occurs, the result of .split() would be the rest of the CIGAR string, which
            #would include both numbers and letters (M, D, I, etc.)
            if num_left_softclipped_bases_pos_strand.isdigit():
                true_position = given_position - int(num_left_softclipped_bases_pos_strand) 
                
            else:
                true_position = given_position
        
        #If no soft-clipping occurs on + strand, SAM given position is the true starting position
        else:
            true_position = given_position    
    
    #Negative strand logic
    else:
        #Initializing dictionary of counts of CIGAR operators present in alignment
        dict_cigar_characters = {"D":0, "N":0, "I":0, "S":0, "M":0, "X":0 , "=":0, "H":0, "P":0}
        
        #Utilizing regex pattern to save list of character tuples to "cigar_string_character_tuples_list" variable.
        cigar_string_character_tuples_list = cigar_string_regex_tuple.findall(cigar_string)
        
        #If the negative strand has right-side soft-clipping, slice the list of operator tuples to set count of right-clipped bases
        #to "sum_softclipped_bases_neg_strand" variable.
        if cigar_string[-1] == "S":
            sum_softclipped_bases_neg_strand = int(cigar_string_character_tuples_list[-1][0])
        
        #For every operator tuple in the CIGAR string, set dictionary key to operator character and set value to count
        for operator_count_pair in cigar_string_character_tuples_list:
            dict_cigar_characters[operator_count_pair[1]] = dict_cigar_characters[operator_count_pair[1]] + int(operator_count_pair[0])

        #Negative strand true starting position is the sum of the SAM given position (minus -1 because it makes more sense to me),
        #the number of right-sided soft-clipped bases, the number of matches/mismatches, the number of deletions,
        #and the number of gaps.
        true_position = given_position - 1 + sum_softclipped_bases_neg_strand + dict_cigar_characters["M"] + dict_cigar_characters["D"] + dict_cigar_characters["N"]
 
    return true_position


#Main loop

pos_strand_umi_true_position_set = set()
neg_strand_umi_true_position_set = set()

chromosome_label = "1"
count_of_reads_per_chromosome = 0

with open(args.outfile, "w") as output_sam:

    with open(args.file,"r") as input_sam:

        unknown_UMI_count = 0
        duplicate_count = 0
        unique_read_count = 0
        header_linecount = 0


        for line in input_sam:

            #Header line logic
            if line.startswith("@"):
                
                output_sam.write(line)
                header_linecount += 1

            #Alignment line logic
            else:
                split_line = line.split('\t')
                UMI = split_line[0].split(':')[-1]
                
                #If the current line's chromosome does not match the currently iterated chromosome label, empty the UMI sets for each strand
                if chromosome_label != split_line[2]:
                    pos_strand_umi_true_position_set = set()
                    neg_strand_umi_true_position_set = set()
        
                    print(chromosome_label,count_of_reads_per_chromosome, sep = "\t")
                    chromosome_label = split_line[2]
                    count_of_reads_per_chromosome = 0
                
                #If the current line's UMI is in the overall set of UMIs, determine the strand and true starting position of the read
                if UMI in set_of_UMIs:
                    strand = determine_strand(int(split_line[1]))
                    true_position = determine_true_starting_position(split_line[5], int(split_line[3]), strand)
        
                    #If on the positive strand, 
                    if strand == '+': 
                        
                        #If the current line's UMI and true starting postion pair has not been encountered yet, the read is not duplicated and will be written to file.
                        #The UMI-true position pair will then be added to the set
                        if(UMI,true_position) not in pos_strand_umi_true_position_set:
                            output_sam.write(line)
                            pos_strand_umi_true_position_set.add((UMI,true_position))
                            unique_read_count += 1
                            count_of_reads_per_chromosome += 1
                        
                        #Otherwise, a duplicate has been encountered
                        else:
                            duplicate_count += 1
                    
                    #If on the negative strand, and a duplicate is not encountered, write line to file, add UMI-true position pair to negative set
                    elif (UMI,true_position) not in neg_strand_umi_true_position_set:
                        output_sam.write(line)
                        neg_strand_umi_true_position_set.add((UMI,true_position))
                        unique_read_count += 1
                        count_of_reads_per_chromosome += 1
                    
                    #Otherwise a negative strand duplicate is encountered and counted
                    else:
                        duplicate_count += 1
                
                #If read's UMI is not known from original set, read is not written to file.
                else:
                    unknown_UMI_count += 1   

print("\nHeader lines =", header_linecount, "\nUnique reads =", unique_read_count, "\nNumber of unknown UMIs =", unknown_UMI_count, "\nNumber of duplicates =", duplicate_count)      
     