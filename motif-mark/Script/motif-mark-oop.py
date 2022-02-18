#!/usr/bin/env python
from ast import pattern
from itertools import product
import re
import cairo
import math
from IPython import display
from IPython.display import SVG, display, Image
import re


IPUAC_dict = {'A':['A'],'C':['C'], 'G':['G'],'T':['T'],'U':['U'],
                'W':['A','T'],'S':['G','C'],'M':['A','C'],'K':['G','T'],
                'R':['A','G'],'Y':['C','T'], 'B':['C','G','T'],
                'D':['A','G','T'],'H':['A','C','T'],'V':['A','C','G'],
                'N':['A','C','T','G']}


class Sequence:
    def __init__(self, fasta_file_name):
        self.file = fasta_file_name
        self.sequence = self.oneline_fasta(fasta_file_name)
        self.intron1, self.exon, self.intron2 = self.find_exon_introns(self.sequence)

    def oneline_fasta(self, file):
        '''
        Takes a Fasta file where the sequences are split across multiple lines, and concatenates the sequence lines to a single string. Returns the sequence line.
        Input: fasta_file_handle
        Output: sequence
        '''
        file_handle = open(file, "rt")
        all_lines = file_handle.read()
        record_array = all_lines.split('>')
        record_list = []
        for record in record_array[1:]:
            record_line_array = record.split('\n')  
            seq = record_line_array[1:]
            seq_line = ('').join(seq)
            return(seq_line)
    
    def find_exon_introns(self, sequence):
        '''
        This function takes as its input a DNA sequence with introns lowercase and exons UPPERCASE. It identifies the positions of the first intron, the exon
        and the second intron.
        Input: tttgtttccacagAAAAACCTCTTCAGGCACTGGTGCCGAGGACCCTAGgtatgactcacctgt"
        Output: ([0, 12], [13, 48], [49, 63])
        '''
        intron = re.compile("[actgn]+")
        exon = re.compile("[ACTGN]+")
        exon_positions =  [[m.start(),m.end()] for m in exon.finditer(sequence)]
        intron_positions = [[m.start(),m.end()] for m in intron.finditer(sequence)]
        intron1_start_end = [intron_positions[0][0], (intron_positions[0][1]-1)] 
        exon_start_end = [exon_positions[0][0], (exon_positions[0][1]-1)]
        intron2_start_end = [intron_positions[1][0], (intron_positions[1][1]-1)] 
        return intron1_start_end, exon_start_end, intron2_start_end    


class Motifs:
    def __init__(self, motif_file_name):
        # things that only have to happen one time 
        self.file = motif_file_name
        self.motif_array = self.init_motif_array(motif_file_name)
        self.motif_options_dict = self.find_all_motif_options()
        
    def __str__(self):
        return str(self.motif_array)

    def init_motif_array(self, file):
        '''
        Create an array from the motifs in the file 
        Input: motif_file_handle
        Output: [motif1, motif2, motif3]
        '''
        motif_file = open(self.file, "rt")
        all_motif_lines = motif_file.read()
        motif_array = all_motif_lines.split('\n')
        motif_file.close()
        return motif_array

    def get_motif_array(self):
        return self.motif_array
        
    def find_all_motif_options(self):
        '''
        Takes as input a single ambiguous motif and generates all possible motifs from it. Returns a dictionary with {motif: [[motif1.1, motif 1.2, motif 1.3],[motif2.1, motif 2.2, motif 2.3]}
        Input: [motif1, motif2, motif3]
        Output: [[motif1.1, motif 1.2, motif 1.3],[motif2.1, motif 2.2, motif 2.3], ...]
        '''
        motif_dict = {}
        for motif in self.motif_array:
            motif = motif.upper() # make all letters uppercase
            all_characters = []
            for character in motif:
                one_position_possible_chars = IPUAC_dict[character] # create an array for each position in the motif that holds every possible character that could be at that position according to its IUPAC
                all_characters.append(one_position_possible_chars)

            possible_strings_iters = list(product(*all_characters)) # create strings of all possible cut sites
            possible_strings = []
            for i in range(0, len((possible_strings_iters))):
                possible_strings.append(''.join(possible_strings_iters[i]))   
            motif_dict[motif] = possible_strings
        return motif_dict

    def get_all_motif_options(self):
        return self.motif_options_dict 


motifs = Motifs("Fig_1_motifs.txt")
sequence = Sequence("Figure_1.fasta")
seq = sequence.oneline_fasta("Figure_1.fasta")

    
#def Find_Motifs():
    # look at all motifs in the dictionary - 
motif_dict = motifs.get_all_motif_options()

motif_occurence_dict = {}
for motif in motif_dict.keys(): # for each motif
    motif_occurence_dict[motif] = []
    for possible_motif in motif_dict[motif]: # for each motif option (because motifs are ambiguous)
        motif_pattern_exon = re.compile(possible_motif)
        motif_patter_intron = re.compile(possible_motif.lower())
        motif_pattern_intron_occurence = [[m.start(),m.end()] for m in motif_patter_intron.finditer(seq)]
        if len(motif_pattern_intron_occurence) != 0:
            motif_occurence_dict[motif].append(motif_pattern_intron_occurence)
print(motif_occurence_dict)

def show_svg(file):
    display(SVG(filename=file))

#create the coordinates to display your graph
width, height = 1000, 500

surface = cairo.SVGSurface("line_and_rectangle.png", width, height) # create png with w/h dimensions
context = cairo.Context(surface) # create the surface on which to draw
context.set_line_width(3)
context.move_to(100,25) # origin of the line, left-most point
context.line_to(500,25) # right-most point
context.stroke() 

# draw a rectangle
context.rectangle(300,13,140,25) # (x, y, width, height)
context.fill() 
surface.write_to_png("line_and_rectangle.png")



