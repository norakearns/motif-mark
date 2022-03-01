#!/usr/bin/env python
# conda activate /Users/norakearns/miniconda3/envs/motif
from fileinput import filename
from itertools import product
import re
from turtle import width
import cairo
import math
from IPython import display
from IPython.display import SVG, display, Image
import re
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="A program to input coverage limit")
    parser.add_argument("-f", help="fasta file", required=True, type=str)
    parser.add_argument("-m", help="motifs file", required=True, type=str)
    return parser.parse_args()

args = get_args()
fasta_filename = args.f 
motif_filename = args.m


IPUAC_dict = {'A':['A'],'C':['C'], 'G':['G'],'T':['T'],'U':['U'],
                'W':['A','T'],'S':['G','C'],'M':['A','C'],'K':['G','T'],
                'R':['A','G'],'Y':['C','T'], 'B':['C','G','T'],
                'D':['A','G','T'],'H':['A','C','T'],'V':['A','C','G'],
                'N':['A','C','T','G']}

def get_names(file):
    '''
    Takes a Fasta file and returns an array of just the sequence IDs
    Input: fasta_file_handle
    Output: list of sequence IDs
    '''
    file_handle = open(file, "rt")
    names_list = []
    lines = file_handle.readlines()
    for line in lines:
        if line.startswith(">"):
            names_list.append(line.strip("\n")[1:])
    return names_list


def get_fasta_records(file):
    '''
    Takes a Fasta file where the sequences are split across multiple lines, and concatenates the sequence lines to a single string. Returns the sequence line.
    Input: fasta_file_handle
    Output: sequence
    '''
    file_handle = open(file, "rt")
    record_list = []
    for line in file_handle:
        all_lines = file_handle.read()
        record_array = all_lines.split('>')
        record0_line_array = record_array[0].split('\n') 
        record0_seq_line = ('').join(record0_line_array)
        record_list.append(record0_seq_line)
        for record in record_array[1:]:
            record_line_array = record.split('\n')  
            seqs = record_line_array[1:]
            seq_line = ('').join(seqs)
            record_list.append(seq_line)
    return record_list

class Sequence:
    def __init__(self, sequence):
        self.sequence = sequence
        self.exon, self.seq_length = self.find_exon_introns(self.sequence)
    
    def find_exon_introns(self, sequence):
        '''
        This function takes as its input a DNA sequence with introns lowercase and exons UPPERCASE. It identifies the positions of the first intron, the exon
        and the second intron.
        Input: tttgtttccacagAAAAACCTCTTCAGGCACTGGTGCCGAGGACCCTAGgtatgactcacctgt"
        Output: ([0, 12], [13, 48], [49, 63])
        '''
        intron = re.compile("[actgn]+")
        seq_len = len(sequence)
        exon = re.compile("[ACTGN]+")
        exon_positions =  [[m.start(),m.end()] for m in exon.finditer(sequence)]
        intron_positions = [[m.start(),m.end()] for m in intron.finditer(sequence)]
        intron1_start_end = [intron_positions[0][0], (intron_positions[0][1]-1)] 
        exon_start_end = [exon_positions[0][0], (exon_positions[0][1]-1)]
        intron2_start_end = [intron_positions[1][0], (intron_positions[1][1]-1)] 
        return exon_start_end, seq_len    

    
    def get_motif_occurence_dict(self, motif_dict):
        '''
        This function takes as its input a fasta sequence and the dictionary of all possible motifs, and returns a dictionary of each motif and where it is located in the sequence.
        Input: {"YGCY":["TGCT","CGCT",...]...}, and "tttgtttccacagAAAAACCTCTTCAGGCACTGGTGCCGAGGACCCTAGgtatgactcacctgt"
        Output: {{"YGCY":[[25,29], [31,35]}}
        '''
        motif_occurence_dict = {}
        for motif in motif_dict.keys(): # for each motif
            motif_occurence_dict[motif] = []
            for possible_motif in motif_dict[motif]: # for each motif option (because motifs are ambiguous)
                motif_pattern = re.compile(possible_motif)
                #motif_patter_intron = re.compile(possible_motif.lower())
                motif_pattern_occurence = [[m.start(),m.end()] for m in motif_pattern.finditer((self.sequence).upper())]
                if len(motif_pattern_occurence) != 0:
                    motif_occurence_dict[motif].append(motif_pattern_occurence)
        return motif_occurence_dict


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

            possible_strings_iters = list(product(*all_characters)) # create all possible strings
            possible_strings = []
            for i in range(0, len((possible_strings_iters))):
                possible_strings.append(''.join(possible_strings_iters[i]))   
            motif_dict[motif] = possible_strings
        return motif_dict

    def get_all_motif_options(self):
        return self.motif_options_dict 



class Plot:
    def __init__(self, plot_width, plot_height, file_name):
        self.my_motif_colors = {}
        self.width = plot_width
        self.height = plot_height 
        self.file = file_name
        self.surface = cairo.SVGSurface(self.file, self.width, self.height) # create png with w/h dimensions
        self.context = cairo.Context(self.surface)

    def set_motif_colors(self, motif_string, colors): # setter
        self.my_motif_colors[motif_string] = colors

    def get_colors_for_motif(self, motif_string): # getter
        color_r = self.my_motif_colors[motif_string][0]
        color_g = self.my_motif_colors[motif_string][1]
        color_b = self.my_motif_colors[motif_string][2]
        return color_r, color_g, color_b
    
    def plot_motifs(self, motif_occurence_dict, sequence_pointer, seq, plot_index):
        exon = sequence_pointer.find_exon_introns(seq)[0]
        exon_start = exon[0]
        exon_len = (exon[1]-exon[0])
        seq_len = sequence_pointer.find_exon_introns(seq)[1]
        self.context.set_source_rgb(0,0,0)
        self.context.move_to(50,(100 + (200*(plot_index-1))))
        self.context.set_font_size(20)
        self.context.show_text(names_list[plot_index-1])
        self.context.set_line_width(3)
        self.context.set_source_rgb(0,0,0)
        self.context.move_to(100, (200*plot_index)) # origin of the line, left-most point
        self.context.line_to((100+seq_len),(200*plot_index)) # right-most point
        self.context.stroke() 
        self.context.rectangle((100+exon_start),(150 + (200*(plot_index-1))),exon_len,(100)) # (x, y, width, height)
        self.context.set_source_rgba(0,0,0, 0.1)
        self.context.fill()


        for motif in motif_occurence_dict.keys():
        # YCGY, GCAUG, CATAG, YYYYYYYYYY
            motif_blocks = [] # position array like this: [[79, 89], [80, 90], [76, 86], [78, 88], [75, 85], [77, 87], [74, 84], [73, 83]]
            for location_array in motif_occurence_dict[motif]:
                for location in location_array:
                    start = location[0]
                    length = (location[1]-location[0])
                    motif_blocks.append(location)
                    color_r, color_g, color_b = self.get_colors_for_motif(motif)
                    self.context.set_source_rgba(color_r, color_g, color_b, 0.5)
                    self.context.rectangle((100+start),(150 + (200*(plot_index-1))),length,(100)) # (x, y, width, height)
                    self.context.fill()
                    self.context.stroke()
    
        # Create Key
        self.context.rectangle(850, 100, 20, 20) # (x, y, width, height)
        self.context.set_source_rgba(0,0.2,0.8, 0.5) # blue
        self.context.fill()

        self.context.move_to(880,120)
        self.context.set_font_size(20)
        self.context.set_source_rgb(0,0,0)
        self.context.show_text("CATAG")

        self.context.rectangle(850, 140, 20, 20) # (x, y, width, height)
        self.context.set_source_rgba(0.2,0.9,0.8, 0.5) # green
        self.context.fill()

        self.context.move_to(880,160)
        self.context.set_font_size(20)
        self.context.set_source_rgb(0,0,0)
        self.context.show_text("YYYYYYYYYY")

        self.context.rectangle(850, 180, 20, 20) # (x, y, width, height)
        self.context.set_source_rgba(1,0.6,0.3, 0.5) # orange
        self.context.fill()

        self.context.move_to(880,200)
        self.context.set_font_size(20)
        self.context.set_source_rgb(0,0,0)
        self.context.show_text("YGCY")

        self.context.rectangle(850, 220, 20, 20) # (x, y, width, height)
        self.context.set_source_rgba(0,0,0,0.1) # orange
        self.context.fill()

        self.context.move_to(880,238)
        self.context.set_font_size(20)
        self.context.set_source_rgb(0,0,0)
        self.context.show_text("EXON")
    
        self.surface.write_to_png("Figure_1.png")


motifs_pointer = Motifs(motif_filename)
motif_dict = motifs_pointer.get_all_motif_options()

record_list = get_fasta_records(fasta_filename)
names_list = get_names(fasta_filename)

sequence1 = Sequence(record_list[0])
sequence2 = Sequence(record_list[1])
sequence3 = Sequence(record_list[2])
sequence4 = Sequence(record_list[3])

MOdict1 =  sequence1.get_motif_occurence_dict(motif_dict)
MOdict2 = sequence2.get_motif_occurence_dict(motif_dict)
MOdict3 = sequence3.get_motif_occurence_dict(motif_dict)
MOdict4 =  sequence4.get_motif_occurence_dict(motif_dict)

myplot = Plot(1200, 1200, "Figure_1.png")
myplot.set_motif_colors("GCUAG", [1,0.8,0.2 ]) # yellow
myplot.set_motif_colors("CATAG", [0,0.2,0.8 ]) # blue
myplot.set_motif_colors("YYYYYYYYYY", [0.2,0.9,0.8]) # green
myplot.set_motif_colors("YGCY", [1,0.6,0.3]) # orange

myplot.plot_motifs(MOdict1, sequence1, record_list[0], 1)
myplot.plot_motifs(MOdict2, sequence2, record_list[1], 2)
myplot.plot_motifs(MOdict3, sequence3, record_list[2], 3)
myplot.plot_motifs(MOdict4, sequence4, record_list[3], 4)


