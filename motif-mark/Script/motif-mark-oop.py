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


IPUAC_dict = {'A':['A'],'C':['C'], 'G':['G'],'T':['T'],'U':['U'],
                'W':['A','T'],'S':['G','C'],'M':['A','C'],'K':['G','T'],
                'R':['A','G'],'Y':['C','T'], 'B':['C','G','T'],
                'D':['A','G','T'],'H':['A','C','T'],'V':['A','C','G'],
                'N':['A','C','T','G']}


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
        #self.intron1, self.exon, self.intron2 = self.find_exon_introns(self.sequence)
    
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
                motif_pattern_exon = re.compile(possible_motif)
                motif_patter_intron = re.compile(possible_motif.lower())
                motif_pattern_intron_occurence = [[m.start(),m.end()] for m in motif_patter_intron.finditer(self.sequence)]
                if len(motif_pattern_intron_occurence) != 0:
                    motif_occurence_dict[motif].append(motif_pattern_intron_occurence)
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
    def __init__(self):
        self.my_motif_colors = {}
    #     self.width = plot_width
    #     self.height = plot_height 
    #     surface = cairo.SVGSurface("line_and_rectangle.png", width, height)
    #     self.file = file_name
    #     self.surface, self.context = cairo.SVGSurface(self.file, self.width, self.height) # create png with w/h dimensions
    #     self.context = self.create_context(self.file, self.width, self.height) # create the surface on which to draw

    # def create_context(self, file, width, height):
    #     surface = cairo.SVGSurface(file, width, height)
    #     context = cairo.Context(surface)
    #     return context

    def set_motif_colors(self, motif_string, colors): # setter
        self.my_motif_colors[motif_string] = colors

    def get_colors_for_motif(self, motif_string): # getter
        color_r = self.my_motif_colors[motif_string][0]
        color_g = self.my_motif_colors[motif_string][1]
        color_b = self.my_motif_colors[motif_string][2]
        return color_r, color_g, color_b
    
    def plot_motifs(self, motif_occurence_dict, sequence_pointer, seq, plot_index, context):
        exon = sequence_pointer.find_exon_introns(seq)[0]
        exon_start = exon[0]
        exon_len = (exon[1]-exon[0])
        seq_len = sequence_pointer.find_exon_introns(seq)[1]
        context.set_line_width(3)
        context.set_source_rgb(0,0,0)
        context.move_to(100, (200*plot_index)) # origin of the line, left-most point
        context.line_to((100+seq_len),(200*plot_index)) # right-most point
        context.stroke() 
        context.rectangle((100+exon_start),(150 + (200*(plot_index-1))),exon_len,(100)) # (x, y, width, height)
        context.set_source_rgba(0,0,0, 0.1)
        context.fill()


        for motif in motif_occurence_dict.keys():
        # YCGY, GCAUG, CATAG, YYYYYYYYYY
            motif_blocks = [] # position array like this: [[79, 89], [80, 90], [76, 86], [78, 88], [75, 85], [77, 87], [74, 84], [73, 83]]
            for location_array in motif_occurence_dict[motif]:
                for location in location_array:
                    start = location[0]
                    length = (location[1]-location[0])
                    motif_blocks.append(location)
                    color_r, color_g, color_b = self.get_colors_for_motif(motif)
                    context.set_source_rgb(color_r, color_g, color_b)
                    context.rectangle((100+start),(150 + (200*(plot_index-1))),length,(100)) # (x, y, width, height)
                    context.fill()
                    #context.stroke()



motifs_pointer = Motifs("Fig_1_motifs.txt")
motif_dict = motifs_pointer.get_all_motif_options()

record_list = get_fasta_records("Figure_1.fasta")
sequence1 = Sequence(record_list[0])

sequence2 = Sequence(record_list[1])
sequence3 = Sequence(record_list[2])
sequence4 = Sequence(record_list[3])

MOdict1 =  sequence1.get_motif_occurence_dict(motif_dict)
MOdict2 = sequence2.get_motif_occurence_dict(motif_dict)
MOdict3 = sequence3.get_motif_occurence_dict(motif_dict)
MOdict4 =  sequence4.get_motif_occurence_dict(motif_dict)

height = 1200
width = 1200

surface = cairo.SVGSurface("line_and_rectangle.svg", width, height)
context = cairo.Context(surface)

myplot = Plot()
myplot.set_motif_colors("GCUAG", [255,255,0]) # yellow
myplot.set_motif_colors("CATAG", [0,0,255]) # blue
myplot.set_motif_colors("YYYYYYYYYY", [255,0,0]) # red
myplot.set_motif_colors("YGCY", [255,0,255]) # fuschia

myplot.plot_motifs(MOdict1, sequence1, record_list[0], 1, context)
myplot.plot_motifs(MOdict2, sequence2, record_list[1], 2, context)
myplot.plot_motifs(MOdict3, sequence3, record_list[2], 3, context)
myplot.plot_motifs(MOdict4, sequence4, record_list[3], 4, context)









#surface.write_to_png("line_and_rectangle.png")



# for motif in motif_occurence_dict.keys():
#     motif_blocks = []
#     for location_array in motif_occurence_dict[motif]:
#         for location in location_array:
#             start = location[0]
#             length = location[1]
#             motif_blocks.append(location)


#sequence = Sequence("Figure_1.fasta")

# what do I actually need to know?
# I need to know the exon start position and end position, and the sequence length.
 # if motif == "GCAUG":
                    #     context.set_source_rgb(255,255,0) # yellow
                    # elif motif == "CATAG":
                    #     context.set_source_rgb(0,0,255) # blue
                    # elif motif == "YYYYYYYYYY":
                    #     context.set_source_rgb(255,0,0) # reds
                    # else:
                    #     context.set_source_rgb(255,0,255) # fuschia
# create gorup class to group things together
# offset gets stored inside the group