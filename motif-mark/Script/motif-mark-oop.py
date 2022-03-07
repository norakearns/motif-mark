#!/usr/bin/env python
import itertools
from itertools import product
import re
import cairo
import argparse
from IPython import display

def get_args():
    parser = argparse.ArgumentParser(description="A program to input coverage limit")
    parser.add_argument("-f", help="fasta file", required=True, type=str)
    parser.add_argument("-m", help="motifs file", required=True, type=str)
    return parser.parse_args()

# args = get_args()
fasta_filename = "Figure_1.fasta" # -f file carrying the sequences in fasta format
motif_filename = "Fig_1_motifs.txt" # -m file carrying the motifs
out_filename = fasta_filename.split('.')[0] + ".png"

# Dictionary which holds all the possible nucleotides represented by a specific IUPAC characters
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
        self.sequence = sequence # to create the sequence object you pass in one of the fasta sequences
        self.exon, self.seq_length = self.find_exon_introns(self.sequence)
    
    def find_exon_introns(self, sequence):
        '''
        This function takes as its input a DNA sequence with introns lowercase and exons UPPERCASE. It identifies the positions of the exon and the sequence length
        Input: tttgtttccacagAAAAACCTCTTCAGGCACTGGTGCCGAGGACCCTAGgtatgactcacctgt"
        Output: [13, 48], 100
        '''
        seq_len = len(sequence)
        exon = re.compile("[ACTGN]+") # exons are represented by uppercase letters
        exon_positions =  [[m.start(),m.end()] for m in exon.finditer(sequence)] # find the start and stop of the exon
        exon_start_end = [exon_positions[0][0], (exon_positions[0][1]-1)]
        return exon_start_end, seq_len    

    
    def get_motif_occurence_dict(self, motif_dict):
        '''
        This function takes as its input a fasta sequence and the dictionary of all possible motifs, and returns a dictionary of each motif and where it is located in the sequence.
        Input: {"YGCY":["TGCT","CGCT",...]...}, and "tttgtttccacagAAAAACCTCTTCAGGCACTGGTGCCGAGGACCCTAGgtatgactcacctgt"
        Output: {"YGCY":[[25,29], [31,35]...}
        '''
        motif_occurence_dict = {}
        for motif in motif_dict.keys(): # for each motif in the dictionary created in the motif class (which has as keys all the motifs, and values all the possible motifs an ambiguous motif could represent)
            motif_occurence_dict[motif] = [] # store the motif as a key in motif_occurence_dict, and set the value to an empty list
            for possible_motif in motif_dict[motif]: # for each motif option in motif_dict (because motifs are ambiguous)
                motif_pattern = re.compile(possible_motif) # create a regular expression pattern
                motif_pattern_occurence = [[m.start(),m.end()] for m in motif_pattern.finditer((self.sequence).upper())] # search for occurrences of it in the gene sequence. NOTE: the motifs have been made uppercase, so sequence.upper is used
                if len(motif_pattern_occurence) != 0: # if there is an occurence of the motif inside the gene
                    motif_occurence_dict[motif].append(motif_pattern_occurence) # add it to the value in the dictionary
        return motif_occurence_dict 


class Motifs:
    def __init__(self, motif_file_name):
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
        Output: {motif1: [[motif1.1, motif 1.2, motif 1.3], motif2: [motif2.1, motif 2.2, motif 2.3]}
        '''
        motif_dict = {}
        for motif in self.motif_array:
            motif = motif.upper() # make all letters uppercase
            all_characters = []
            for character in motif:
                one_position_possible_chars = IPUAC_dict[character] # create an array for each position in the motif that holds every possible character that could be at that position according to its IUPAC
                all_characters.append(one_position_possible_chars) 
            possible_strings_iters = list(product(*all_characters)) # itertools product function creates all possible combinations for each ambiguous motif, but it returns them as lists (C, C, T, G, ...)
            possible_strings = []
            for i in range(0, len((possible_strings_iters))):
                possible_strings.append(''.join(possible_strings_iters[i])) # use join to turn the output of itertools product into strings  
            motif_dict[motif] = possible_strings 
        return motif_dict

    def get_all_motif_options(self):
        return self.motif_options_dict 

class Plot:
    def __init__(self, plot_width, plot_height, file_name):
        self.my_motif_colors = {} # dictionary which holds motifs and their rgb values
        self.width = plot_width 
        self.height = plot_height 
        self.file = file_name # file name for the png file
        self.surface = cairo.SVGSurface(self.file, self.width, self.height) # create png with w/h dimensions
        self.context = cairo.Context(self.surface) # create context on which to draw
        
    def set_motif_colors(self, motif_string, colors): 
        '''
        This function takes as its input a motif string and a set of rgb colors, and creates a dictionary with the motif string as key and rgb numbers as values
        Input: "GCUAG", [1,0.8,0.2 ]
        Output: {"GCUAG": [1,0.8,0.2 ], ...}
        '''
        self.my_motif_colors[motif_string] = colors

    def get_colors_for_motif(self, motif_string): 
        '''
        This function takes as input a motif string and gets the rgb values for that motif string from the my_motif_colors dictionary and store them as color_r, color_b, color_g
        Input: a motif
        Outpu: 0.1, 0.6, 0.8 (some RGB colors)
        '''
        color_r = self.my_motif_colors[motif_string][0]
        color_g = self.my_motif_colors[motif_string][1]
        color_b = self.my_motif_colors[motif_string][2]
        return color_r, color_g, color_b
    
    def create_background_color(self): # set the backgroun color to white
        '''
        This function creates a white rectangle that is the height and width of the surface because cairo's default context background is black/transparent
        '''
        self.context.rectangle(0,0,self.width,self.height)
        self.context.set_source_rgb(1, 1, 1)
        self.context.fill()
    
    def plot_motifs(self, motif_occurence_dict, sequence_obj, seq, plot_index):
        '''
        This function creates the motif plots 
        Input: motif_occurence_dict (a dictionary of each motif and where it occurs in a sequence), a sequence object, a gene sequence, the index of the gen sequence in the gene sequence array)
        Output: A plot which creates color-coded blocks of the motifs on a sequence, which is represented as a black line
        '''
        exon = sequence_obj.find_exon_introns(seq)[0] # get the exon start and stop from the sequence object
        exon_start = exon[0] # the start is the first item in the exon list
        exon_len = (exon[1]-exon[0]) # the length of the exon is the second item - first item in the exon list
        seq_len = sequence_obj.find_exon_introns(seq)[1] # get the sequence length from the sequence object
        self.context.set_source_rgb(0,0,0) # set the color of the text labels to be black
        self.context.move_to(50,(100 + (200*(plot_index-1)))) # start the text at 50,100, and increment this start by 200 for every subsequent sequence
        self.context.set_font_size(20) 
        self.context.show_text(names_list[plot_index-1])
        self.context.set_line_width(3) # width of the line representing the sequence
        self.context.set_source_rgb(0,0,0) # color of sequence line = black
        self.context.move_to(100, (200*plot_index)) # origin of the line, left-most point
        self.context.line_to((100+seq_len),(200*plot_index)) # right-most point is sequence length - start position
        self.context.stroke() 
        self.context.rectangle((100+exon_start),(150 + (200*(plot_index-1))),exon_len,(100)) # Create the exon box (x, y, width, height)
        self.context.set_source_rgba(0,0,0, 0.1) # make it transparent
        self.context.fill()
        
        

        for motif in motif_occurence_dict.keys():
        # for each motif : YCGY, GCAUG, CATAG, YYYYYYYYYY
            motif_blocks = [] # position array like this: [[79, 89], [80, 90], [76, 86], [78, 88], [75, 85], [77, 87], [74, 84], [73, 83]]
            for location_array in motif_occurence_dict[motif]: # for the locations of each motif
                for location in location_array: 
                    start = location[0] # get the start of the motif
                    length = (location[1]-location[0]) # get the length of the motif
                    motif_blocks.append(location)
                    color_r, color_g, color_b = self.get_colors_for_motif(motif) # get the color assigned to that motif 
                    self.context.set_source_rgba(color_r, color_g, color_b, 0.5) 
                    self.context.rectangle((100+start),(150 + (200*(plot_index-1))),length,(100)) # create a rectangle (x, y, width, height) of the right size and color
                    self.context.fill()
                    self.context.stroke()

    def exon_key(self):
        self.context.rectangle(850, 60, 20, 20) # (x, y, width, height)
        self.context.set_source_rgba(0,0,0, 0.1) # blue
        self.context.fill()
        
        self.context.move_to(880,80)
        self.context.set_font_size(20)
        self.context.set_source_rgb(0,0,0) # black
        self.context.show_text("EXON")
            
    def create_key(self, motif, n):
        self.context.rectangle(850, (100 + 40*n), 20, 20) # (x, y, width, height)
        color_r, color_g, color_b = self.get_colors_for_motif(motif)
        self.context.set_source_rgba(color_r, color_g, color_b, 0.5) # blue
        self.context.fill()
        
        self.context.move_to(880,(120 + 40*n))
        self.context.set_font_size(20)
        self.context.set_source_rgb(0,0,0) # black
        self.context.show_text(motif)

        self.surface.write_to_png(out_filename)


motifs_obj = Motifs(motif_filename) # create motif object
motif_list = motifs_obj.get_motif_array() # ["YGCY","GCAUG","CATAG","YYYYYYYYYY"]
motif_dict = motifs_obj.get_all_motif_options()

record_list = get_fasta_records(fasta_filename)
names_list = get_names(fasta_filename)

myplot = Plot(1200, 1200, out_filename)
myplot.create_background_color()

# Dictionary of color lists. Length of the list is dependent on length of the motif list.
color_list_dict = { 4: [[1,0.6,0.3],[1,0.8,0.2 ],[0,0.2,0.8 ], [0.2,0.9,0.8]], 5: [[1,0.6,0.3],[1,0.8,0.2 ],[0,0.2,0.8 ], [0.2,0.9,0.8],[0.5,0.3,1]],
                       6: [[1,0.6,0.3],[1,0.8,0.2 ],[0,0.2,0.8 ], [0.2,0.9,0.8],[0.5,0.3,1], [0.7,0,0.1]],
                       7: [[1,0.6,0.3],[1,0.8,0.2 ],[0,0.2,0.8 ], [0.2,0.9,0.8],[0.5,0.3,1], [0.7,0,0.1], [0.5,0.7,0.6]],
                       8: [[1,0.6,0.3],[1,0.8,0.2 ],[0,0.2,0.8 ], [0.2,0.9,0.8],[0.5,0.3,1], [0.7,0,0.1], [0.5,0.7,0.6], [0.1,0.7,0.6]],
                       9: [[1,0.6,0.3],[1,0.8,0.2 ],[0,0.2,0.8 ], [0.2,0.9,0.8],[0.5,0.3,1], [0.7,0,0.1], [0.5,0.7,0.6], [0.1,0.7,0.6], [0.5,0.7,0.5]],
                       10: [[1,0.6,0.3],[1,0.8,0.2 ],[0,0.2,0.8 ], [0.2,0.9,0.8],[0.5,0.3,1], [0.7,0,0.1], [0.5,0.7,0.6], [0.1,0.7,0.6], [0.5,0.7,0.5], [0.5,0.3,0.4]] }

len_motif_list = len(motif_list)
color_list = color_list_dict[len_motif_list] # get a list of colors the length of the list of motifs so that each motif has a color

for i in zip(motif_list, color_list):
    myplot.set_motif_colors(i[0].upper(), i[1]) 
    
    
n = 0
for i in record_list:
    n += 1
    sequence = Sequence(i) # create sequence object
    MOdict = sequence.get_motif_occurence_dict(motif_dict)
    myplot.plot_motifs(MOdict, sequence, i, n)

myplot.exon_key()
    
for i in range(len(motif_list)):
    n = i
    motif = motif_list[i].upper()
    myplot.create_key(motif, n)
    