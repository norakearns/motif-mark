# OOCA_1

**Objective:** Develop a python script to plot protein binding motifs on an image of an exon and flanking introns

**Input files:** fasta files, motif file

**Output files:** one svg file for each fasta entry that visualizes introns as a line and exons as a large rectangle on that line, and small different colored rectangles as different motifs. 

1) Read in fasta file

def Read_fasta(fasta_file_handle):
'''
A function which takes as its input a fasta and creates a tuple for the file:
[name, sequence]

Input: fasta_file_handle
Output: [name, sequence] arrays
'''

2) Read in motif file

def Read_motif(motif_file_handle):
'''
A function which takes as its input a file of motifs (10 bases each, one motif per line) and outputs a list of motifs

Input: motif_file_handle
Output: array of motifs [motif1, motif2, motif3...]
'''

3) Generate all possible motifs

def Create_motifs(array of motifs):
'''
A function which takes as its input the initial array of motifs and then identifies ambiguous motifs and creates all possible motifs from those ambiguous motifs;
Input: [motif1, motif2, motif3]
Output: [[motif1.1, motif 1.2, motif 1.3],[motif2.1, motif 2.2, motif 2.3], ...]
'''

**class Sequence:**

    arguments: input string from Fasta, list of all possible motifs 
    
    attributes: 
    - length
    - intron 1 positions: start/end position of lowercase letters: exon 1
    - exons positions: start/end position of UPPERCASE letters
    - intron 2 positions: start/end position of lowercase letters: exon 2
    - motif locations

    methods:
    - find length
    - use regex to find positions of introns and exons
        - return intron 1 positions
        - return exons positions
        - return intron 2 positions
    - find motifs
    '''
    Input: takes as input an array of motifs (including all their possibilities)
    marches across string and each time it finds a motif it adds the positions to a dictionary
    Output: dictionary of motifs and position {motif1: [[start, stop], [start, stop]}, ...}
    '''

**class Plot:**

    argument: a Sequence object

    attributes: 
    - labels

    methods:
    - get locations of motifs from sequence
    - plot a line representing the sequence
    - plot a rectangle representing the exon
    - plot different colored rectangles for each motif


**Working code to generate a line and a rectangle, not at the origin, using pycairo**
'''
import cairo
from IPython import display

width, height = 1000, 500

#create the coordinates to display your graph
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
'''

**The image that your code created**

![line_rectangle_png](/Users/norakearns/bioinformatics/Bi625/motif-mark/line_and_rectangle.png)