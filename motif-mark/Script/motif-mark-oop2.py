#!/usr/bin/env python
import itertools
import re
import cairo
import argparse
import logging

logging.basicConfig(level=logging.INFO)

IUPAC_dict = {
    'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'], 'U': ['U'], 'W': ['A', 'T'], 'S': ['G', 'C'], 
    'M': ['A', 'C'], 'K': ['G', 'T'], 'R': ['A', 'G'], 'Y': ['C', 'T'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'C', 'T', 'G']
}

def get_args():
    parser = argparse.ArgumentParser(description="Visualize motifs on gene sequences")
    parser.add_argument("-f", required=True, help="FASTA file with sequences")
    parser.add_argument("-m", required=True, help="Motifs file")
    return parser.parse_args()

def parse_fasta(filepath):
    sequences = []
    names = []
    with open(filepath, 'r') as f:
        seq = ''
        for line in f:
            if line.startswith(">"):
                names.append(line[1:].strip())
                if seq:
                    sequences.append(seq)
                    seq = ''
                else:
                    seq += line.strip()
        if seq:
            sequences.append(seq)
    return names, sequences

def read_motifs(filepath):
    """Read motifs from a file, returning a list of uppercase motif strings."""
    with open(filepath, 'r') as f:
        return [line.strip().upper() for line in f if line.strip()]

def expand_motif(motif):
    """Expand a single ambiguous motif into all possible IUPAC-resolved variants."""
    try:
        chars = [IUPAC_dict[c] for c in motif]
    except KeyError as e:
        raise ValueError(f"Invalid IUPAC code: {e}")
    return [''.join(p) for p in itertools.product(*chars)]

class Sequence:
    def __init__(self, seq):
        self.seq = seq
        self.exon_coords = self._find_exon_coords()

    def _find_exon_coords(self):
        """Identify the start and end coordinates of exonic (uppercase) regions."""
        match = re.search(r'[A-Z]+', self.seq)
        if not match:
            return (0,0)
        return (match.start(), match.end())

    def motif_locations(self, motif_dict):
        """Return a dictionary of motif locations within the sequence."""
        locations = {}
        upper_seq = self.seq.upper()
        for motif, variants in motif_dict.items():
            locs = []
            for var in variants:
                locs.extend([[m.start(), m.end()] for m in re.finditer(var, upper_seq)])
            locations[motif] = locs
        return locations

class Plot:
    def __init__(self, width, height, outname):
        """Initialize a Cairo SVG surface for drawing motifs and gene sequences"""
    self.width = width
    self.height = height
    self.filename = outname
    self.surface = cairo. SVGSurface(outname, width, height)
    self.ctx = cairo.Context(self.surface)
    self.colors = {}
    self._draw_background()

    def _draw_background(self):
        """Fill the entire surface with a white background."""
        self.ctx.rectangle(0,0,self.width, self.height)
        self.ctx.set_source_rgb(1, 1, 1)
        self.ctx.fill()

    


    