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

    def assign_colors(self, motifs, palette):
        """Assign RGB colors to each motif using a predefined color palette."""
        for motif, color in zip(motifs, palette):
            self.colors[motif] = color

    def draw(self, names, sequences, motif_dict):
        """Draw gene sequences, exons, and motif locations for all input sequences."""
        for idx, (name, seq) in enumerate(zip(names, sequences)):
            y_offset = 200 * (idx + 1)
            seq_obj = Sequence(seq)
            self.ctx.set_source_rgb(0, 0, 0)
            self.ctx.move_to(50, y_offset - 100)
            self.ctx.set_font_size(20)
            self.ctx.show_text(name)

            self.ctx.set_line_width(3)
            self.ctx.move_to(100, y_offset)
            self.ctx.line_to(100 + len(seq), y_offset)
            self.ctx.stroke()

            exon_start, exon_end = seq_obj.exon_coords
            self.ctx.set_source_rgba(0, 0, 0, 0.1)
            self.ctx.rectangle(100 + exon_start, y_offset - 50, exon_end - exon_start, 100)
            self.ctx.fill()

            motif_locs = seq_obj.motif_locations(motif_dict)
            for motif, locs in motif_locs.items():
                color = self.colors[motif]
                self.ctx.set_source_rgba(*color, 0.5)
                for start, end in locs:
                    self.ctx.rectangle(100 + start, y_offset - 50, end - start, 100)
                    self.ctx.fill()

    def draw_legend(self, motifs):
        """Draw a color legend for motifs and exon regions on the SVG plot."""
        for i, motif in enumerate(motifs):
            y = 100 + 40 * i
            self.ctx.rectangle(850, y, 20, 20)
            self.ctx.set_source_rgba(*self.colors[motif], 0.5)
            self.ctx.fill()
            self.ctx.move_to(880, y + 15)
            self.ctx.set_font_size(20)
            self.ctx.set_source_rgb(0, 0, 0)
            self.ctx.show_text(motif)

        self.ctx.rectangle(850, 60, 20, 20)
        self.ctx.set_source_rgba(0, 0, 0, 0.1)
        self.ctx.fill()
        self.ctx.move_to(880, 75)
        self.ctx.show_text("EXON")

        self.surface.finish()

def main():
    """Main execution function to run the motif visualization pipeline."""
    args = get_args()
    names, sequences = parse_fasta(args.f)
    motifs = read_motifs(args.m)
    motif_dict = {m: expand_motif(m) for m in motifs}
    out_file = args.f.rsplit('.', 1)[0] + '.svg'

    palette = [
        [1, 0.6, 0.3], [1, 0.8, 0.2], [0, 0.2, 0.8], [0.2, 0.9, 0.8], [0.5, 0.3, 1],
        [0.7, 0, 0.1], [0.5, 0.7, 0.6], [0.1, 0.7, 0.6], [0.5, 0.7, 0.5], [0.5, 0.3, 0.4]
    ][:len(motifs)]

    plot = Plot(1200, 200 * len(sequences) + 100, out_file)
    plot.assign_colors(motifs, palette)
    plot.draw(names, sequences, motif_dict)
    plot.draw_legend(motifs)

if __name__ == "__main__":
    main()



    