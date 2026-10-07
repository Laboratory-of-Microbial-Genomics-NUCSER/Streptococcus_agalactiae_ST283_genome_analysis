#!/usr/bin/env python3
import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# Argument Parser

parser = argparse.ArgumentParser(
    description="Remove recombination regions from a core genome alignment using Gubbins GFF output."
)

parser.add_argument(
    "-i",
    "--input",
    nargs=2,
    required=True,
    metavar=("GFF_FILE", "ALIGNMENT_FILE"),
    help="Gubbins recombination_predictions.gff file AND core genome alignment (FASTA/ALN)"
)

parser.add_argument(
    "-o",
    "--output",
    required=True,
    help="Output FASTA file with recombination removed"
)

args = parser.parse_args()

gff_file = args.input[0]
alignment_file = args.input[1]
output_file = args.output

# STEP 1: Read recombination coordinates from GFF file
recomb_pos = set()

with open(gff_file) as gff:
    for line in gff:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 5:
            continue

        start = int(parts[3]) - 1  # convert 1-based → 0-based
        end = int(parts[4])
        recomb_pos.update(range(start, end))


# STEP 2: Read alignment
alignment = AlignIO.read(alignment_file, "fasta")
alignment_length = alignment.get_alignment_length()

keep_positions = [i for i in range(alignment_length) if i not in recomb_pos]

# STEP 3: Build recombination-free alignment
new_alignment = MultipleSeqAlignment([])

for record in alignment:
    new_seq = ''.join(record.seq[i] for i in keep_positions)
    new_record = SeqRecord(Seq(new_seq), id=record.id, description="")
    new_alignment.append(new_record)

# STEP 4: Save output
AlignIO.write(new_alignment, output_file, "fasta")
print(f"\nRecombination-removed alignment saved to: {output_file}\n")

