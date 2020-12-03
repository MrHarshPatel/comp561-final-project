"""Usage: python3 {unaligned_fasta_file}
"""
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio.Align import AlignInfo
import sys

# Align using muscle.
muscle_exe = "../third_party/muscle.exe"
in_file = sys.argv[1]
out_aligned_file = sys.argv[1].split('.')[0] + '-aligned.fasta'
muscle_cline = MuscleCommandline(muscle_exe, input=in_file, out=out_aligned_file)
stdout, stderr = muscle_cline()

# Get consensus sequence and put in separate file.
alignment = AlignIO.read(out_aligned_file, 'fasta')
summary_align = AlignInfo.SummaryInfo(alignment)

out_consensus_file = sys.argv[1].split('.')[0] + '-consensus.txt'
with open(out_consensus_file, 'w+') as f:
    f.write(str(summary_align.dumb_consensus()))
