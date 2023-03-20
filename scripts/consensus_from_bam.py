import sys
import pysam
import operator
from collections import defaultdict
from Bio import SeqIO

## Bam file 
bam_file_name = sys.argv[1]
samfile = pysam.AlignmentFile(bam_file_name, "rb")

## Fasta file
ref = [str(r.seq) for r in SeqIO.parse(sys.argv[2], 'fasta')][0]

sequence = defaultdict(lambda: "N")

for pileupcolumn in samfile.pileup(sys.argv[3]): ## Contig name
	bases = {"A":0,"C":0,"T":0,"G":0}
	for pileupread in pileupcolumn.pileups:
		if not pileupread.is_del and not pileupread.is_refskip:
			base = pileupread.alignment.query_sequence[pileupread.query_position]
			if base in bases:
				bases[base] += 1
	# print(bases)
	best_base = max(bases.items(), key=operator.itemgetter(1))[0]
	if bases[best_base] >= 1: # At least one base - can change to 2 or 3.
		sequence[pileupcolumn.reference_pos] = best_base

i = 0
final = ''
for b in ref:
	final += sequence[i]
	i += 1

print(">" + bam_file_name.split(".bam")[0] + "_consensus")
print(final)
