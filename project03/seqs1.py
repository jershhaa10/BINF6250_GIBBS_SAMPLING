'''
for testing with just chr1 derived forward sequences
'''

import bamnostic 
bam = bamnostic.AlignmentFile("data/SRR9090854.subsampled_5pct.bam", "rb")
seqs1 = []
for read in bam.fetch('1', 0, 1000000):
    if getattr(read, "is_unmapped", True) or getattr(read, "is_secondary", True) or\
    getattr(read, "is_supplementary", True) or getattr(read, "is_reverse", True): 
        continue
    seq = getattr(read, "query_sequence", None)
    seqs1.append(seq)
#print(len(seqs1))
#now in other scripts can import seqs1 from seqs1 for testing
bam.close()