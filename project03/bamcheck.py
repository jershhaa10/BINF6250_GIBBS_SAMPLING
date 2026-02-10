import bamnostic 

# bam check
bam = bamnostic.AlignmentFile("data/SRR9090854.subsampled_5pct.bam", "rb")

print(type(bam.header))
print("header sections:", bam.header.keys())
sq = bam.header.get("SQ", [])
print("n_references:", bam.nreferences, ", len(sq):", len(sq))
print("first 3 SQ records:", sq[:3])
print("first 3 contig lengths:", bam.lengths[:3])

print("Reference contigs:", bam.references)


bam.close()