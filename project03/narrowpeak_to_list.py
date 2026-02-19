#!/usr/bin/env python3
"""
Read a MACS3 *_peaks.narrowPeak file and a reference FASTA, then writes two newline-delimited sequence files suitable for loading into Python lists:
    - window_seqs.txt: summit-centered windows (len=WINLEN, default 101). Intended as the motif-search input (`seqs`) to the Gibbs sampler.
    - flank_seqs.txt: fixed-length background per peak (len=2*FLANK, default 100), formed by concatenating the FLANK bases immediately left and right of the window. Intended as the background-training input (`bg_sample`) for a 0th-order background model, reducing noise contamination from the motif window itself.
"""

import argparse
import pandas as pd # to read narrowpeak as table
import pysam

FLANK = 50
WINLEN = 101

def get_window_and_flanks(
    fa: pysam.FastaFile, chrom: str, summit0: int,
    winlen: int = 101, flank: int = 50
    ) -> tuple[str, str] | None:

    """
    fetch one extended region around summit and cut it into:
    - window: centered on summit, length winlen
    - flanks: left flank + right flank outside window, total length 2*flank

    coordinates are 0-based half-open.

    skips if near contig ends.
    """

    if winlen % 2 == 0:
        raise ValueError("winlen must be odd to be summit-centered cleanly")
    win_half = winlen // 2

    radius = win_half + flank

    start0 = summit0 - radius
    end0 = summit0 + radius + 1 # half-open; pysam and BED-style format 

    # skip if near contig ends
    if start0 < 0:
        return None
    try:
        chrom_len = fa.get_reference_length(chrom)
    except Exception:
        return None
    if end0 > chrom_len:
        return None

    # defenses

    try:
        ext = fa.fetch(chrom, start0, end0).upper()
    except Exception:
        return None

    expected_len = 2 * radius + 1

    if len(ext) != expected_len:
        return None

    # layout is [left_flank][window][right_flank]
    # lengths are flank, winlen, flank

    left = ext[:flank]
    window = ext[flank:flank + winlen]
    right = ext[flank + winlen:]

    if len(left) != flank or len(window) != winlen or len(right) != flank:
        return None

    flanks = left + right

    return window, flanks

def main():
    # parse args
    ap = argparse.ArgumentParser()

    ap.add_argument("--narrowpeak", default="macs3_out/raw_tlen30-1000_peaks.narrowPeak", \
        help="MACS *_peaks.narrowPeak file path; default is output of sort_and_callpeaks.sh")

    ap.add_argument("--ref_fasta", default="ref/grch37_human_g1k_v37/human_g1k_v37.fasta", \
        help="reference fasta file path; default is output of get_grch37_human_g1k_v37.sh")

    ap.add_argument("--top", type=int, default=500, help="top N peaks by qValue (descending)")

    ap.add_argument("--window_out", default="window_seqs.txt", help="output path for motif windows")

    ap.add_argument("--flank_out", default="flank_seqs.txt", help="output path for background flank seqs")

    args = ap.parse_args()


    # reads macs *.narrowPeak into a DataFrame
    peaks = pd.read_csv(
        args.narrowpeak, 
        sep="\t", 
        header=None,
        names=["chrom","start","end","name","score","strand","signalValue","pValue","qValue","peak"]
    )

    # take only peaks with a defined summit
    peaks = peaks[peaks["peak"] != -1].copy()

    # select N most confident peaks by qValue
    peaks = peaks.sort_values("qValue", ascending=False).head(args.top)

    # prepare ref
    fa = pysam.FastaFile(args.ref_fasta)
    refs = set(fa.references) # lists contigs

    # loops over peaks, gets motif windows and flanks
    window_seqs = []
    flank_seqs = []

    for row in peaks.itertuples(index=False):
        chrom = str(row.chrom)
        if chrom not in refs:
            #chrom not in reference contigs
            continue

        summit0 = int(row.start) + int(row.peak)

        wf = get_window_and_flanks(
            fa, 
            chrom, 
            summit0,
            winlen=WINLEN, 
            flank=FLANK
            )

        if wf is None:
            continue

        w, fl = wf
        window_seqs.append(w)
        flank_seqs.append(fl)

    print(f"top peaks requested: {args.top}; after peak!=-1 filter: {len(peaks)}")
    print(f"kept windows: {len(window_seqs)} (len={WINLEN})")
    print(f"kept flanks (concatenated left+right): {len(flank_seqs)} (len={2*FLANK})")

    with open(args.window_out, "w") as f:
        for s in window_seqs:
            f.write(s + "\n")

    with open(args.flank_out, "w") as f:
        for s in flank_seqs:
            f.write(s + "\n")

if __name__ == "__main__":
    main()