import pandas as pd
import requests
import time
import os
from chipseq_preprocessing import PeaksProcessing

def extract_motif_positions(peaks_df, flank_size = 50):
    motif_positions_df = pd.DataFrame()
    motif_positions_df['chrom'] = peaks_df['chrom']
    motif_positions_df['start_pos'] = peaks_df['summit_pos'] - flank_size
    motif_positions_df['end_pos'] = peaks_df['summit_pos'] + flank_size

    return motif_positions_df

def extract_background_positions(peaks_df, flank_size = 50):
    background_positions_df = pd.DataFrame()
    background_positions_df['chrom'] = peaks_df['chrom']
    background_positions_df['start_pos_lf'] = peaks_df['start_pos']
    background_positions_df['end_pos_lf'] = peaks_df['summit_pos'] - flank_size - 1
    background_positions_df['start_pos_rf'] = peaks_df['summit_pos'] + flank_size + 1
    background_positions_df['end_pos_rf'] = peaks_df['end_pos']

    return background_positions_df

def get_seq_b37(chrom, start, end):
    # Ensembl uses GRCh37 (equivalent to b37 which is Broad's Institute version)
    # Resource: https://grch37.rest.ensembl.org/documentation/info/sequence_region
    url = f"https://grch37.rest.ensembl.org/sequence/region/human/{chrom}:{start}..{end}:1"
    headers = {"Content-Type": "text/plain"}

    # Resource: https://requests.readthedocs.io/en/latest/user/quickstart/
    response = requests.get(url, headers=headers)
    if response.ok:
        return response.text.upper()
    return None

def extract_motif_sequences(motif_df):
    sequences = []
    print("\nStarting to extract motif sequences...")
    for counter, (idx, row) in enumerate(motif_df.iterrows()):
        chrom = row['chrom']
        seq_start = row['start_pos']
        seq_end = row['end_pos']
        seq = get_seq_b37(chrom, seq_start, seq_end)
        if seq:
            sequences.append(seq)

        time.sleep(0.1)
        if counter % 10 == 0:
            print(f"Extracted {counter + 1}/{len(motif_df)}")
    print("\nFinished extracting motif sequences")
    return sequences

def extract_background_sequences(background_df):
    print("\nStarting to extract background sequences...")
    background_sequences = []
    for counter, (idx, row) in enumerate(background_df.iterrows()):
        chrom = row['chrom']
        seq_start_lf = row['start_pos_lf']
        seq_end_lf = row['end_pos_lf']
        seq_start_rf = row['start_pos_rf']
        seq_end_rf = row['end_pos_rf']
        seq_lf = get_seq_b37(chrom, seq_start_lf, seq_end_lf)
        seq_rf = get_seq_b37(chrom, seq_end_rf, seq_end_rf)
        if seq_lf and seq_rf:
            seq = seq_lf + seq_rf
            background_sequences.append(seq)

        time.sleep(0.1)  # Rate limiting
        if counter % 10 == 0:
            print(f"Extracted {counter + 1}/{len(background_df)}")
    print("\nFinished extracting background sequences")
    return background_sequences

def extract_sequences():
    peaks_df = PeaksProcessing()

    motif_df = extract_motif_positions(peaks_df)
    background_df = extract_background_positions(peaks_df)

    motif_sequences = extract_motif_sequences(motif_df)
    background_sequences = extract_background_sequences(background_df)

    return motif_sequences, background_sequences

if __name__ == "__main__":
    motif_sequences, background_sequences = extract_sequences()
    motif_seq_df = pd.DataFrame({'sequence': motif_sequences})
    motif_seq_df.to_csv("motif_sequences.csv", index=False, header=False)
    print(f"\nMotif sequences saved to {os.path.abspath('motif_sequences.csv')}")

    background_seq_df = pd.DataFrame({'sequence': background_sequences})
    background_seq_df.to_csv("background_sequences.csv", index=False, header=False)
    print(f"Background sequences saved to {os.path.abspath('background_sequences.csv')}")