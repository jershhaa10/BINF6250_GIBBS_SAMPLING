import pandas as pd
import numpy as np

def load_peaks(narrowPeaks_file):
    """Load MACS3 narrowPeak file into df"""
    peaks = pd.read_csv(narrowPeaks_file, sep="\t", header=None,
                        names=["chrom", "start_pos", "end_pos", "name", "score", "strand", "signalValue", "pValue",
                               "qValue", "peak"])
    peaks['summit_pos'] = peaks['start_pos'] + peaks['peak']
    return peaks

def filter_peaks_by_peak_width(peaks_df, min_width = 100, max_width = 1000):
    """Filter peaks by removing abnormally narrow or wide peaks"""
    print(f"\nFiltering peaks by peak width where width is between {min_width} and {max_width}")
    print(f"Original Number of Peaks: {len(peaks_df)}")
    filtered_peaks = peaks_df[peaks_df["end_pos"] - peaks_df["start_pos"] >= min_width]
    filtered_peaks = filtered_peaks[filtered_peaks["end_pos"] - filtered_peaks["start_pos"] <= max_width]
    print(f"Number of Peaks after filtering: {len(filtered_peaks)}")
    return filtered_peaks

def filter_peaks_by_signalValue(peaks_df, threshold):
    """Filter peaks """
    print(f"\nFiltering peaks by signal value")
    print(f"Original Number of Peaks: {len(peaks_df)}")
    filtered_peaks = peaks_df[peaks_df["signalValue"] >= threshold]
    print(f"Number of Peaks after filtering: {len(filtered_peaks)}")
    return peaks_df[peaks_df["signalValue"] > threshold]

def PeaksProcessing(filter_by_width = True, min_width = 100, max_width = 1000, filter_by_signal = True, q = 0.75):
    print(f"\nStarting Processing Chip-seq Peaks data")
    # Create Df for peaks data
    peaks_path = "macs3_output/p53_peaks.narrowPeak"
    peaks = load_peaks(peaks_path)

    # Filter peaks by peaks width and by signalValue (summit height)
    if filter_by_width:
        peaks = filter_peaks_by_peak_width(peaks, min_width, max_width)
    if filter_by_signal:
        quantile_75 = np.quantile(peaks['signalValue'], q)
        quantile_75 = np.round(quantile_75, 0)
        peaks = filter_peaks_by_signalValue(peaks, threshold=quantile_75)

        print(f"\n{peaks.head()}\nChip-seq Peaks data has been processed.")

    return peaks




