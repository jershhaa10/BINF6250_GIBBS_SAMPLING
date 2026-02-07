import numpy as np
from seq_ops import reverse_complement
from typing import List, Tuple, Optional


TRANS_TAB = str.maketrans('ACGT', '0123')


def build_pfm(sequences: List[str], length: int) -> np.ndarray:
    """
    Build a Position Frequency Matrix (PFM) from a list of sequences.

    Args:
        sequences (List[str]): List of sequence strings.
        length (int): Size of the PFM to build.

    Returns:
        np.ndarray: PFM with dimensions 4 x length.
    """
    # Returns a 2D array of given shape (4 x length) filled with zeros
    pfm = np.zeros((4, length), dtype=np.int32)

    # Create a 1D array of length 256 filled with zeros
    # Goal: Create a faster lookup table
    base_to_index = np.zeros(256, dtype=np.int8)
    # `ord` takes a character and returns its unicode code point A = 65
    base_to_index[ord('A')] = 0 # ord('A') = 65, so base_to_index[65] = 0
    base_to_index[ord('C')] = 1 # ord('C') = 67, so base_to_index[65] = 1
    base_to_index[ord('G')] = 2 # ord('G') = 71, so base_to_index[71] = 2
    base_to_index[ord('T')] = 3 # ord('T') = 84, so base_to_index[84] = 3

    # Skipping k-mers containing 'N'
    # It would allow for a cleaner model but potential data loss... In Future Dev, could impute
    # Addition: N is also in sequences will need to skip those
    base_to_index[ord('N')] = -1

    # for each seq
    for seq in sequences:
        # The encode() method encodes the string, using the specified encoding. If no encoding is specified, UTF-8 will be used.
        # The numpy.frombuffer() function is used to interpret a buffer as a 1-dimensional NumPy array.


        # Based on notes above my understanding is that it is getting k-mer, encoding it to bytes and interpreting those
        # bytes as integers dtype = np.int8
        # This pairs with the lookup table
        seq_array = np.frombuffer(seq[:length].encode(), dtype=np.int8)
        # get the indices for the nucleotide
        indices = base_to_index[seq_array]
        # If N is in a sequence, skip it entirely
        if -1 in indices:
            continue
        else:
            # Performs in-place addition on specific elements in array
            # pfm is the target array we want to modify
            # Example: ACAGTGG
            # indices = [0,1,0,2,3,2,2]
            # np.arrange(length) = [0,1,2,3,4,5,6] if k = 7
            #so when it is added, by 1 it will be pfm[0,0] +=1, pfm[1,1] += 1,pfm[0,2] += 1
            np.add.at(pfm, (indices, np.arange(length)), 1)

    return pfm


def build_pwm(pfm: np.ndarray) -> np.ndarray:
    """
    Build a Position Weight Matrix (PWM) from a Position Frequency Matrix (PFM).

    Args:
        pfm (np.ndarray): PFM with dimensions 4 x length.

    Returns:
        np.ndarray: PWM with dimensions 4 x length.
    """
    p = 0.25
    bg = 0.25
    
    sums = np.sum(pfm, axis=0) + 4 * p
    pwm = np.log2((pfm + p) / sums[:, np.newaxis].T) - np.log2(bg)
    
    return pwm


def score_kmer(seq: str, pwm: np.ndarray) -> float:
    """
    Score a k-mer using a Position Weight Matrix (PWM).

    The k-mer length is expected to be the same as the PWM length.

    Args:
        seq (str): k-mer to score.
        pwm (np.ndarray): PWM for scoring.

    Returns:
        float: PWM score for the k-mer.

    Raises:
        ValueError: If the k-mer and PWM are different lengths.
    """
    if len(seq) != len(pwm[0]):
        raise ValueError('K-mer and PWM are different lengths!')

    indices = np.fromiter((int(nt) for nt in seq.translate(TRANS_TAB)), dtype = int)
    return np.sum(pwm[indices, np.arange(len(indices))])


def pfm_ic(pfm: np.ndarray) -> float:
    """
    Calculate the information content of a Position Frequency Matrix (PFM).

    Args:
        pfm (np.ndarray): Position Frequency Matrix with dimensions 4 x length.

    Returns:
        float: The information content of the PFM.
    """
    p = 0.25
    sums = np.sum(pfm, axis=0) + 4 * p
    fij = (pfm + p) / sums
    ic = np.sum(2 + np.sum(fij * np.log2(fij), axis=0, where=fij > 0))
    return ic
