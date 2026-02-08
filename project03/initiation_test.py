import random
import numpy as np

BASES = ("A", "C", "G", "T")
BASE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}

def is_valid_kmer(seq): # boolean helper 
    return all(ch in BASE_INDEX for ch in seq)

def build_pfm(candidates, k):
    '''
    input: list of candidate kmers, length k 
    output: pfm matrix; how many of my currently selected k-mers have base b at position j?
    '''
    pfm = np.zeros((4, k), dtype=int) # initialize pfm array
    for c in candidates:
        for j, ch in enumerate(c): 
            # j increments position by 1, ch goes to next character of string
            pfm[BASE_INDEX[ch], j] += 1
    return pfm


def build_pwm(pfm, bg, pseudocount=0.25):
    '''
    input: pfm matrix, bg base probs (0th order), pseudocount 
    output: pwm matrix
    '''
    k = pfm.shape[1] # non-hardcoded k here
    pwm = np.zeros_like(pfm, dtype=float) # make zero array with same shape as pfm

    # get count totals for columns
    # smoothing with pseudocount to guard against 0-occurence
    # should be == number of reads if no Ns
    col_totals = pfm.sum(axis=0) + 4 * pseudocount

    for b in BASES: # this is the part we fix if we want to go 1st order
        bi = BASE_INDEX[b]
        bg_b = bg[b]  # for 0th order
    
        for j in range(k):
            # for each position in the motif,
            
            # turn counts to probabilities
            p_motif = (pfm[bi, j] + pseudocount) / col_totals[j]

            # convert to log-odds vs background
            pwm[bi, j] = np.log2(p_motif / bg_b)

    return pwm


def build_pwm_vectorized(pfm, bg, pseudocount) -> np.ndarray:
    """
    PWM[b,j] = log2( P_motif(b|j) / P_bg(b) )
    pfm: (4, k) counts
    bg:  dict like {"A":pA,"C":pC,"G":pG,"T":pT}
    """
    # this one's simpler but it's just for 0th order

    # bg is a (4, 1) column vector in A,C,G,T order
    bg_vec = np.array([bg[b] for b in BASES], dtype=float).reshape(4, 1)

    # column sums is (1, k)
    col_totals = pfm.sum(axis=0, keepdims=True) + 4 * pseudocount

    # motif probabilities is (4, k)
    p_motif = (pfm.astype(float) + pseudocount) / col_totals

    # log-odds conversion
    pwm = np.log2(p_motif) - np.log2(bg_vec)

    return pwm

def init_sites_and_matrices(reads, k, bg, seed, pseudocount=0.25):
    """
    initiation:
    - choose random start positions for each read
    - build initial PFM/PWM from chosen k-mers

    input: reads, k, bg, seed, motif pseudocount
    output: starts, candidates, pfm, pwm
    """
    rng = random.Random(seed)

    starts = []
    candidates = []

    max_tries = 50

    for read in reads:
        read=read.upper()
        L = len(read)
        if L < k:
            raise ValueError("All reads must have length >= k")
        
        for _ in range(max_tries):
            start = rng.randrange(0, L - k + 1)
            kmer = read[start:start + k]
            if is_valid_kmer(kmer):
                starts.append(start)
                candidates.append(kmer)
                break 
        else:
            raise ValueError("could not sample a valid initial k-mer from read")
        

    pfm = build_pfm(candidates, k)
    pwm = build_pwm(pfm, bg, pseudocount=pseudocount)

    return starts, candidates, pfm, pwm

if __name__ == "__main__":
    from seqs1 import seqs1 as reads
    from bg_markov import build_markov_0th_order

    k = 9
    seed = 123
    bg_pseudocount = 1
    motif_pseudocount = 0.25

    bg = build_markov_0th_order(reads, bg_pseudocount)

    starts, candidates, pfm, pwm = \
        init_sites_and_matrices(reads, k, bg, seed, motif_pseudocount)

    #print("start positions (first 5):", starts[:5])
    #print("initial candidates (first 5):", candidates[:5])
    #print("PFM:\n", pfm)
    #print("PWM:\n", pwm)