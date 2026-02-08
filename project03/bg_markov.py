BASES = ("A", "C", "G", "T")

# 0th order MC (just iid?) 

def build_markov_0th_order(reads:list, pseudocount:float):
    base_counts = {b: 0 for b in BASES}

    for read in reads:
        if not read:
            continue
        for ch in read.upper():
            if ch in BASES:
                base_counts[ch] += 1

    total = sum(base_counts[b] for b in BASES) + pseudocount * 4
    return {b: (base_counts[b] + pseudocount) / total for b in BASES}

# 1st order MC

def build_bg_markov_1st_order(reads:list, pseudocount:float):
    init_counts = {b: 0 for b in BASES}
    trans_counts = {a: {b: 0 for b in BASES} for a in BASES} 

    # for each base, count how many times it's the "next" base, 
    # incl. when no previous base (it is the initial base in a read)
    for read in reads:
        read = read.upper()
        prev = None #previous valid base
        for ch in read:
            if ch not in BASES:
                prev = None  # this breaks the chain at N
                continue
            if prev is None:
                init_counts[ch] += 1
            else:
                trans_counts[prev][ch] += 1
            prev = ch

    # init probabilities P( x1 = base)
    init_total = sum(init_counts[b] for b in BASES) + pseudocount * 4 # total initial bases 
    init_probs = {b: (init_counts[b] + pseudocount) / init_total for b in BASES} # dict of probs of initial base

    # transition probabilities P(next = b | prev = a)
    trans_probs = {}
    for cur in BASES:
        row_total = sum(trans_counts[cur][nxt] for nxt in BASES) + pseudocount * 4 # all transitions from cur
        trans_probs[cur] = {nxt: (trans_counts[cur][nxt] + pseudocount) / row_total for nxt in BASES}

    return {"init": init_probs, "trans": trans_probs}

# we could just ignore init since it doesn't make sense.

if __name__ == "__main__":
    from seqs1 import seqs1
    reads = seqs1
    
    bg_1st = build_bg_markov_1st_order(reads, pseudocount=1)
    bg_1st_init = bg_1st["init"]
    bg_1st_trans = bg_1st["trans"]
    bg_0th = build_markov_0th_order(reads, pseudocount=1)

    #print(bg_0th )
    #print(bg_1st_init)
    #print(bg_1st_trans)