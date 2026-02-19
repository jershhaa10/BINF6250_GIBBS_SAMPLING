
# Introduction
The project implements the Gibbs Sampling algorithm to identify the shared DNA motif regions across multiple sequences. The algorithm repeatedly updates the motif positions using probability-based sampling to improve the motif model. Over several iterations, the motif becomes stable and represents a conserved pattern in the sequences.

# Assumptions
* We assumed that the main Gibbs Sampling function would receive a list of sequences.
* We assumed that the motif is present in all of the sequences

# Input
Motif Sequences: We used `macs3` to derive peaks from the 3.2 million Chip-Seq fragments. We filtered the peak data, and retrieved possible motif sequences by flanking the summit of each peak by a default of 50 bps on each side. As the Chip-seq data was annotated against b37, a variant of GRCH37, we utilized GRCh37 as the reference genome to get the sequence information.

# Helper Functions


# Pseudocode
```
GibbsMotifSampler function
    Validate all the input parameters
      If k <= 0
         Raise an error
      If sequence list is empty
         Raise an error
    To build the background model
      Compute the overall base frequencies from the sequences
      Apply pseudocounts

    Initialization of motif positions and matrix
      Set random seed generation with seed
      Pick a start position in each of the sequence
      Build initial start list 
      Build initial PFM counts using motif_pseudocount

    Print the statement "burning in..."
    Set the parameters
    prev_pwm as None
    n as number of sequences

    Create a loop for sweeps - from 1 to max_sweeps:
      Do n updates per sweep, where each update changes the motif position in the sequence
        Randomly select a sequence index i
        Update the sequence
          Remove its current motif from the PFM
          For every possible k-mer in the sequence i
              Score the k-mers by comparing the current motif PWM to background
          Convert the current score into probability
          Randomly start a new position with those probabilities in hand
          Update the start[i] with the new position
          Add the newly chosen k-mer into the PFM
      Check if the sweep is within burn in range
        Yes, then continue sweep
      If the burn in ends
        Print "done burning in.."
      If the current sweep is not the multiple of check_every_sweep after burn in
        Continue

      Check for convergence
      Build the PWM from the current model, converting counts to probabilities

      If the prev PWM is None
          Store the current PWM as prev PWM for comparison (initially)
          Continue to the next sweep

      Calculate the difference between current PWM and prev PWM
      Calculate the Frobenius norm (magnitude) of the difference matrix and store as numerator
      Calculate the Frobenius norm (magnitude) of the prev PWM and store as denonimator
      Calculate the relative change
              If denominator > 0, else infinity

      Update the current PWM as prev PWM
      Calculate the information content for the motif
      Print the no of sweeps completed, IC and the relative change

      If relative change < threshold value
        Break

    Return PFM
```

# Successes
We were able to implement the Gibbs sampling algorithm, and we felt more confident in our approach and results when we included peaks calling. We did work on this individually, but were able to consolidate our thoughts and approaches smoothly.

# Struggles
Some struggles we had with this project includes conceptually understanding the Gibbs sampling algorithm. We also struggled with processing the dataset due to its size, so we eneded up subsampling and never running the full file. This was alleviated when we used peak calling. We also struggled with integrating the functions into our workflow initially.

# Personal Reflections
## Group Leader
Aaronie Jersha Jenyfred:  Chantera and Linh were great to work with. We were able to discuss and plan out the implementation together since that was the hardest part. We were able to align and validate our logic and at the same time have different viewpoints over the same script, as to how to write it with a flow. 

## Other members
Chantera Lazard: My team were great and they challenged me to think of the algorithm mathematically. We were able to learn from each other especialy as this algorithm was challenging for us all. I also respected the ambition in my teammates to go above and beyond in understanding the algorithm and code.

Ngoc Linh Nguyen: 

# Generative AI Appendix
None
