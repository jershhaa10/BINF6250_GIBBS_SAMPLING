# Introduction
Description of the project

# Pseudocode

```
GibbsMotifSampler(...)
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
Description of the team's learning points

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
Aaronie Jersha Jenyfred: Chantera and Linh were both amazing to work with. We were able to co-ordinate well with each other and clarify doubts and debugging challenges. For me, as already mentioned, this was definitely one of the hardest projects to work with. 
## Other member
Other members' reflections on the project

# Generative AI Appendix
As per the syllabus
