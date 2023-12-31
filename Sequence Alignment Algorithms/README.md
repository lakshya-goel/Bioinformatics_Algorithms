
# Sequence Alignment

## Needleman-Wunsch algorithm
This algorithm performs global sequence alignment, where the entire sequences are aligned. It considers all possible alignments and calculates an optimal alignment score using dynamic programming. The algorithm guarantees finding the best alignment but can be computationally expensive for long sequences.

## Smith-Waterman algorithm
This algorithm is similar to Needleman-Wunsch but is used for local sequence alignment. It identifies the best alignment within a subset of the sequences rather than aligning the entire sequences. It also uses dynamic programming and is useful for identifying regions of similarity between sequences.

## Screenshots
![Needleman-wunsch](Images/needleman-wunsch.png)
![Smith-waterman](Images/smith-waterman.png)

