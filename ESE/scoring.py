import sys
import pysam
import csv
import numpy as np
from Bio.Seq import Seq
from Bio import motifs
import cProfile

## Bring in arguments
# file = sys.argv[1]
# output = sys.argv[2]
# reference_genome = sys.argv[3]

## Define ESE PWM arrays (Source: ESEFinder)
# SF2/ASF
SRSF1 = [[-1.14,1.37,-0.21,-1.58],[0.62,-1.1,0.17,-0.5],
        [-1.58,0.73,0.48,-1.58],[1.32,0.33,-1.58,-1.13],
        [-1.58,0.94,0.33,-1.58],[-1.58,-1.58,0.99,-1.13],
        [0.62,-1.58,-0.11,0.27]]

# SF2/ASF (IgM-BRCA1)
SRSF1_igM = [[-1.58,1.55,-1.35,-1.55],[0.15,-0.53,0.44,-0.28],
            [-0.97,0.79,0.41,-1.28],[0.74,0.33,-0.98,-0.92],
            [-1.19,0.72,0.51,-1.09],[-0.75,-0.62,1.03,-0.52],
            [0.43,-0.99,0,0.2]]

# SC35
SRSF2 = [[-0.88,-1.16,0.87,-1.18],[0.09,-1.58,0.45,-0.2],
        [-0.06,0.95,-1.36,0.38],[-1.58,1.11,-1.58,0.88],
        [0.09,0.56,-0.33,-0.2],[-0.41,0.86,-0.05,-0.86],
        [-0.06,0.32,-1.36,0.96],[0.23,-1.58,0.68,-1.58]]

# SRp40
SRSF5 = [[-0.13,0.56,-1.58,0.92],[-1.58,0.68,-0.14,0.37],
        [1.28,-1.12,-1.33,0.23],[-0.33,1.24,-0.48,-1.14],
        [0.97,-0.77,-1.58,0.72],[-0.13,0.13,0.44,-1.58],
        [-1.58,-0.05,0.8,-1.58]]

# SRp55
SRSF6 = [[-0.66,0.39,-1.58,1.22],[0.11,-1.58,0.72,-1.58],
        [-0.66,1.48,-1.58,0.07],[0.11,-1.58,0.72,-1.58],
        [-1.58,-1.58,0.21,1.02],[0.61,0.98,-0.79,-1.58]]

## Define ESS PWM arrays (Adapted from: 10.1186/s12915-016-0279-9)
hnRNPA1 = [[-0.13,-0.51,-0.51,0.75],[-0.74,-0.19,-0.92,0.99],
          [1.92,-5.97,-3.11,-3.72],[-1.61,-4.64,1.78,-2.41],
          [-2.13,-5.64,1.80,-1.88],[0.05,-1.49,0.92,-0.48],
          [1.94,-6.38,-4.97,-3.27],[-4.16,-4.16,1.95,-4.97]]

## Define thresholds
SRSF1_threshold = 1.956
SRSF1_igM_threshold = 1.867
SRSF2_threshold = 2.383
SRSF5_threshold = 2.670
SRSF6_threshold = 2.676

NUCLEOTIDE_TO_INDEX = dict(zip("ACGT", range(4))) | dict(zip("acgt", range(4)))

def x(sequence, pssm_list):
    for pssm_matrix, threshold in pssm_list:
        a = find_pssm_matches_with_numpy([list(i) for i in zip(*pssm_matrix)], sequence, threshold)
        print([(i, i + len(pssm_matrix), 1, sequence[i:i+len(pssm_matrix)]) for i in a])

def find_pssm_matches_with_numpy(pssm_matrix, sequence, threshold):
    """Return every index in the +1 strand wit a PSSM score above threshold.

    Adapted from dnachisel

    Parameters:
    -----------

    pssm
      A matrix whose rows give the frequency motif of ATGC (in this order).

    sequence
      A string representing a DNA sequence.

    threshold
      Every index with a score above this threshold will be returned.
    """
    len_pattern = len(pssm_matrix[0])

    # If sequence is small, use normal python to avoid numpy overhead
    nucl_indices = [NUCLEOTIDE_TO_INDEX[n] for n in sequence]
    print(nucl_indices, len(sequence), len_pattern)
    for i in range(len(sequence) - len_pattern + 1):
        print(sequence[i : len_pattern + i], np.choose(nucl_indices[i : len_pattern + i], pssm_matrix).sum())

    if len(sequence) < 60:
        nucl_indices = [NUCLEOTIDE_TO_INDEX[n] for n in sequence]
        return [
            i
            for i in range(len(sequence) - len_pattern + 1)
            if np.choose(nucl_indices[i : len_pattern + i], pssm_matrix).sum()
            >= threshold
        ]
    
    # If sequence is large, use Numpy for speed. tested experimentally

    nucl_indices = np.array([NUCLEOTIDE_TO_INDEX[n] for n in sequence], dtype="uint8")
    len_pattern = len(pssm_matrix[0])
    scores = np.array(
        [
            np.choose(nucl_indices[k : len_pattern + k], pssm_matrix).sum()
            for k in range(len(sequence) - len_pattern)
        ]
    )
    return np.nonzero(scores >= threshold)[0]

def main():
    test_seq0 = 'cttgttataggtggtcc'
    test_seq0 = 'cttgttataggc'
    test_seq1 = 'ggaaaggttactcaga'    
    # test_seq1 = 'ggaaagggtactcaga'    
    test_seq2 = 'cttgttataggtggtccaggaagtggaaagggtactcaga'

    x(test_seq1, [(SRSF1, SRSF1_threshold), (SRSF1_igM, SRSF1_igM_threshold),
                  (SRSF2, SRSF2_threshold), (SRSF5, SRSF5_threshold), (SRSF6, SRSF6_threshold)])

if __name__ == "__main__":
    main()
