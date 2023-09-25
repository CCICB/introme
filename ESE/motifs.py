from Bio.motifs.matrix import PositionSpecificScoringMatrix
import numpy as np
from variants import VariantContext, StrandDirection
from typing import Optional
import math

class RBPsplice(PositionSpecificScoringMatrix):
    def __init__(self, alphabet, values, name, *, threshold: Optional[float]=None):
        super().__init__(alphabet, values)
        self.threshold = threshold
        self.name = name

    def calculate_variant(self, variant: VariantContext):
        if variant.strand_direction == StrandDirection.FORWARD:
            return round(
                np.max(self._calculate(variant.alt_sequence_fwd(self.length)))
                - np.max(self._calculate(variant.ref_sequence_fwd(self.length)))
            , 3)
        elif variant.strand_direction == StrandDirection.REVERSE:
            return round(
                np.max(self._calculate(variant.alt_sequence_rev(self.length)))
                - np.max(self._calculate(variant.ref_sequence_rev(self.length)))
            , 3)
        raise NotImplementedError("BOTH/NONE strand not supported currently")

    # Warning i changed use_threshold to False from True
    def _calculate(self, sequence, *, use_threshold=True):
        l = super().calculate(sequence)
        if use_threshold and self.threshold is not None:
            return np.where(l >= self.threshold, l, 0.0)
        else:
            return l

    @classmethod
    def from_2D_list(cls, arr: list[list[float]], name: str, *, arr_by_base=True, threshold: Optional[float]=None):
        """
        Convert PSSM in 2D list form into a Bio.motifs.matrix.PSSM compatible format and instantiate it
        For an exemplary 2bp motif:
        "by base" means 1 list per base letter, ie 4 lists total: [[0.1, 0.6], [-0.4, 2.3], [0.4, 0.1], [-0.9, 1.2]]
        "by position" means 1 list of len 4 per base pair postion, ie [[0.1, -0.4, 0.4, -0.9], [0.6, 2.3, 0.1, 1.2]]
        """
        if not arr_by_base:
            # Transpose
            arr = [list(i) for i in zip(*arr)]
        assert len(arr) == 4, "Check whether you really gave arr by base letter?"

        return cls(["A", "C", "G", "T"], {"A": arr[0], "C": arr[1], "G": arr[2], "T": arr[3]}, name, threshold=threshold)
    
    @classmethod
    def from_RCRUNCH(cls, arr: list[list[float]], name: str, *, arr_by_base=True, threshold: Optional[float]=None):
        """
        Convert IC matrix in 2D list from into a Bio.motifs.matrix.PSSM compatible format and instantiate it.
        For an exemplary 2bp motif:
        "by base" means 1 list per base letter, ie 4 lists total: [[0.1, 0.6], [-0.4, 2.3], [0.4, 0.1], [-0.9, 1.2]]
        "by position" means 1 list of len 4 per base pair postion, ie [[0.1, -0.4, 0.4, -0.9], [0.6, 2.3, 0.1, 1.2]]
        """
        if arr_by_base:
            # Transpose
            arr = [list(i) for i in zip(*arr)]
        assert len(arr[0]) == 4, "Check whether you really gave arr by postion?"

        # IC matrix is just frequency * information content of that position
        ic_matrix = arr
        freq_matrix = []
        # Iterate through each position
        for ic_by_basepair in ic_matrix:
            ic_by_basepair = [x + (0/100 * sum(ic_by_basepair)) for x in ic_by_basepair]
            total_ic_in_position = sum(ic_by_basepair)
            freq_matrix.append([base_letter_ic / total_ic_in_position for base_letter_ic in ic_by_basepair])

        print(np.round(np.array(freq_matrix), 5))
        # TODO: add functionality to have non-uniform background
        bg = [0.25, 0.25, 0.25, 0.25]
        log_odds = []
        # Iterate through each position
        for freq_by_basepair in freq_matrix:
            log_odds.append([math.log(base_letter_freq / bg[i], 2) for i, base_letter_freq in enumerate(freq_by_basepair)])

        # for row in arr:
        #     print([round(x, 2) for x in row])
        # print()
        # for row in log_odds:
        #     print([round(x, 2) for x in row])

        # Transpose into by base form
        log_odds = [list(i) for i in zip(*log_odds)]
        return cls(["A", "C", "G", "T"], {"A": log_odds[0], "C": log_odds[1], "G": log_odds[2], "T": log_odds[3]}, name, threshold=threshold)
    

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
            [0.43,-0.99,0.0,0.2]]

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
#TODO: scrape/reverse engineer from ESEfinder and get threshold
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

#TODO: obsolete
_RBPmotifs = [(RBPsplice.from_2D_list(mot, name, arr_by_base=False, threshold=thr))
             for mot, name, thr in [
                 (SRSF1, "SRSF1", SRSF1_threshold), (SRSF1_igM, "SRSF1_igM", SRSF1_igM_threshold),
                (SRSF2, "SRSF2", SRSF2_threshold), (SRSF5, "SRSF5", SRSF5_threshold), (SRSF6, "SRSF6", SRSF6_threshold)]#,
                # (hnRNPA1, "hnRNPA1", None)]
            ]

if __name__ == "__main__":
    # a = RBPsplice.from_2D_list(SRSF1, "SRSF1", arr_by_base=False, threshold=SRSF1_threshold)
    # print(a.length)
    # print(type(scores:=a.calculate("GACTACACACTAC")))

    # for i in scores:
    #     print(i)

    # print("SRSF1 motif:")
    # print(RBPmotifs)
    # print(a.__str__())

    NIPBL = [
        [0.0001994784622723167, 1.9949841011854392, 0.0001994784622723167, 0.0001994784622723167],
        [0.06152484191379193, 0.2458965607924698, 0.36880427981616504, 6.760228756597289e-05],
        [0.029382570353052955, 0.11743342643906694, 0.14678371180107158, 0.029382570353052955],
        [0.00011513443947882394, 0.3141212912300754, 0.8374533724371217, 0.00011513443947882394],
        [0.00010511383648029525, 0.382341068813426, 0.6690075236624872, 0.00010511383648029525],
        [0.0001994784622723167, 1.9949841011854392, 0.0001994784622723167, 0.0001994784622723167],
        [0.00019947846227231673, 0.00019947846227231673, 1.9949841011854395, 0.00019947846227231673]
    ]

    DROSHA = [
        [0.00019947846227231673, 0.00019947846227231673, 1.9949841011854395, 0.00019947846227231673],
        [0.00010923782761831944, 0.00010923782761831944, 0.7429264656321906, 0.3496702862062405],
        [0.08157492728835188, 0.2069908224272406, 0.2602968302705404, 5.4869797059495445e-05],
        [0.00011412587954497888, 0.31966658860548586, 0.00011412587954497888, 0.8218204586033929],
        [0.00019947846227231673, 0.00019947846227231673, 1.9949841011854395, 0.00019947846227231673],
        [6.943975177383666e-05, 0.05165623134455708, 0.3452822217202254, 0.29766738392890557],
        [0.00019447860498534043, 0.011299206949648281, 1.9338758001137266, 0.00019447860498534043]
    ]

    NKRF = [
        [0.00019947846227231673, 0.00019947846227231673, 1.9949841011854395, 0.00019947846227231673],
        [0.0001994784622723167, 1.9949841011854392, 0.0001994784622723167, 0.0001994784622723167],
        [0.00019947846227231673, 0.00019947846227231673, 1.9949841011854395, 0.00019947846227231673],
        [0.00019947846227231673, 0.00019947846227231673, 1.9949841011854395, 0.00019947846227231673],
        [0.00010019614960339405, 0.458938443643386, 0.5432234446897611, 0.00010019614960339405],
        [0.00019947846227231673, 0.00019947846227231673, 1.9949841011854395, 0.00019947846227231673]
    ]

    # b = RBPsplice.from_RCRUNCH(NIPBL, "NIPBL", arr_by_base=False)
    # c = RBPsplice.from_RCRUNCH(DROSHA, "DROSHA", arr_by_base=False)
    d = RBPsplice.from_RCRUNCH(NKRF, "NKRF", arr_by_base=False)

    # print(b.__str__())
    # print(c.__str__())
    print(d.__str__())
    print(res:=d.calculate("ACCGCGCACCTGGCTTGTTTGTTGCATTTCATAGCAAGTGTCTGATTGCTTCTTTTTTCAGATGTTCACTGCCTTCTTCGGCAGTTCTGTTTTATTGTTATTTATGTTCTCAGTGTTTATTCTTCTTTTCCTTTTGAATGCTATGGCCTTTTAGATACACAGTGACTTTTTCCTTTGTGGCT", use_threshold=False))
    print(sorted(res)[-5:])