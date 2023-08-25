from Bio.motifs.matrix import PositionSpecificScoringMatrix
import numpy as np

class RBPsplice(PositionSpecificScoringMatrix):
    def __init__(self, alphabet, values, name, *, threshold: float=None):
        super().__init__(alphabet, values)
        self.threshold = threshold
        self.name = name

    def calculate(self, sequence, *, use_threshold=True):
        l = super().calculate(sequence)
        if use_threshold and self.threshold is not None:
            return np.where(l >= self.threshold, l, 0.0)
        else:
            return l

    @classmethod
    def from_2D_list(cls, arr: list[list[float]], name: str, *, arr_in_row_form=True, threshold: float=None):
        """
        Convert PSSM in 2D list form into a Bio.motifs.matrix.PSSM compatible format and instantiate it

        "row form" means list of 4 lists eg [[0.1, 0.6], [-0.4, 2.3], [0.4, 0.1], [-0.9, 1.2]]
        "column form" is the transpose of row form eg [[0.1, -0.4, 0.4, -0.9], [0.6, 2.3, 0.1, 1.2]]
        """

        if not arr_in_row_form:
            # Transpose
            arr = [list(i) for i in zip(*arr)]
        # return cls(["A", "C", "G", "T"], {"A": [0, 1, 0], "C": [0, 0, 1], "G": [1, 0, 0], "T": [0, 0, 0]})
        assert len(arr) == 4, "Is arr in row form?"

        return cls(["A", "C", "G", "T"], {"A": arr[0], "C": arr[1], "G": arr[2], "T": arr[3]}, name, threshold=threshold)

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

RBPmotifs = [(RBPsplice.from_2D_list(mot, name, arr_in_row_form=False, threshold=thr))
             for mot, name, thr in [
                 (SRSF1, "SRSF1", SRSF1_threshold), (SRSF1_igM, "SRSF1_igM", SRSF1_igM_threshold),
                (SRSF2, "SRSF2", SRSF2_threshold), (SRSF5, "SRSF5", SRSF5_threshold), (SRSF6, "SRSF6", SRSF6_threshold)]#,
                # (hnRNPA1, "hnRNPA1", None)]
            ]

if __name__ == "__main__":
    a = RBPsplice.from_2D_list(SRSF1, "SRSF1", arr_in_row_form=False, threshold=SRSF1_threshold)
    print(a.length)
    print(type(x:=a.calculate("GACTACACACTAC")))

    for i in x:
        print(i)

    print(RBPmotifs)

