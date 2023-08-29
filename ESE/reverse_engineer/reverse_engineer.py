from collections import defaultdict
from collections import namedtuple
import sys, pprint

# Define the named tuple
SeqScore = namedtuple('SeqScore', ['sequence', 'score'])

def read_input(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # Skip the header line
    lines = lines[1:]
    
    data = []
    for line in lines:
        fields = line.strip().split('\t')
        seqscore = SeqScore(fields[3], float(fields[4]))
        data.append(seqscore)
    
    return data

def differ_by_one(str1, str2):
    # Check if two strings differ by exactly one character
    if len(str1) != len(str2):
        return None
    
    diff_count = 0
    diff_index = None
    for i, (c1, c2) in enumerate(zip(str1, str2)):
        if c1 != c2:
            diff_count += 1
            diff_index = i
            if diff_count > 1:
                return None
    
    if diff_count == 1:
        return diff_index
    else:
        return None

def calculate_pssm(data):
    pssm = defaultdict(lambda: defaultdict(float))
    
    for i in range(len(data) - 1):
        seq1, score1 = data[i]
        print(tuple([seq1, score1]))
        seq2, score2 = data[i + 1]
        
        # Find the position where the sequences differ
        for pos in range(len(seq1)):
            if seq1[pos] != seq2[pos]:
                diff_pos = pos
                break
        
        # Calculate the score difference
        score_diff = score2 - score1
        
        # Update the PSSM
        pssm[diff_pos][seq2[diff_pos]] = score_diff
    
    return pssm

def differences_by_index(tuple_list: list):    
    # Your list of tuples
    # tuple_list = [
    #     ('AGGGGG', -0.318),
    #     ('GGGGGG', -0.318),
    #     ('AAAAAA', -3.782),
    #     ('AAAAAU', -4.3),
    #     ('AAAAUU', -3.078)
    # ]

    # Sort and uniq the list by the first element of each tuple
    tuple_list = list(set(tuple_list))
    tuple_list.sort(key=lambda x: x.sequence)
    pprint.pprint(tuple_list)

    d = defaultdict(set)
    # Find tuples where the first element differs by one and one character only
    for i in range(len(tuple_list)):
        for j in range(i+1, len(tuple_list)):
            diff_index = differ_by_one(tuple_list[i].sequence, tuple_list[j].sequence)
            if diff_index is not None:
                # print(f"{tuple_list[i]} and {tuple_list[j]} differ by one character at index {diff_index}.")
                # print(f"{tuple_list[i].sequence[diff_index]} is higher than {tuple_list[j].sequence[diff_index]} by "
                #       f"{tuple_list[i].score - tuple_list[j].score:.3f}")
                d[diff_index].add(tuple(i for i in (
                    tuple_list[i].sequence[diff_index],
                    tuple_list[j].sequence[diff_index],
                    round(tuple_list[i].score - tuple_list[j].score, 3)
                    )))
                
    return d

def convert_dict_to_list(d: dict):
    converted_data = []

    for key, values in d.items():
        for value in values:
            concat_str = f"{value[0]}_{key}-{value[1]}_{key}={value[2]}"
            converted_data.append(concat_str)

    return converted_data

def make_equation_from_seqscore(seq_score: SeqScore):
    # Create a list where each item is a character from the sequence appended with "_" and its index
    char_index_list = [f'{char}_{i}' for i, char in enumerate(seq_score.sequence)]

    # Join the list with "+" and append "=" and the score
    result_string = '+'.join(char_index_list) + f'={seq_score.score}'

    return result_string

def main():
    data = read_input(sys.argv[1])
    diff_by_index = differences_by_index(data)

    print()
    pprint.pprint(diff_by_index)

    equations = convert_dict_to_list(diff_by_index)
    for seq_score in data:
        equations.append(make_equation_from_seqscore(seq_score))

    print()
    pprint.pprint(equations)

    from sympy import Eq, solve
    from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application

    # eqs = ['2w + x + 4y + 3z = 5',
    #     'w - 2x + 3z = 3',
    #     '3w + 2x - y + z = -1',
    #     '4x - 5z = -3']
    eqs = equations
    equations.append('N_0 = 0')
    equations.append('N_1 = 0')
    equations.append('N_2 = 0')
    equations.append('N_3 = 0')
    equations.append('N_4 = 0')
    equations.append('N_5 = 0')

    transformations=(standard_transformations + (implicit_multiplication_application,))
    eqs_sympy = [Eq(parse_expr(e.split('=')[0], transformations=transformations),
                    parse_expr(e.split('=')[1], transformations=transformations))
                for e in eqs]
    sol = solve(eqs_sympy)
    # if not dict (list instead), then no solutions
    for k, v in sol.items():
        sol[k] = round(v, 3)
    print(sol)

    # Initialize empty lists for 'A', 'C', 'G', and 'U'
    list_A = []
    list_C = []
    list_G = []
    list_U = []

    # Iterate through the dictionary
    for key, value in sol.items():
        key = str(key)
        # Skip keys that start with 'N'
        if key.startswith('N'):
            continue
        # Add values to the appropriate list
        elif key.startswith('A'):
            list_A.append((key, value))
        elif key.startswith('C'):
            list_C.append((key, value))
        elif key.startswith('G'):
            list_G.append((key, value))
        elif key.startswith('U'):
            list_U.append((key, value))

    # Sort the lists based on the keys
    list_A.sort(key=lambda x: x[0])
    list_C.sort(key=lambda x: x[0])
    list_G.sort(key=lambda x: x[0])
    list_U.sort(key=lambda x: x[0])

    # Extract the values from the sorted lists
    list_A_values = [x[1] for x in list_A]
    list_C_values = [x[1] for x in list_C]
    list_G_values = [x[1] for x in list_G]
    list_U_values = [x[1] for x in list_U]

    PSSM = [list_A_values, list_C_values, list_G_values, list_U_values]
    pprint.pprint(PSSM)
    import os
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import motifs
    import math


    mot = motifs.RBPsplice.from_2D_list(PSSM, "A1_neuBG")
    for seq_score in data:
        if 'N' in seq_score.sequence:
            continue
        recovered_score = mot.calculate(seq_score.sequence.replace('U', 'T'))
        assert math.isclose(recovered_score, seq_score.score, abs_tol=0.0001), [seq_score.sequence, recovered_score, seq_score.score]
    print("Sanity tests worked!")

if __name__ == "__main__":
    main()
