import pprint as pp
from typing import TextIO
import hashlib

def parse_rbp_data(file_obj: TextIO):
    text = file_obj.read()

    # Split the data by each PWM
    sections = [section.strip("//") for section in text.strip().split("////")]

    rbp_data = {}
    for section in sections:
        # Split the section into lines
        lines = section.strip().split("\n")
        
        # Check if there's actual content in the section
        if len(lines) <= 2:  # Only RNA binding protein name and ACGT line
            continue
        
        # The first line is the RNA binding protein name
        rbp_name = lines[0]
        
        # The values start from the third line
        values = []
        for line in lines[2:]:
            # Check if the line contains numeric values
            if any(char.isdigit() for char in line):
                values.append(list(map(float, line.split("\t"))))
        
        keyname = f"{rbp_name}_{len(values)}"
        if keyname not in rbp_data:
            rbp_data[keyname] = values
        else:
            # Check if the old and new motif are the same
            current_hash = hash_2d_list(values)
            previous_hash = hash_2d_list(rbp_data[keyname])
            if current_hash != previous_hash:
                print(f"Detected two different score sets for RBP: {keyname}")
                assert f"{keyname}_b" not in rbp_data
                rbp_data[f"{keyname}_b"] = values
    return rbp_data

def hash_2d_list(two_d_list):
    """Compute a hash for a 2D list."""
    m = hashlib.md5()
    for row in two_d_list:
        m.update(''.join(str(e) for e in row).encode('utf-8'))
    return m.hexdigest()

if __name__ == "__main__":
    import sys
    import json
    with open(sys.argv[1], 'r') as file, open(sys.argv[2], 'w+') as out:
        RBPs = parse_rbp_data(file)
        # pp.pprint(RBPs, out)
        out.write(json.dumps(RBPs, indent=4))




