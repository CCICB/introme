import pandas as pd
from io import StringIO
import sys

data = """
CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	A1_Hazeem	A1_Hazeem_july	A1_neuBG	A1_winBG	SRSF1	SRSF1_igM	SRSF2	SRSF5	SRSF6
chr1	3244043	.	T	C	.	PASS	strand=+;gene_id=ENST00000270722.10;location=Intronic;MFASS_delta_index=0.00845268;classification=Normal	0	0	0	0	0	0	0	0	0
chr1	3244065	.	C	T	.	PASS	strand=+;gene_id=ENST00000270722.10;location=Intronic;MFASS_delta_index=0.0138251;classification=Normal	0	0	0	0	0	0	0	0	0
chr1	3244069	.	C	T	.	PASS	strand=+;gene_id=ENST00000270722.10;location=Intronic;MFASS_delta_index=0.0137124;classification=Normal	0	0	0	0	0	0	0	0	0
chr1	3244093	.	C	G	.	PASS	strand=+;gene_id=ENST00000270722.10;location=Exonic;MFASS_delta_index=0.0134846;classification=Normal	4.537	4.062	5.26	3.868	-1.577	-4.153	-2.625	0	0
"""

data = open(sys.argv[1], 'r')
output = open(sys.argv[2], 'w+')

# Read the data into a pandas DataFrame
df = pd.read_csv(data, sep='\t')

# Function to expand INFO column
def expand_info(row):
    info_items = row['INFO'].split(';')
    info_dict = {}
    for item in info_items:
        key, value = item.split('=')
        info_dict[key] = value
    return pd.Series(info_dict)

# Apply the function to each row and join the resulting DataFrame with the original one
df_info_expanded = df.apply(expand_info, axis=1)
df_expanded = pd.concat([df, df_info_expanded], axis=1).drop(columns=['INFO'])

df_expanded.to_csv(output, encoding='utf-8', index=False, sep='\t')
