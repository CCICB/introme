import sys
import pysam
import csv
from Bio.Seq import Seq

## Bring in arguments
file = sys.argv[1]
output = sys.argv[2]
reference_genome = sys.argv[3]


## Define ESE PWM arrays (Source: ESEFinder)
SRSF1 = [[-1.14,1.37,-0.21,-1.58],[0.62,-1.1,0.17,-0.5],
		[-1.58,0.73,0.48,-1.58],[1.32,0.33,-1.58,-1.13],
		[-1.58,0.94,0.33,-1.58],[-1.58,-1.58,0.99,-1.13],
		[0.62,-1.58,-0.11,0.27]]

SRSF1_igM = [[-1.58,1.55,-1.35,-1.55],[0.15,-0.53,0.44,-0.28],
			[-0.97,0.79,0.41,-1.28],[0.74,0.33,-0.98,-0.92],
			[-1.19,0.72,0.51,-1.09],[-0.75,-0.62,1.03,-0.52],
			[0.43,-0.99,0,0.2]]

SRSF2 = [[-0.88,-1.16,0.87,-1.18],[0.09,-1.58,0.45,-0.2],
		[-0.06,0.95,-1.36,0.38],[-1.58,1.11,-1.58,0.88],
		[0.09,0.56,-0.33,-0.2],[-0.41,0.86,-0.05,-0.86],
		[-0.06,0.32,-1.36,0.96],[0.23,-1.58,0.68,-1.58]]

SRSF5 = [[-0.13,0.56,-1.58,0.92],[-1.58,0.68,-0.14,0.37],
		[1.28,-1.12,-1.33,0.23],[-0.33,1.24,-0.48,-1.14],
		[0.97,-0.77,-1.58,0.72],[-0.13,0.13,0.44,-1.58],
		[-1.58,-0.05,0.8,-1.58]]

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

## Define functions
# Score all sequence options for the sequence
def seq_scan(seq, motif, strand):
	highest_score = 0
	if strand == "+":
		start = 0
		limit = 1
	elif strand == "-":
		start = 1
		limit = 2
	else:
		start = 0
		limit = 2

	for n in range(start,limit):
		# Swap to reverse strand if strand is negative or both
		if n == 1:
			seq = seq.reverse_complement()

		if motif == "SRSF1":
			for i in range(1,len(seq)-1):
				split = seq_split(seq,7,i)
				if split != None:
					score = score_seq(split,"SRSF1")
					if score > highest_score and score > SRSF1_threshold:
						highest_score = score
		elif motif == "SRSF1_igM":
			for i in range(1,len(seq)-1):
				split = seq_split(seq,7,i)
				if split != None:
					score = score_seq(split,"SRSF1_igM")
					if score > highest_score and score > SRSF1_igM_threshold:
						highest_score = score
		elif motif == "SRSF2":
			for i in range(0,len(seq)):
				split = seq_split(seq,8,i)
				if split != None:
					score = score_seq(split,"SRSF2")
					if score > highest_score and score > SRSF2_threshold:
						highest_score = score
		elif motif == "SRSF5":
			for i in range(1,len(seq)-1):
				split = seq_split(seq,7,i)
				if split != None:
					score = score_seq(split,"SRSF5")
					if score > highest_score and score > SRSF5_threshold:
						highest_score = score
		elif motif == "SRSF6":
			for i in range(2,len(seq)-2):
				split = seq_split(seq,6,i)
				if split != None:
					score = score_seq(split,"SRSF6")
					if score > highest_score and score > SRSF6_threshold:
						highest_score = score
		elif motif == "hnRNPA1":
			for i in range(2,len(seq)-2):
				split = seq_split(seq,6,i)
				if split != None:
					score = score_seq(split,"hnRNPA1")
					if score > highest_score:
						highest_score = score

	return highest_score

# Split the longer sequence into chuncks to be scored
def seq_split(seq, length, start):
	for x in range(start,len(seq)-length):
		return seq[x:x+length]

# Apply the position weight scoring matrix to the sequence
def score_seq(seq, motif):
	strength = 0
	if motif == "SRSF1":
		for pos in range(0,7):
	 		strength += SRSF1[pos][seq_to_array(seq[pos])]
	 	return strength
	elif motif == "SRSF1_igM":
		for pos in range(0,7):
	 		strength += SRSF1_igM[pos][seq_to_array(seq[pos])]
	 	return strength
	elif motif == "SRSF2":
		for pos in range(0,8):
	 		strength += SRSF2[pos][seq_to_array(seq[pos])]
	 	return strength
	elif motif == "SRSF5":
		for pos in range(0,7):
	 		strength += SRSF5[pos][seq_to_array(seq[pos])]
	 	return strength
	elif motif == "SRSF6":
		for pos in range(0,6):
	 		strength += SRSF6[pos][seq_to_array(seq[pos])]
		return strength
	elif motif == "hnRNPA1":
		for pos in range(0,6):
	 		strength += hnRNPA1[pos][seq_to_array(seq[pos])]
		return strength

# Convert the base to the corresponding position in the array
def seq_to_array(base):
	if base == "A":
		return int(0)
	elif base == "C":
		return int(1)
	elif base == "G":
		return int(2)
	elif base == "T":
		return int(3)

################################

with open(file) as read:
	with open(output, 'w') as write:
		writer = csv.writer(write, delimiter='\t')
	  	reader = csv.reader(read, delimiter='\t')

	  	# Loop through each entry
	  	for row in reader:
	  		# Extract necessary variables from file
	  		chromosome = str(row[0])

			# Write header line
			if chromosome[0] == "#":
				row.extend(['SRSF1','SRSF1_igM','SRSF2','SRSF5','SRSF6', 'hnRNPA1'])
				writer.writerow(row)
				continue

			# Extract relevant information from annotated tsv
	  		position = int(row[1])
	  		ref = str(row[3])
	  		alt = str(row[4])
	  		strand = str(row[8])

	  		# Extract reference sequence
	  		# Skip insertions longer than 30bp
	  		if len(alt) > 30:
	  			row.extend(['.','.','.','.','.','.'])
	  			writer.writerow(row)
	  			continue
	  		# If SNV, extract 7bp either side of reference base
			elif len(ref) == 1:
				genome_ref8 = pysam.faidx(reference_genome, chromosome + ":" + str(position-7) + "-" + str(position+7)).split()
				ref8 = Seq(genome_ref8[1])
				alt8 = ref8[:7] + alt + ref8[len(ref8)-7:]
			# If indel, extract 7bp either side of the reference sequence
			else:
				genome_ref8 = pysam.faidx(reference_genome, chromosome + ":" + str(position-7) + "-" + str(position+len(ref)+6)).split()
				ref8 = Seq(genome_ref8[1])
				alt8 = ref8[:8] + ref8[len(ref8)-7:]

			# Score each ref and alt sequence, find the difference and append to the original row
			row.append(seq_scan(ref8, "SRSF1", strand) - seq_scan(alt8, "SRSF1", strand))
			row.append(seq_scan(ref8, "SRSF1_igM", strand) - seq_scan(alt8, "SRSF1_igM", strand))
			row.append(seq_scan(ref8, "SRSF2", strand) - seq_scan(alt8, "SRSF2", strand))
			row.append(seq_scan(ref8, "SRSF5", strand) - seq_scan(alt8, "SRSF5", strand))
			row.append(seq_scan(ref8, "SRSF6", strand) - seq_scan(alt8, "SRSF6", strand))
			row.append(seq_scan(ref8, "hnRNPA1", strand) - seq_scan(alt8, "hnRNPA1", strand))

			# Write the new row to file
			writer.writerow(row)
