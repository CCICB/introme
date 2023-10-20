import csv, sys

# Define the path to the CSV file
file_path = sys.argv[1]

# Lists to store the headers and their corresponding values from the second row
headers = []
values = []

# Open the CSV file and read its content
with open(file_path, 'r') as file:
    reader = csv.reader(file)
    
    # Extract headers
    headers = [col.strip() for col in next(reader)]
    
    # Extract values from the second row
    values = [col.strip() for col in next(reader)]

# Get columns with a value of '1' in the second row
columns_with_1 = [header for header, value in zip(headers, values) if int(value) <= 10]

# Print the result
print(columns_with_1)
