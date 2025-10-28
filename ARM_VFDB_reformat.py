import os

# Function to process a single TSV file
def process_tsv(file_path):
    # Read the TSV file and filter data
    filtered_data = {}
    with open(file_path, 'r') as file:
        next(file)  # Skip header
        for line in file:
            fields = line.strip().split('\t')
            gene_fraction = float(fields[3])
            # Filter out rows with gene fraction <= 99
            if gene_fraction > 99:
                gene = fields[1]
                hits = int(fields[2])  # Convert Hits column to integer
                filtered_data[gene] = hits

    return filtered_data

# Function to process all TSV files in a folder
def process_folder(folder_path):
    # Dictionary to store filtered data for each sample
    sample_data = {}
    
    # Iterate through each file in the folder
    for file_name in os.listdir(folder_path):
        if file_name.endswith('.tsv'):
            file_path = os.path.join(folder_path, file_name)
            # Process the TSV file
            sample_hits = process_tsv(file_path)
            # Store hits for this sample
            sample_data[file_name] = sample_hits
    
    return sample_data

# Function to generate the output format
def generate_output(sample_data):
    # Extract gene names from filtered data
    genes = set()
    for data in sample_data.values():
        genes.update(data.keys())

    # Generate output header
    output = "Gene\t" + '\t'.join(sample_data.keys()) + '\n'

    # Generate output rows
    for gene in sorted(genes):
        row = [gene]
        for sample_hits in sample_data.values():
            hits = str(sample_hits.get(gene, '0'))  # Convert hits to string
            row.append(hits)
        output += '\t'.join(row) + '\n'

    return output

# Main function
def main():
    folder_path = '/home/user/SR/DU-DMS/DMS-Results/AMR-DMS/'  # Update with the folder path containing TSV files
    sample_data = process_folder(folder_path)
    output = generate_output(sample_data)
    
    # Write output to a TSV file
    output_file_path = '/home/user/SR/DU-DMS/DMS-Results/AMR-DMS/99_output.tsv'  # Update with the desired output file path
    with open(output_file_path, 'w') as output_file:
        output_file.write(output)

# Execute the main function
if __name__ == "__main__":
    main()
