import re
import os

def extract_trim_stats(file_content):
    # Define the regex pattern to extract relevant information
    pattern = re.compile(r'Input Read Pairs: (\d+) Both Surviving: (\d+ \([\d.]+%)\) Forward Only Surviving: (\d+ \([\d.]+%)\) Reverse Only Surviving: (\d+ \([\d.]+%)\) Dropped: (\d+ \([\d.]+%)\)')

    # Search for the pattern in the file content
    match = pattern.search(file_content)

    if match:
        # Extract matched groups
        input_pairs = match.group(1)
        both_surviving = match.group(2)
        forward_surviving = match.group(3)
        reverse_surviving = match.group(4)
        dropped = match.group(5)

        # Return formatted result
        return f"{input_pairs}, {both_surviving}, {forward_surviving}, {reverse_surviving}, {dropped},"

    return None

def process_txt_file(file_path):
    with open(file_path, 'r') as file:
        file_content = file.read()

    # Extract trim stats from the file content
    trim_stats = extract_trim_stats(file_content)

    if trim_stats:
        # Extract filename without extension
        filename = os.path.splitext(os.path.basename(file_path))[0]
        return f"{filename}, {trim_stats}"
    else:
        return f"Error: No match found in {file_path}"

def main():
    # Get all .txt files in the current directory
    file_paths = [file for file in os.listdir() if file.endswith(".txt")]

    # Output file
    output_file_path = "result_output.txt"

    with open(output_file_path, 'w') as output_file:
        # Write header
        output_file.write("File_name, Input Read Pairs, Both Surviving, Forward Only Surviving, Reverse Only Surviving, Dropped,\n")

        # Process each file and write the result to the output file
        for file_path in file_paths:
            result = process_txt_file(file_path)
            output_file.write(result + '\n')

    print(f"Results written to {output_file_path}")

if __name__ == "__main__":
    main()

