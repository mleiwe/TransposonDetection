import argparse
import csv
from collections import Counter

def read_fastq_file(filename):
    sequences = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        for i in range(1, len(lines), 4):
            sequences.append(lines[i].strip())
    return sequences

def search_sequences(sequences, target_sequence):
    matched_indices = []
    for idx, seq in enumerate(sequences):
        if target_sequence in seq:
            matched_indices.append(idx)
    return matched_indices

def extract_subsequences(sequences, matched_indices, flank_length):
    extracted_subseq = []
    for idx in matched_indices:
        if len(sequences[idx]) >= flank_length * 2:
            left_flank = sequences[idx][flank_length - 25:flank_length]
            right_flank = sequences[idx][flank_length:flank_length + 25]
            subseq = sequences[idx][flank_length - 5:flank_length + 5]
            extracted_subseq.append((left_flank, right_flank, subseq))
    return extracted_subseq

def calculate_percentage_abundance(matches):
    total_matches = len(matches)
    counter = Counter(matches)
    result = []
    for match, count in counter.items():
        percentage = (count / total_matches) * 100
        result.append((match, count, percentage))
    return result

def write_to_csv(filename, data, headers):
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(headers)
        for row in data:
            writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(description="Process FASTQ file and sequences.")
    parser.add_argument("fastq_file", help="Path to the FASTQ file")
    parser.add_argument("right_end", help="15 bp sequence to match (RightEnd)")
    parser.add_argument("left_end", help="15 bp sequence to match (LeftEnd)")
    args = parser.parse_args()

    sequences = read_fastq_file(args.fastq_file)
    right_end_indices = search_sequences(sequences, args.right_end)
    left_end_indices = search_sequences(sequences, args.left_end)

    right_end = [sequences[idx] for idx in right_end_indices]
    left_end = [sequences[idx] for idx in left_end_indices]

    right_end_trim = extract_subsequences(right_end, list(range(len(right_end))), flank_length=15)
    left_end_trim = extract_subsequences(left_end, list(range(len(left_end))), flank_length=15)

    right_end_tsd = [(seq[0][-5:], seq[1][:5]) for seq in right_end_trim]
    left_end_tsd = [(seq[0][-5:], seq[1][:5]) for seq in left_end_trim]

    write_to_csv("RightEndTrim.csv", right_end_trim, ["LeftFlank", "RightFlank", "Subsequence"])
    write_to_csv("LeftEndTrim.csv", left_end_trim, ["LeftFlank", "RightFlank", "Subsequence"])
    write_to_csv("RightEndTSD.csv", right_end_tsd, ["LeftTSD", "RightTSD"])
    write_to_csv("LeftEndTSD.csv", left_end_tsd, ["LeftTSD", "RightTSD"])

    right_end_percentage = calculate_percentage_abundance([seq[2] for seq in right_end_trim])
    left_end_percentage = calculate_percentage_abundance([seq[2] for seq in left_end_trim])

    write_to_csv("RightEndData.csv", right_end_percentage, ["Subsequence", "Count", "Percentage"])
    write_to_csv("LeftEndData.csv", left_end_percentage, ["Subsequence", "Count", "Percentage"])
    
    write_to_csv("RightEnd.csv", [[seq] for seq in right_end], ["Sequence"])
    write_to_csv("LeftEnd.csv", [[seq] for seq in left_end], ["Sequence"])

if __name__ == "__main__":
    main()
