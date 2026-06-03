import sys


def extract_accessions(hmm_file, num_headers):
    # Output file name: same base name with "_headers.txt"
    output_file = hmm_file.rsplit('.', 1)[0] + "_headers.txt"
    print(f"Writing first {num_headers} accessions to {output_file}")

    count = 0
    with open(hmm_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith(">>"):
                count += 1
                # Strip leading '>>' and whitespace
                header = line[2:].strip()
                # First whitespace-delimited token
                first_token = header.split()[0]
                # If UniProt-style (tr|ACC|ID), extract the accession
                if "|" in first_token:
                    parts = first_token.split("|")
                    if len(parts) >= 2:
                        accession = parts[1]
                    else:
                        accession = first_token
                else:
                    accession = first_token
                outfile.write(accession + "\n")

                if count >= num_headers:
                    break

    print(f"Extracted {count} accessions into {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python hmm_header_get.py <hmmsearch_output.txt> <num_headers>")
        sys.exit(1)
    try:
        num_headers = int(sys.argv[2])
    except ValueError:
        print("Error: <num_headers> must be an integer")
        sys.exit(1)

    extract_accessions(sys.argv[1], num_headers)
