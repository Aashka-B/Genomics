# Function to extract reads from the input file
def read_contig_file(file_path):
    """
    Reads an input fasta file of contigs and outputs a dictionary of contig names
    and their corresponding sequences.
    :param file_path:
    :return:
    """

    # Initializing variables
    contig_names = []
    contig_sequences = []
    current_sequence = ""

    with open(file_path) as file:
        for line in file:
            line = line.strip()  # Remove leading/trailing whitespace
            if not line:
                continue  # Skip empty lines

            if line[0] == ">":
                contig_names.append(line[1:].split("|")[0])  # Extract contig name
                if current_sequence:
                    contig_sequences.append(current_sequence)
                current_sequence = ""
            else:
                current_sequence += line

        if current_sequence:
            contig_sequences.append(current_sequence)

        contig_dict = dict(zip(contig_names, contig_sequences))

        return contig_dict


def scoreMotif(sequence):
    """
    Defines the scoring matrix and performs score calculation.
    :param sequence:
    :return:
    """
    # Ensure the sequence is 13 bases long
    if len(sequence) != 13:
        raise ValueError("Sequence must be 13 bases long")

    # Scoring matrix for each base at each position
    base_idx = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    motif = [[.5, .5, .5, .5, 0, 0, 0, 0, 0, 2, -99, -99, .5],  # A
             [0, 0, 0, 0, 0, 0, 0, 0, 0, -99, 2, -99, 0],  # T
             [0, 0, 0, 0, 0, 0, 0, 0, 0, -99, -99, -99, 0],  # C
             [.5, .5, .5, .5, 0, 0, 0, 0, 0, .5, -99, 2, 0]  # G
             ]

    # Calculate the score
    score = 0
    for i, base in enumerate(sequence):
        if base in base_idx:  # Check if the base is one of A, T, C, or G
            base_score = motif[base_idx[base]][i]
            if base_score != -99:  # Check if the position is not to be ignored
                score += base_score
        else:
            raise ValueError(f"Invalid base '{base}' in sequence at position {i + 1}")

    return score


def scoreCodon(codon):
    """
    Score a single codon based on a simplified scoring logic.
    :param codon: A 3-base codon sequence.
    :return: Score of the codon.
    """
    # Simplified scoring logic
    score = 0
    if codon == 'ATG':  # Start codon example
        score = 1
    elif codon == 'GTG':  # Stop codons example
        score = -1
    return score


def scanSeq(contigs_dict):
    potentialStartPos = {}
    ORFlengths = {}
    ORFSeqs = {}

    start_codons = ['ATG', 'GTG']
    stop_codons = ['TAA', 'TAG', 'TGA']

    for contig_name, sequence in contigs_dict.items():
        unique_ORFs = set() # Initialize an empty set to track previously identified ORFs for each contig
        potentialStartPos[contig_name] = []
        ORFlengths[contig_name] = []
        ORFSeqs[contig_name] = []

        # Iterate through the sequence, considering each 13-base segment as a potential motif
        for idx in range(len(sequence) - 12):
            motif = sequence[idx:idx+13]
            # Directly score the motif, assuming it includes a start codon
            if scoreMotif(motif) > 7.25:
                # Identify the start codon within the motif
                max_score = -1
                for start_codon_idx in range(idx, min(idx + 13, len(sequence) - 2), 3):
                    potential_start_codon = sequence[start_codon_idx:start_codon_idx+3]
                    if potential_start_codon in start_codons:
                        codon_score = scoreCodon(potential_start_codon)
                        if codon_score > max_score:
                            max_score = codon_score
                            max_score_start_codon_idx = start_codon_idx
                            # Once a start codon is found, search for the nearest stop codon
                            for j in range(max_score_start_codon_idx + 3, len(sequence) - 2, 3):
                                if sequence[j:j+3] in stop_codons:
                                    orf_seq = sequence[max_score_start_codon_idx:j+3]
                                    break
                                elif j + 3 >= len(sequence) - 2:  # If no stop codon is found
                                    orf_seq = sequence[max_score_start_codon_idx:]
                            orf_length = len(orf_seq)
                            if orf_length >= 60:
                                # Create a unique key for the ORF based on its start position and length
                                orf_key = (max_score_start_codon_idx, orf_length)
                                # Check if the ORF has already been identified
                                if orf_key not in unique_ORFs:
                                    unique_ORFs.add(orf_key)  # Mark this ORF as identified
                                    potentialStartPos[contig_name].append(max_score_start_codon_idx)
                                    ORFlengths[contig_name].append(orf_length)
                                    ORFSeqs[contig_name].append(orf_seq)
                            break  # Proceed to search for the next motif after finding an ORF

    return potentialStartPos, ORFlengths, ORFSeqs


def identifyORFs(input_file_path):
    # Load contigs from the input file
    contigs = read_contig_file(input_file_path)

    # Specify the output file name directly
    output_file_path = "ORFs.fa"

    # Open the output file to write the ORFs
    with open(output_file_path, 'w') as output_file:
        for contig_name, sequence in contigs.items():
            # Retrieve the ORF data for the current contig
            potentialStartPos, ORFlengths, ORFSeqs = scanSeq({contig_name: sequence})

            # Iterate through the identified ORFs for this contig
            for i, orf_seq in enumerate(ORFSeqs[contig_name], start=1):
                orf_start = potentialStartPos[contig_name][i - 1] + 1  # Convert to 1-based index
                orf_length = ORFlengths[contig_name][i - 1]

                # Write the ORF information to the output file
                output_file.write(f"> {contig_name}_ORF{i}|Length {orf_length}|at position {orf_start}\n{orf_seq}\n")



infile = "/Users/aashka.b/ORFid/spaceSeq.fa"
identifyORFs(infile)
