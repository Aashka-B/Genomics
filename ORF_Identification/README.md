# ORF Identification Tool

## Description
This Python script identifies Open Reading Frames (ORFs) within DNA sequences stored in FASTA format. It utilizes scoring matrices for motifs and codons to find potential start sites and distinguishes between start and stop codons to accurately identify ORFs. The output is a FASTA file containing the identified ORFs, including their lengths and positions.

## Functions
### read_contig_file()
- This function reads DNA sequences from a FASTA file.
- The function processes the file to extract contig names (sequence identifiers) and their corresponding sequences, storing them in a dictionary where keys are the contig names, and values are the sequences.

### scoreMotif()
- Scores a given sequence (motif) based on a predefined scoring matrix for each base at specific positions.
- The function checks if the sequence is exactly 13 bases long and contains valid nucleotide bases (A, T, C, G).
- Each position in the sequence is scored according to how well it matches the expected pattern for an ORF start site.
- The scores are used to evaluate the likelihood of a sequence being a functional start site for transcription.

### scoreCodon()
- This function assigns specific scores for start codons based on the existing scoring matrix for motifs.

### scanSeq()
- Scans sequences for potential ORFs based on the scores from scoreMotif and scoreCodon.
- It looks for high-scoring motifs as potential ORF start sites, then checks for valid start codons within those motifs.
- Once a start codon is identified, it searches for the nearest stop codon to delineate the ORF.
- The function keeps track of unique ORFs identified within each contig, storing their start positions, lengths, and sequences.

### identifyORFs()
- The main function that ties everything together.
- It reads the input FASTA file to get the contig sequences, then uses scanSeq to identify ORFs within those sequences.
- Finally, it writes the identified ORFs to an output file, ORFs.fa, formatting each ORF's information in FASTA format.

## Installation
Ensure you have Python 3.6 or later installed on your system. No additional dependencies are required, as the script uses the Python Standard Library.

## Usage
To use this script, follow these steps:

1. Place your DNA sequence file in FASTA format in a known directory.
2. Run the script from the command line or terminal with the path to your FASTA file as an argument:

    ```bash
    python orf_identifier.py /path/to/your/input_file.fasta
    ```

3. The script will output a file named `ORFs.fa` in the same directory as the input file, containing the identified ORFs.

## Dependencies
- Python 3.6 or higher

## Contributing
Contributions to this project are welcome. To contribute:

1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Commit your changes.
4. Push to the branch.
5. Submit a pull request.

Please ensure your code adheres to the project's coding standards and include comments explaining your changes.

## License
This project is open-sourced under the MIT License. See the LICENSE file for more details.

## Contact
For questions or feedback regarding this project, please contact Aashka Bhowmick at aashka.bhowmick@gmail.com.
