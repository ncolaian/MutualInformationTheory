#!/usr/bin/env python3

"""MIT-run.py

Command-line helper to compute Mutual Information Theory (MIT) scores
between positions in aligned sequences (MSA files). The script provides
utilities to read FASTA MSAs, build position-identity matrices, compute
pair-frequency matrices and calculate a Mutual Information score for
position pairs.

This module is intended to be run as a small CLI tool (see `main`). It
expects two MSAs (or one MSA used for both inputs) and outputs a CSV with
per-position MIT scores.

*NOTE* 
This has a known issue where the amino acid usage bias exists - meaning if a position utilizes 
more diverse amino acids it will tend to have a higher MIT score overall regardless of true co-evolutionary signal.

There are plans to address this in future versions - and can also be mitigated at the analysis stage with z-score normalization.
"""

import os
import sys
import argparse
import pandas as pd
from typing import Dict
from math import log

#some full lists of amino acids
AMINO_ACIDS = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", 
          "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

#all the nucleic acids
NUCLEIC_ACIDS = ["A", "C", "G", "T"]

def get_args() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Calculate MIT scores for Amino Acid sequences or Nucleic Acid sequences"
    )
    parser.add_argument("-m1", "--msa1", required=True, type=str, help="Multiple Sequence Alignment file in FASTA format")
    parser.add_argument("-m2", "--msa2", required=False, type=str, help="Multiple Sequence Alignment file in FASTA format")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output path to write the final MIT scores to. File will be named mit_results.csv")
    parser.add_argument("-t", "--type", required=True, type=str, choices=["A", "N"], default="A", help="Type of sequence: A for Amino Acid, N for Nucleic Acid. Default is A")

    return parser.parse_args()

def get_sequences_from_fasta(fasta_file: str) -> list:
    """Read sequences from a FASTA file and return a list of sequences.

    The function handles multi-line FASTA entries and returns sequence
    strings in the same order they appear in the file. Description lines
    beginning with '>' are ignored.

    Parameters
    ----------
    fasta_file : str
        Path to a FASTA file.

    Returns
    -------
    list[str]
        List of sequence strings (no header lines).
    """
    sequences = []
    with open(fasta_file, 'r') as f:
        seq = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                #check if there was a sequence built out
                if seq:
                    sequences.append(seq)
                    seq = ""
            else:
                seq += line
        if seq:
            sequences.append(seq)
    
    return sequences

def create_position_identity_matrix(sequences: list) -> Dict[int, list]:
    """Build a position identity matrix from aligned sequences.

    Parameters
    ----------
    sequences : list[str]
        Aligned sequences (all sequences should have the same length).

    Returns
    -------
    Dict[int, list]
        Mapping from 1-based position to a list of characters observed at
        that position across all sequences.
    """
    identity_dict = {}
    for seq in sequences:
        for pos, char in enumerate(seq, start=1):
            if pos in identity_dict:
                identity_dict[pos].append(char)
            else:
                identity_dict[pos] = [char]

    return identity_dict

def calculate_identity_pair_frequency_and_MIT(identity_dict1: Dict[int, list], identity_dict2: Dict[int, list], valid_chars: list) -> pd.DataFrame:
    """Compute pair-frequency matrices and MIT scores for position pairs.

    For every position in `identity_dict1` and every position in
    `identity_dict2`, this function counts co-occurring character pairs
    (e.g. 'A-R') across genomes, converts counts to probabilities, and
    calculates an MIT score for the position pair. The result is returned
    as a pandas.DataFrame with columns 'Position_MSA1', 'Position_MSA2',
    and 'MIT_Score'.

    Parameters
    ----------
    identity_dict1, identity_dict2 : Dict[int, list]
        Position identity matrices produced by `create_position_identity_matrix`.
    valid_chars : list[str]
        List of allowed characters (amino-acids or nucleotides).

    Returns
    -------
    pandas.DataFrame
        One row per (pos1,pos2) with the computed MIT score.
    """
    # Here you want to get the sequence that comes from the same genome and the fasta files should already be filtered and ordered as such
    # Initialize a dictionary to hold the pair frequencies
    pair_probabilities = {f"{a}-{b}": 0 for a in valid_chars for b in valid_chars}
    single_probabilities = {f"{a}": 0 for a in valid_chars}

    #initialize a list of dictionaries to hold the results
    results = []
    #go through each position within the sequence inside the associated genome and count up the pairings
    for pos in identity_dict1:
        for pos2 in identity_dict2:
            #create a counting mechanism for all pairs at these positions across all genomes
            pair_counts = pair_probabilities.copy()
            #counting mechanism for the single amino acids for each position
            seq1_pb = single_probabilities.copy()
            seq2_pb = single_probabilities.copy()

            #loop through each amino acid/nucleic acid pairs between the two positions
            #genome position = the matched genome from each sequence - since these should already be ordered
            n = len(identity_dict1[pos])
            #the length of the arrays within each position should be the same because there should be the same number of sequences
            if n != len(identity_dict2[pos2]):
                print(f"Error: The number of sequences at position {pos} in MSA1 does not match the number of sequences at position {pos2} in MSA2.")
                sys.exit(1)

            for genome_pos in range(n):
                #go in each dict and pull the amino acid pair from the same genome
                base1 = identity_dict1[pos][genome_pos]
                base2 = identity_dict2[pos2][genome_pos]

                #check if one of the characters is a gap or not in the valid characters
                #usually for -, X, N, *, etc.
                if base1 not in valid_chars or base2 not in valid_chars:
                    continue

                pair = f"{base1}-{base2}"

                #increment the counts
                pair_counts[pair] += 1
                seq1_pb[base1] += 1
                seq2_pb[base2] += 1
            
            #now convert counts to probabilities and calculate the MIT score
            mit_score = calculate_mit(
                convert_counts_to_probabilities(pair_counts),
                convert_counts_to_probabilities(seq1_pb),
                convert_counts_to_probabilities(seq2_pb),
                valid_chars
            )

            results.append({"Position_MSA1": pos, "Position_MSA2": pos2, "MIT_Score": mit_score})

    # Convert results to a DataFrame
    df = pd.DataFrame(results)

    return df

def convert_counts_to_probabilities(counts: Dict[str, int]) -> Dict[str, float]:
    """Normalize a counts dict into probabilities.

    Parameters
    ----------
    counts : Dict[str, int]
        Mapping from key to integer counts.

    Returns
    -------
    Dict[str, float]
        Mapping from key to fractional probability. If the total count is
        zero, returns zeros for all keys.
    """
    total = sum(counts.values())
    if total == 0:
        return {k: 0.00 for k in counts}
    return {k: v / total for k, v in counts.items()}       

####this portion calculates the MIT score based on the probabilities for each position identity and pair identity
def calculate_mit(pair_probabilities: Dict[str, float], seq1_probabilities: Dict[str, float], 
                  seq2_probabilities: Dict[str, float], valid_chars: list) -> float:
    """Calculate Mutual Information (MIT) score for two positions.

    Parameters
    ----------
    pair_probabilities : Dict[str, float]
        Joint probabilities for pairs like 'A-R'.
    seq1_probabilities, seq2_probabilities : Dict[str, float]
        Marginal probabilities for each character at the position being investigated.
    valid_chars : list[str]
        List of valid characters; used to iterate over possible pairs.

    Returns
    -------
    float
        MIT score (using log base `len(valid_chars)`). Terms with zero
        probability are skipped (they do not contribute).
    """
    mit_score = 0.0
    for a in valid_chars:
        for b in valid_chars:
            pair_key = f"{a}-{b}"
            #there should be something for each key but just in case, use get with a default of 0.0
            p_ab = pair_probabilities.get(pair_key, 0.0)
            p_a = seq1_probabilities.get(a, 0.0)
            p_b = seq2_probabilities.get(b, 0.0)

            #the log is dependent on if it is amino acids or nucleic acids
            #if there are zeroes, skip them because they don't contribute to the score
            if p_ab > 0 and p_a > 0 and p_b > 0:
                mit_score += p_ab * log( p_ab/(p_a * p_b), len(valid_chars) )
    
    return mit_score

def main():
    """Main entry point: parse arguments, compute MIT scores, and write CSV.

    The function reads input MSAs, constructs position identity matrices,
    computes MIT scores for every position pair, and writes a CSV file to
    the directory specified by the `-o/--output` argument.
    """

    # parse the arguments
    args = get_args()
    #check if the input files exist and pull the sequences
    if not os.path.isfile(args.msa1):
        print(f"Error: The file {args.msa1} does not exist.")
        sys.exit(1)
    else:
        sequences1 = get_sequences_from_fasta(args.msa1)
        position_identies1 = create_position_identity_matrix(sequences1)
        if len(sequences1) == 0:
            print(f"Error: No sequences found in {args.msa1}.")
            sys.exit(1)
    #check if the second MSA file exists and pull the sequences, otherwise use the first MSA for both inputs
    if args.msa2 and not os.path.isfile(args.msa2):
        print(f"Error: The file {args.msa2} was passed and does not exist.")
        sys.exit(1)
    elif args.msa2 is None:
        position_identies2 = position_identies1
        print("No second MSA provided, using the first MSA for both inputs to calculate MIT on itself.")
    else:
        sequences2 = get_sequences_from_fasta(args.msa2)
        position_identies2 = create_position_identity_matrix(sequences2)
        if len(sequences2) == 0:
            print(f"Error: No sequences found in {args.msa2}.")
            sys.exit(1)

    #report the type of sequences being processed
    if args.type == "A":
        valid_chars = AMINO_ACIDS
        print("Processing as Amino Acid sequences.")
    else:
        valid_chars = NUCLEIC_ACIDS
        print("Processing as Nucleic Acid sequences.")

    
    mit_results = calculate_identity_pair_frequency_and_MIT(position_identies1, position_identies2, valid_chars)
    
    #write the results to a CSV file
    os.path.isdir(args.output) or os.makedirs(args.output)
    mit_results.to_csv(f"{args.output}/mit_results.csv", index=False, header=True)

    return 1


if __name__ == "__main__":
    main()