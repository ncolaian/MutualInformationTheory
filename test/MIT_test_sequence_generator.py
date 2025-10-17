#!/usr/bin/env python3
"""
This script generates test sequences for unit tests in the MIT (Mutual Information Theory).
It is designed to be very simple to allow easy verification of correctness and can be used to test
the performance of the MIT calculations single and paired sequence alignments (MSAs).
"""


import os
import sys
import argparse
import pandas as pd
from typing import Dict
from math import log
import random

#some full lists of amino acids
AMINO_ACIDS = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", 
          "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

def get_args() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Creates Amino Acid or Nucleic Acid sequences that have \"coevolving\" positions."
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        required=True,
        help="Output directory for the generated sequences."
    )
    parser.add_argument(
        "-n", "--num_sequences",
        type=int,
        default=1000,
        help="Number of sequences to generate. Default is 1000."
    )
    parser.add_argument(
        "-s", "--single",
        action="store_true",
        help="Generate a single MSA with coevolving positions instead of a pair."
    )

    return parser.parse_args()

def generate_sequences(num_sequences: int, seq_length: 20) -> Dict[str, list]:
    #sequences will be 100 bases long
    sequences = []
    seq = "A"*seq_length
    for i in range(num_sequences):
        sequences.append(seq)

    return sequences

def create_changes_paired(sequence1: list, sequence2: list, positions: list, change_frequency: int, same_character_freq: int, valid_chars: list) -> list:
    """Create changes at specified positions in the sequences.

    Args:
        sequences (list): List of sequences to modify.
        positions (list): List of positions to modify in each sequence.
        change_frequency (int): Frequency of changes (e.g., every Nth sequence).
        same_character_freq (int): Frequency of keeping the same character.
        valid_chars (list): List of valid characters to use for changes.

    Returns:
        list: Modified list of sequences.
    """
    
    #these are static characters to change things to
    seq1_change = random.choice(valid_chars)
    seq2_change = random.choice(valid_chars)

    #ensure there is change
    if seq1_change == "A":
        seq1_change = "L"
        
    if seq2_change == "A":
        seq2_change = "N"
    
    for i in range(len(sequence1)):
        change = random.randint(1, 100)
        #this will change the sequences together
        if change >= change_frequency:
            same_change = random.randint(1, 100)
            #change both positions to the new characters but "consistent characters"
            if same_change >= same_character_freq:
                sequence1[i] = sequence1[i][:positions[0]] + seq1_change + sequence1[i][positions[0]+1:]
                sequence2[i] = sequence2[i][:positions[1]] + seq2_change + sequence2[i][positions[1]+1:]
            else:
                #create random changes at the positions 
                seq1_change = random.choice(valid_chars)
                seq2_change = random.choice(valid_chars)
                sequence1[i] = sequence1[i][:positions[0]] + seq1_change + sequence1[i][positions[0]+1:]
                sequence2[i] = sequence2[i][:positions[1]] + seq2_change + sequence2[i][positions[1]+1:]
                

    return sequence1, sequence2

def create_changes_single(sequence1: list, positions: list, change_frequency: int, same_character_freq: int, valid_chars: list) -> list:
    """Create changes at specified positions in the sequences.

    Args:
        sequences (list): List of sequences to modify.
        positions (list): List of positions to modify in each sequence.
        change_frequency (int): Frequency of changes (e.g., every Nth sequence).
        same_character_freq (int): Frequency of keeping the same character.
        valid_chars (list): List of valid characters to use for changes.

    Returns:
        list: Modified list of sequences.
    """
    
    #these are static characters to change things to
    pos1_change = random.choice(valid_chars)
    pos2_change = random.choice(valid_chars)

    #ensure there is change
    if pos1_change == "A":
        pos1_change = "L"
        
    if pos2_change == "A":
        pos2_change = "N"
    
    for i in range(len(sequence1)):
        change = random.randint(1, 100)
        #this will change the sequences together
        if change >= change_frequency:
            same_change = random.randint(1, 100)
            #change both positions to the new characters but "consistent characters"
            if same_change >= same_character_freq:
                sequence1[i] = sequence1[i][:positions[0]] + pos1_change + sequence1[i][positions[0]+1:]
                sequence1[i] = sequence1[i][:positions[1]] + pos2_change + sequence1[i][positions[1]+1:]
            else:
                #create random changes at the positions 
                pos1_change = random.choice(valid_chars)
                pos2_change = random.choice(valid_chars)
                sequence1[i] = sequence1[i][:positions[0]] + pos1_change + sequence1[i][positions[0]+1:]
                sequence1[i] = sequence1[i][:positions[1]] + pos2_change + sequence1[i][positions[1]+1:]
                

    return sequence1

def write_fasta(sequences: list, output_file: str):
    """Write sequences to a FASTA file.

    Args:
        sequences (list): List of sequences to write.
        output_file (str): Path to the output FASTA file.
    """
    with open(output_file, 'w') as f:
        for i, seq in enumerate(sequences):
            f.write(f">seq{i+1}\n")
            f.write(f"{seq}\n")
    return 1

def main():
    #reproducibility
    random.seed(1994)
    
    args = get_args()

    os.makedirs(args.output, exist_ok=True)

    sequence1 = generate_sequences(args.num_sequences, seq_length=20)

    #generate paired sequences if not single
    if not args.single:
        sequence2 = generate_sequences(args.num_sequences, seq_length=20)

    #I am now going to create changes at specific positions with different 
    #frequencies to test the MIT calculations

    #positions 2 and 5 (index 1 and 4)
    #this will be the highest mutual information positions
    #it will have 50% changes with 100% same character frequency (both positions always change together)
    if args.single:
        sequence1 = create_changes_single(
            sequence1, 
            positions=[1, 4], 
            change_frequency=50, 
            same_character_freq=100,
            valid_chars=AMINO_ACIDS
        )
    else:
        #position 1 in the first sequence and 4 in the second sequence
        sequence1, sequence2 = create_changes_paired(
            sequence1, 
            sequence2,
            positions=[1, 4], 
            change_frequency=50, 
            same_character_freq=100,
            valid_chars=AMINO_ACIDS
        )
    
    #positions 7 and 10 (index 6 and 9)
    #it will have 50% changes with 10% same character frequency (both positions change but rarely to the same character each time)
    if args.single:
        sequence1 = create_changes_single(
            sequence1, 
            positions=[6, 9], 
            change_frequency=50, 
            same_character_freq=10,
            valid_chars=AMINO_ACIDS
        )
    else:
        #first position will change in the first sequence and second position in the second sequence
        sequence1, sequence2 = create_changes_paired(
            sequence1, 
            sequence2,
            positions=[6, 9], 
            change_frequency=50, 
            same_character_freq=10,
            valid_chars=AMINO_ACIDS
        )
    
    #positions 12 and 14 (index 11 and 13)
    #it will have 10% changes with 100% same character frequency (both positions change but rarely to the same character each time)
    if args.single:
        sequence1 = create_changes_single(
            sequence1, 
            positions=[11, 13], 
            change_frequency=10, 
            same_character_freq=100,
            valid_chars=AMINO_ACIDS
        )
    else:
        #first position will change in the first sequence and second position in the second sequence
        sequence1, sequence2 = create_changes_paired(
            sequence1, 
            sequence2,
            positions=[11, 13], 
            change_frequency=10, 
            same_character_freq=100,
            valid_chars=AMINO_ACIDS
        )

    #positions 16 and 19 (index 15 and 18)
    #it will have 10% changes with 30% same character frequency (both positions change but rarely to the same character each time)
    if args.single:
        sequence1 = create_changes_single(
            sequence1, 
            positions=[15, 18], 
            change_frequency=10, 
            same_character_freq=30,
            valid_chars=AMINO_ACIDS
        )
    else:
        #first position will change in the first sequence and second position in the second sequence
        sequence1, sequence2 = create_changes_paired(
            sequence1, 
            sequence2,
            positions=[15, 18], 
            change_frequency=10, 
            same_character_freq=30,
            valid_chars=AMINO_ACIDS
        )

    if args.single:
        write_fasta(sequence1, f"{args.output}/single_MSA.fasta")
    else:
        write_fasta(sequence1, f"{args.output}/paired_MSA_1.fasta")
        write_fasta(sequence2, f"{args.output}/paired_MSA_2.fasta")

    return 0

if __name__ == "__main__":
    main()