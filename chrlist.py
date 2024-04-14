#!/usr/bin/env python

import os
import argparse
from Bio import SeqIO

def parse_arguments():
    parser = argparse.ArgumentParser(description="A script to generate a list of reference sequence lengths")
    parser.add_argument("fasta", help="Sequence fasta file.")
    parser.add_argument("-o", default='./ref.list', help="Output file")
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    file = SeqIO.parse(args.fasta, 'fasta')
    output_file = args.o
    with open(output_file, 'w') as lfile:
        for seq in file:
            lfile.write(seq.id + "\t" + str(len(seq.seq)) + "\n")

if __name__ == "__main__":
    main()

