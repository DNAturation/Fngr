#!/usr/bin/env python3

import argparse
from Bio import SeqIO

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-g', '--genome', required=True, metavar='PATH',
                        help='Path to FASTA-formatted genome')

    parser.add_argument('-f', '--fngr', required=True, metavar='PATH',
                        help='Path to Fngr report')

    parser.add_argument('-t', '--threshold', type=float,
                        default=1.0, metavar='NUM',
                        help='Foreign locus to contig length ratio threshold \
                             for contamination. Values >= this are removed \
                             from the import assembly [1.0]')

    parser.add_argument('-o', '--output', default='', metavar='PATH',
                        help='Output folder [\'\']')

    return parser.parse_args()

def main():

    args = arguments()

if __name__ == '__main__':
    main()
