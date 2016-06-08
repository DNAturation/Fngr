#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import json
import os

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

def locate_contaminant_names(fngr, threshold):

    with open(fngr) as f:
        data = json.load(f)

    return [contig for contig in data
            if float(data[contig][0]['length']['ratio']) >= threshold]

def separate_contigs(genome, contaminants):

    good, bad = {}, {}

    with open(genome) as g:
        for rec in SeqIO.parse(g, 'fasta'):

            if rec.id in contaminants:

                bad[rec.id] = rec.seq
                continue

            good[rec.id] = rec.seq

    return good, bad

def write_output(contigs, name, ext, outpath):

    if not os.access(os.path.join(outpath, name), os.F_OK):
        os.mkdir(os.path.join(outpath, name))

    with open(os.path.join(outpath, name, name + ext),'w') as o:
        seqs = (SeqRecord(contigs[rec], rec, description='') for rec in contigs)
        SeqIO.write(sorted(seqs, key=lambda x: x.id), o, 'fasta')

def main():

    args = arguments()

    name = os.path.splitext(os.path.basename(args.genome))[0]
    contaminants = locate_contaminant_names(args.fngr, args.threshold)
    good, bad = separate_contigs(args.genome, contaminants)

    write_output(good, name, '.fasta', args.output)
    write_output(bad, name, '.contamination', args.output)

if __name__ == '__main__':
    main()
