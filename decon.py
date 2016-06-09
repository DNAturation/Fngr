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

    parser.add_argument('-i', '--in-place', action='store_true',
                        help='Damn the torpedoes! Overwrite the input \
                             genome with the decontaminated version, \
                             and suppress output of the .contamination file. \
                             The action cannot be reversed [off]')

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

def create_outdir(outpath, name):

    subdir = os.path.join(outpath, name)

    if not os.access(subdir, os.F_OK):
        os.mkdir(subdir)

    return subdir

def write_output(contigs, name, outpath):

    with open(os.path.join(outpath, name),'w') as o:
        seqs = (SeqRecord(contigs[rec], rec, description='') for rec in contigs)
        SeqIO.write(sorted(seqs, key=lambda x: x.id), o, 'fasta')

def main():

    args = arguments()

    name = os.path.splitext(os.path.basename(args.genome))[0]
    contaminants = locate_contaminant_names(args.fngr, args.threshold)
    good, bad = separate_contigs(args.genome, contaminants)

    if args.in_place:
        write_output(good, name + '.fasta', os.path.dirname(args.genome))
    else:
        outdir = create_outdir(args.output, name)
        write_output(good, name + '.fasta', outdir)
        write_output(bad, name + '.contamination', outdir)

if __name__ == '__main__':
    main()
