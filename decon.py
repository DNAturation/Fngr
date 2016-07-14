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

    parser.add_argument('-m', '--remove-mobile', action='store_true',
                        help='Remove contigs composed solely  of mobile \
                              elements like plasmids or phages as though \
                              they were contamination [off]')

    return parser.parse_args()

def locate_contaminant_names(fngr, threshold, rm_mob):

    def above_threshold(data, contig):

        value = float(data[contig][0]['length']['ratio'])

        return value >= threshold

    def handle_mobile(data, contig):

        # these are deliberately truncated to increase sensitivity
        mobs = ('plasmid', 'phage', 'transpos', 'extra-chromo',
                'extra chromo', 'mobil')

        if not rm_mob:

            blast_results = data[contig][0]['blast_hits']

            blast_names = [b['align_title'] for b in blast_results]

            kraken = list(data[contig][0]['read_classification'].keys())

            hits = [o.lower() for o in kraken + blast_names]

            return not any(m in hits for m in mobs)

        else:
            True

    def is_removable(data, contig, rm_mob):

        thresh = above_threshold(data, contig)
        mob = handle_mobile(data, contig, rm_mob)

        return thresh and mob

    with open(fngr) as f:
        data = json.load(f)

    return [contig for contig in data if is_removable(data, contig, rm_mob)]

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
        seqs = (SeqRecord(contigs[s], s, description='') for s in contigs)
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
