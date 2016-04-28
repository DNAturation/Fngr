#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict
from itertools import groupby
from multiprocessing import cpu_count
import argparse
import os
import subprocess
import sys

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--organism', type = str.lower, required = True,
                        help = 'Most precise taxonomic description matching \
                                the target organism')

    parser.add_argument('--kraken-database', required = True,
                        help = 'Path to Kraken DB')

    parser.add_argument('--nt-database', default = None,
                        help = 'Path to NCBI `nt` database')

    parser.add_argument('--threshold', type = int, default = 100,
                        help = 'Number of consecutive pseudoreads not \
                                matching `organism` required identify \
                                the region as being of foreign origin [100]')

    parser.add_argument('--fragment', type = int, default = 250,
                        help = 'Size of pseudoreads to be generated [250]')

    parser.add_argument('--cores', type = int, default = cpu_count(),
                        help = 'Number of CPU cores to use [all]')

    parser.add_argument('assembly', help = 'FASTA formatted file \
                                            or \'-\' to read from stdin')

    return parser.parse_args()

def generate_pseudoreads(filepath, fragment_size):

    def read_genome(filepath):

        o = sys.stdin if filepath == '-' else open(filepath, 'r')

        with o as f:
            for contig in SeqIO.parse(f, 'fasta'):
                yield contig.id, contig.seq

    def fragment_contig(seq, fragment_size):

        for start in range(len(seq) - fragment_size + 1):
            end = start + fragment_size
            yield start, seq[start:end]

    pseudoreads = {}

    for contig_name, contig_seq in read_genome(filepath):
        pseudoreads[contig_name] = dict(fragment_contig(contig_seq,
                                                        fragment_size))

    return pseudoreads

def format_query(filepath, fragment_size):

    def format_contig(genome_name, contig_name):
        def f(pair):

            start, seq = pair

            title = '|'.join(['>' + genome_name, contig_name, str(start)])

            return '\n'.join([title, str(seq), ''])
        return f

    out = []

    genome_name = os.path.splitext(os.path.basename(filepath))[0]
    pseudoreads = generate_pseudoreads(filepath, fragment_size)

    pseudoread_counts = {key: len(pseudoreads[key]) for key in pseudoreads}

    for contig_name in pseudoreads:

        format_pseudoread = format_contig(genome_name, contig_name)

        out += [format_pseudoread(read) for
                read in pseudoreads[contig_name].items()]

        return ''.join(out), pseudoread_counts

def classify(queries, cores, db):

    kraken = ('kraken',
              '--threads', str(cores),
              '--db', db,
              '--fasta-input',
              '/dev/fd/0')

    krak_trans = ('kraken-translate', '--db', db)

    kraken_out = subprocess.check_output(kraken, input = queries,
                                         universal_newlines = True)

    translated = subprocess.check_output(krak_trans, input = kraken_out,
                                         universal_newlines = True)

    return translated

def parse_phylogeny(translated):

    calls = defaultdict(dict)

    for line in translated.strip().split('\n'):

        header, phylo = line.strip().split('\t')
        genome, contig, start = header.split('|')

        calls[contig][start] = phylo.lower().split(';')

    return calls

def determine_origin(calls, counts, root):

    origins = {}

    for contig in calls:
        classified = [-1 for i in range(counts[contig])]  # init unclassified

        for s in calls[contig]:

            start = int(s)

            # 0 is possibly foreign, 1 matches root taxonomy
            # unclassified reads are not in calls, so are left -1
            classified[start] = int(root in calls[contig][s])

        origins[contig] = classified

    return origins

def identify_foreign(origins, threshold, fragment_size):

    def process(call, vals):

        seq_len = sum(1 for x in vals) + fragment_size

        flagged = call is not 1 and seq_len >= (threshold + fragment_size - 1)

        return flagged, seq_len

    processed = 0

    foreign_indices = defaultdict(list)

    for contig in origins:
        for call, values in groupby(origins[contig]):

            flagged, seq_len = process(call, values)

            if flagged:

                seq_indices = (processed, processed + seq_len)

                foreign_indices[contig].append(seq_indices)

            processed += seq_len

    return foreign_indices

def main():

    args = arguments()

    query, pseudoread_counts = format_query(args.assembly, args.fragment)

    classifications = classify(query, args.cores, args.db)

    parsed = parse_phylogeny(classifications)

    origins = determine_origin(parsed, pseudoread_counts, args.organism)

    foreign_indices = identify_foreign(origins, args.threshold, args.fragment)
if __name__ == '__main__':
    main()
