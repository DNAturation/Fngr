#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict
from io import StringIO
from multiprocessing import cpu_count
import argparse
import os
import subprocess
import sys

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--organism', required = True,
                        help = 'Most precise taxonomic description matching \
                                the target organism')

    parser.add_argument('--db', required = True, help = 'Path to Kraken DB')

    parser.add_argument('--threshold', type = int, default = 100,
                        help = 'Number of consecutive pseudoreads not \
                                matching `organism` required identify \
                                the region as being of foreign origin [100]')

    parser.add_argument('--fragment', type = int, default = 250,
                        help = 'Size of pseudoreads to be generated [250]')

    parser.add_argument('--cores', type = int, default = cpu_count(),
                        help = 'Number of CPU cores to use [all]')

    parser.add_argument('assembly', help = 'FASTA formatted file to analyze')

    return parser.parse_args()

def generate_pseudoreads(filepath, fragment_size):

    def read_genome(filepath):

        with open(filepath, 'r') as f:
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

def format_query(filepath):

    def format_contig(genome_name, contig_name):
        def f(pair):

            start, seq = pair

            title = '|'.join(['>' + genome_name, contig_name, start])

            return '\n'.join([title, seq, ''])
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
              '--threads', cores,
              '--db', db,
              '--fasta-input', StringIO(queries))

    krak_trans = ('kraken-translate', '--db', db)

    kraken_out = subprocess.check_output(kraken)

    translated = subprocess.check_output(krak_trans, input = kraken_out)

    return translated

def parse_results(translated):

    calls = defaultdict(dict)

    for line in translated:
        header, phylo = line.strip().split('\t')
        genome, contig, start = header.split('|')

        calls[contig][start] = phylo.split(';')

    return calls

def determine_origin(calls, counts, root):

    origin = {}

    for contig in calls:
        classified = [-1 for i in range(counts[contig])]  # init unclassified

        for start in calls[contig]:

            # 0 is possibly foreign, 1 matches root taxonomy
            # unclassified reads are not in calls, so are left -1
            classified[start] = int(root in calls[contig][start])

        origin[contig] = classified

    return origin

def identify_foreign(origin, threshold):
    pass

def reporter():
    pass

def main():

    args = arguments()

if __name__ == '__main__':
    main()
