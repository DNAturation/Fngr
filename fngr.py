#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict
from functools import partial
from io import StringIO
from itertools import groupby
from math import sqrt
from multiprocessing import cpu_count
import argparse
import concurrent.futures
import os
import reporter
import subprocess
import sys

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--organism', type = str.lower, required = True,
                        help = 'Least precise taxonomic description matching \
                                the target organism')

    parser.add_argument('--kraken-database', required = True,
                        help = 'Path to Kraken database')

    parser.add_argument('--nt-database', default = None,
                        help = 'Path to NCBI `nt` database')

    parser.add_argument('--threshold', type = int, default = 100,
                        help = 'Number of consecutive pseudoreads not \
                                matching `organism` required identify \
                                the region as being of foreign origin [100]')

    parser.add_argument('--fragment', type = int, default = 250,
                        help = 'Size of pseudoreads to be generated [250]')

    parser.add_argument('--genome-name', default=None,
                        help = 'Override the default genome name \
                                [infer from filename or \'stdin\']')

    parser.add_argument('--cores', type = int, default = cpu_count(),
                        help = 'Number of CPU cores to use [all]')

    parser.add_argument('assembly', help = 'FASTA formatted file \
                                            or \'-\' to read from stdin')

    return parser.parse_args()

def handle_input(filepath, genome_name):
    """Handles file input from either a specified filepath or stdin.

    If `filepath` is a dash (-), stdin is used,
    otherwise `filepath` is taken to be exactly that.

    The variable `name` is used as an identifier for
    pseudoread FASTA in format_query(). It is taken from the filepath basename,
    overridden by `genome_name`, or given a placeholder.
    """

    if genome_name:
        name = genome_name
    elif filepath != '-':
        name = os.path.splitext(os.path.basename(filepath))[0]
    else:
        name = 'stdin'

    f = sys.stdin if filepath == '-' else open(filepath)

    with f as o:
        return o.read(), name


def generate_pseudoreads(handle, fragment_size):
    """Reads in the genome from the `handle` file-like object,
    and cuts into pseudoreads of length `fragment_size`.

    Returns a dictionary of pseudoreads for each contig.
    """

    def read_genome(f):

        for contig in SeqIO.parse(StringIO(f), 'fasta'):
            yield contig.id.replace('|',';'), contig.seq

    def fragment_contig(seq, fragment_size):

        for start in range(len(seq) - fragment_size + 1):
            end = start + fragment_size
            yield start, seq[start:end]

    pseudoreads = {}

    for contig_name, contig_seq in read_genome(handle):
        pseudoreads[contig_name] = dict(fragment_contig(contig_seq,
                                                        fragment_size))
    return pseudoreads

def format_query(handle, genome_name, fragment_size):
    """Returns a multifasta formatted string of pseudoreads
    and a dictionary containing the number of pseudoreads corresponding
    to each contig.

    Example FASTA header for a pseudoread:

    >NCTC11168|contig_0001|42

    The genome name, contig name, and pseudoread start position within the
    contig are given separated by pipe characters.
    """

    def format_contig(genome_name, contig_name):
        def f(pair):

            start, seq = pair

            title = '|'.join(['>' + genome_name, contig_name, str(start)])

            return '\n'.join([title, str(seq), ''])
        return f

    out = []

    pseudoreads = generate_pseudoreads(handle, fragment_size)

    pseudoread_counts = {key: len(pseudoreads[key]) for key in pseudoreads}

    for contig_name in pseudoreads:
        format_pseudoread = format_contig(genome_name, contig_name)

        out += [format_pseudoread(read) for
                read in pseudoreads[contig_name].items()]

    return ''.join(out), pseudoread_counts

def translate(kraken_out_chunk, db):
    """Runs kraken-translate - called in parallel by classify()"""

    krak_trans = ('kraken-translate', '--db', db)

    return subprocess.check_output(krak_trans,
                                   input=kraken_out_chunk,
                                   universal_newlines=True).strip()

def classify(queries, cores, db):
    """Runs kraken and kraken-translate, which are assumed to be in $PATH"""

    def chunk_kraken_out(kraken_out, frac):
        """Generator which chunks the output of kraken
        into sqrt(args.cores) equal-ish sized pieces
        for each kraken-translate process.
        """
        split_kraken_out = kraken_out.splitlines()
        length = len(split_kraken_out)
        chunk = length // frac

        for i in range(0, length, chunk):
            yield '\n'.join(split_kraken_out[i:i + chunk])

    frac = round(sqrt(cores))
    kraken = ('kraken',
              '--threads', str(cores),
              '--db', db,
              '--fasta-input',
              '/dev/fd/0')

    kraken_out = subprocess.check_output(kraken, input = queries,
                                         universal_newlines = True)

    with concurrent.futures.ProcessPoolExecutor(frac) as executor:

        translated = executor.map(partial(translate, db=db),
                                  chunk_kraken_out(kraken_out, frac))

    return '\n'.join(translated)

def parse_phylogeny(translated):
    """Parses the phylogeny generated by kraken-translate in classify().

    Returns a dictionary in which each positively classified pseudoread
    has a list containing its taxonomy ordered left-to-right,
    generic-to-specific.
    """

    calls = defaultdict(dict)

    for line in translated.splitlines():

        header, phylo = line.strip().split('\t')
        genome, contig, start = header.split('|')

        calls[contig][start] = phylo.lower().split(';')

    return calls

def determine_origin(calls, counts, root):
    """Determines whether each pseudoread originated in the host genome
    or is novel to this species.

    A pseudoread is considered to originate in the host genome if
    the most specific kraken classification for that read if is equal to or
    has as an ancestor `root`.

    Returns a dictionary of potential pseudoread origins for each contig.
    For each position:

     1 -- matches `root` or better
     0 -- positive match to an organism not sharing `root` with host
    -1 -- unclassified / no positive match in the kraken database
    """

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
    """Identifies genomic regions which may be of foreign origin.

    A region is putatively of foreign origin if `threshold` consecutive
    pseudoreads were not postively identified as belonging to the host species
    by determine_origin().

    Returns a dictionary of start and stop positions for each region in each
    contig possessing potentially foreign DNA.
    """

    def process(call, vals):
        """Returns whether a region is foreign, as well as its length"""

        seq_len = sum(1 for x in vals) + fragment_size - 1

        flagged = call is not 1 and seq_len >= (threshold + fragment_size - 1)

        return flagged, seq_len

    foreign_indices = defaultdict(list)

    for contig in origins:

        processed = 0

        for call, values in groupby(origins[contig]):
            flagged, seq_len = process(call, values)

            if flagged:

                seq_indices = (processed, processed + seq_len)

                foreign_indices[contig].append(seq_indices)

            processed += seq_len - fragment_size + 1

    return foreign_indices

def main():

    args = arguments()

    handle, genome_name = handle_input(args.assembly, args.genome_name)

    query, pseudoread_counts = format_query(handle, genome_name, args.fragment)

    classifications = classify(query, args.cores, args.kraken_database)

    phylogeny = parse_phylogeny(classifications)

    origins = determine_origin(phylogeny, pseudoread_counts, args.organism)

    foreign_indices = identify_foreign(origins, args.threshold, args.fragment)

    report = reporter.Reporter(handle, foreign_indices, phylogeny,
                               args.fragment, args.nt_database, args.cores)

    report.report()
if __name__ == '__main__':
    main()
