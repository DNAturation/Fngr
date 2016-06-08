#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Blast import NCBIXML
from collections import defaultdict
from decimal import Decimal
from io import StringIO
import json
import subprocess
import sys

class Reporter(object):

    def __init__(self, handle, foreign_indices, phylogeny, fragment_size,
                 nt_path, cores, top, fast):

        self.genome = self._load_genome(handle)
        self.foreign_indices = foreign_indices
        self.nt_path = nt_path
        self.top = top
        self.phylogeny = phylogeny
        self.fragment_size = fragment_size
        self.cores = cores
        self.fast = fast

    def _load_genome(self, handle):
        """Returns a dictionary representation of the fasta file.
        contig name      => key
        contig sequence  => value
        """

        g = {contig.id.replace('|', ';'): str(contig.seq)
             for contig in SeqIO.parse(StringIO(handle), 'fasta')}

        return g

    def _blast(self, seq):
        """Performs a megablast search of `seq` in a blast database
        (presumably `nt`) if one is provided.

        Returns a list of the sequences names for the `self.top` best hits.
        """

        if self.nt_path:

            blastn = ('blastn',
                      '-db', self.nt_path,
                      '-outfmt', '5',
                      '-num_threads', str(self.cores))
            hits = []

            blastn_out = subprocess.check_output(blastn, input = seq,
                                                 universal_newlines = True)

            result = NCBIXML.read(StringIO(blastn_out))

            for aln in result.alignments:
                if len(hits) < self.top:
                    hits.append({'align_title': aln.title,
                                 'identities': aln.hsps[0].identities,
                                 'length': aln.hsps[0].align_length,
                                 'evalue': aln.hsps[0].expect
                                 })
                else:
                    break

            out = hits

        else:
            out = []

        return out

    def _parse_foreign_phylo(self, contig, start, stop):
        """For an index pair in self.foreign_indices,
        return a dictionary of the proportions of foreign-origin reads.
        """

        classifications = defaultdict(Decimal)

        stop = stop - self.fragment_size

        total = 1 + stop - start
        read_frac = Decimal('1.0') / Decimal(str(total))

        reads = (str(x) for x in range(start, stop + 1))

        for i in reads:
            try:
                leaf = self.phylogeny[contig][i][-1]
            except KeyError:
                leaf = 'unclassified'

            classifications[leaf] += read_frac

        return {key: float(classifications[key]) for key in classifications}

    def _ddivide(self, a, b):
        """Performs Decimal() division and returns a JSON-friendly float()."""
        return float(Decimal(str(a)) / Decimal(str(b)))

    def _gc_content(self, sequence):
        """Calculates the GC content of the given sequence"""
        gcs = sum(1 for x in sequence.upper() if x in "GC")
        l = len(sequence)

        return self._ddivide(gcs, l)

    def _create_json(self):
        """Formats the JSON containing results for every region
        identified in self.foreign_indices.
        """

        def result_json(contig, index_pair):

            start, stop = index_pair
            seq = self.genome[contig][start:stop]

            length_ratio = self._ddivide(len(seq), len(self.genome[contig]))
            phylo = self._parse_foreign_phylo(contig, start, stop)

            will_blast = (not self.fast) or length_ratio < 1.0

            out = {'index': {'start': start,
                             'stop': stop},

                   'length': {'query': len(seq),
                              'contig': len(self.genome[contig]),
                              'ratio': length_ratio},

                   'gc': {'query': self._gc_content(seq),
                          'contig': self._gc_content(self.genome[contig])},

                   'blast_hits': self._blast(seq) if will_blast else None,

                   'read_classification': phylo,

                   'sequence': seq}

            return out

        r = defaultdict(list)

        for contig in self.foreign_indices:
            for index_pair in self.foreign_indices[contig]:

                r[contig].append(result_json(contig, index_pair))

        return r

    def report(self):
        """Orders the creation of JSON-formatted results,
        and sends them to stdout.
        """

        output = self._create_json()
        json.dump(output, sys.stdout,
                  indent=4, separators=(', ', ': '), sort_keys=True)
