#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Blast import NCBIXML
from collections import defaultdict
from decimal import Decimal
from io import StringIO
import json
import subprocess

class Reporter(object):

    def __init__(self, handle, foreign_indices, phylogeny, fragment_size,
                 nt_path, top = 1):

        self.genome = self._load_genome(handle)
        self.foreign_indices = foreign_indices
        self.nt_path = nt_path
        self.top = top
        self.phylogeny = phylogeny
        self.fragment_size = fragment_size

    def _load_genome(self, handle):

        g = {contig.id: str(contig.seq)
             for contig in SeqIO.parse(StringIO(handle), 'fasta')}

        return g

    def _subseq(self, contig, start_stop):

        start, stop = start_stop
        return self.genome[contig][start:stop]

    def _blast(self, seq):

        if self.nt_path:

            blastn = ('blastn', '-db', self.nt_path, '-outfmt', '5')
            hits = []

            blastn_out = subprocess.check_output(blastn, input = seq,
                                                 universal_newlines = True)

            result = NCBIXML.read(StringIO(blastn_out))

            for aln in result.alignments:
                if len(hits) < top:
                    hits.append(aln.title)
                else:
                    break

            out = hits

        else:
            out = []

        return out

    def _parse_foreign_phylo(self, contig, start, stop):

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

        return float(Decimal(str(a)) / Decimal(str(b)))

    def _gc_content(self, sequence):

        gcs = sum(1 for x in sequence.upper() if x in "GC")
        l = len(sequence)

        return self._ddivide(gcs, l)

    def _create_json(self):

        def result_json(contig, index_pair):

            start, stop = index_pair
            seq = self.genome[contig][start:stop]

            length_ratio = self._ddivide(len(seq), len(self.genome[contig]))
            phylo = self._parse_foreign_phylo(contig, start, stop)

            out = {'index': {'start': start,
                             'stop': stop},

                   'length': {'query': len(seq),
                              'contig': len(self.genome[contig]),
                              'ratio': length_ratio},

                   'gc': {'query': self._gc_content(seq),
                          'contig': self._gc_content(self.genome[contig])},

                   'blast_hits': list(enumerate(self._blast(seq), 1)),

                   'read_classification': phylo,

                   'sequence': seq}

            return out

        r = defaultdict(list)

        for contig in self.foreign_indices:
            for index_pair in self.foreign_indices[contig]:

                r[contig].append(result_json(contig, index_pair))

        return json.dumps(r, separators = (', ', ': '))

    def report(self):

        output = self._create_json()
        print(output)  # report is sent to stdout
