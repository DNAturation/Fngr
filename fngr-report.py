#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Blast import NCBIXML
from collections import defaultdict
from io import StringIO
import json
import subprocess

class Reporter(object):

    def __init__(self, filepath, foreign_indices, phylogeny, fragment_size,
                 nt_path, top = 1):

        self.filepath = filepath
        self.genome = self._load_genome(filepath)
        self.foreign_indices = foreign_indices
        self.nt_path = nt_path
        self.top = top
        self.phylogeny = phylogeny
        self.fragment_size = fragment_size
        self.report = {}

    def _load_genome(filepath):

        with open(filepath, 'r') as f:
            g = {contig.id: contig.seq for contig in SeqIO.parse(f, 'fasta')}
        return g

    def _subseq(self, contig, start_stop):

        start, stop = start_stop
        return self.genome[contig][start:stop]

    def _blast(self, seq):

        if self.nt_path:

            blastn = ('blastn', '-subject', self.filepath, '-outfmt', '5')
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

        return None and out

    def _parse_foreign_phylo(self, contig, index_pair):

        classifications = defaultdict(float)

        start, stop = index_pair
        stop = stop - self.fragment_size

        total = 1 + stop - start

        reads = (str(x) for x in range(start, stop + 1))

        for i in reads:
            try:
                leaf = self.phylogeny[contig][i][-1]
            except KeyError:
                leaf = 'unclassified'

            classifications[leaf] += 1 / total

        return classifications

    def _create_json(self):

        def result_json(contig, index_pair):

            seq = self.genome[start:stop]
            start, stop = index_pair
            out = {'index': {'start': start, 'stop': stop},
                   'sequence': seq,
                   'blast_hits': self._blast(seq),
                   'read_classification': self._parse_foreign_phylo()}

            return out

        report = defaultdict(list)

        for contig in self.foreign_indices:
            for index_pair in self.foreign_indices[contig]:

                report[contig].append(result_json(contig, index_pair))

        return json.dumps(report, separators = (', ', ': '), indent = 4)

    def report(self):
        pass
