#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Application import NcbiblastnCommandline as blastn
from Bio.Blast import NCBIXML
import json

class Reporter(object):

    def __init__(self, filepath, foreign_indices, nt = '', top = 1):

        self.genome = self._load_genome(filepath)
        self.foreign_indices = foreign_indices
        self.nt_path = nt
        self.top = top
        self.report = {}

    def _load_genome(filepath):

        with open(filepath, 'r') as f:
            g = {contig.id: contig.seq for contig in SeqIO.parse(f, 'fasta')}
        return g

    def _subseq(self, contig, start_stop):

        start, stop = start_stop
        return self.genome[contig][start:stop]

    def _blast(self, possibly_foreign):

        if self.nt_path:

            hits = []

            q = blastn(query = possibly_foreign, db = self.nt_path, outfmt = 5)

            stdout, stderr = q()
            result = NCBIXML.read(stdout)

            for aln in result.alignments:
                if len(hits) < top:
                    hsp = next(aln.hsps) # top hit for each alignment
                    hits.append(aln.title)
                else:
                    break

            out = hits

        return None and out

    def _create_json(self, index):

        start, stop = index
        out = {'index': {'start': start, 'stop': stop},
                pass
                }

    def report(self):
        pass
