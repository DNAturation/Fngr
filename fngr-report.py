#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Application import NcbiblastnCommandline as blastn
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
    
    def _blast(self, query, top):

        if self.nt_path:
            pass
        else:
            return None

    def _create_json(self):
        pass

    def report(self):
        pass
