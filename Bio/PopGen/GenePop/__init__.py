# Copyright 2007 by Tiago Antao.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to work with GenePop.

See http://wbiomed.curtin.edu.au/genepop/ , the format is documented
here: http://wbiomed.curtin.edu.au/genepop/help_input.html .

Classes:
Record           Holds GenePop data.
RecordParser     Parses a GenePop record (file) into a Record object.

_Scanner         Scans a GenePop record.
_RecordConsumer  Consumes GenePop data to a Record object.

Inspired on MedLine Code.

"""
from types import *

from Bio import File
from Bio.ParserSupport import *


class Record:
    """Holds information from a GenePop record.

    Members:
    marker_len         The marker lenght (2 or 3 digit code per allele).    
    
    comment_line       Comment line.

    loci_list          List of loci names.
    
    populations        List of population data.
    
    populations has one element per population. Each element is itself
    a list of individuals, each individual is a pair composed by individual
    name and a list of alleles (2 per marker): Example
    [
        [
            ('Ind1', [(1,2),    (3,3), (200,201)],
            ('Ind2', [(2,None), (3,3), (None,None)],
        ],
        [
            ('Other1', [(1,1),  (4,3), (200,200)],
        ]
    ]

    
    """
    def __init__(self):
        self.marker_len      = 0
        self.comment_line    = ""
        self.loci_list       = []
        self.populations     = []

    def __str__(self):
        rep  = [self.comment_line + '\n']
        rep.append('\n'.join(self.loci_list) + '\n')
        for pop in self.populations:
            rep.append('Pop\n')
            for indiv in pop:
                name, markers = indiv
                rep.append(name)
                rep.append(',')
                for marker in markers:
                    rep.append(' ')
                    for al in marker:
                        if al == None:
                            al = '0'
                        aStr = str(al)
                        while len(aStr)<self.marker_len:
                            aStr = "".join(['0', aStr])
                        rep.append(aStr)
                rep.append('\n')
        return "".join(rep)
    

class RecordParser(AbstractParser):
    """Parses GenePop data into a Record object.

    """
    def __init__(self):
        self._scanner = _Scanner()
        self._consumer = _RecordConsumer()

    def parse(self, handle):
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

class _Scanner:
    """Scans a GenePop record.
    
    There is only one record per file.
    
    """

    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in a GenePop unit record for scanning.  handle is a file-like
        object that contains a Genepop record.  consumer is a
        Consumer object that will receive events as the report is scanned.

        """
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)


        consumer.start_record()
        
        comment_line = uhandle.readline().rstrip()
        consumer.comment(comment_line)
        
        #We can now have one loci per line or all loci in a single line
        #seperated by either space or comma+space...
        #We will remove all commas on loci... that should not be a problem
        sample_loci_line = uhandle.readline().rstrip().replace(',', '')
        all_loci = sample_loci_line.split(' ')
        if len(all_loci)>1: #This is all loci in one line
            for locus in all_loci:
                consumer.loci_name(locus)
        else:
            consumer.loci_name(sample_loci_line)
        next_line = uhandle.readline().rstrip()
        while next_line.upper()<>'POP':
            consumer.loci_name(next_line)
            next_line = uhandle.readline().rstrip()
        consumer.start_pop()
        first_individual = True
        line = uhandle.readline().rstrip()
        while line<>'':
            if line.upper()=='POP':
                consumer.start_pop()
            else:
                (indiv_name, marker_line) = line.split(',')
                markers = marker_line.replace('\t', ' ').split(' ')
                for i in range(len(markers), 0, -1):
                    if markers[i-1] == '':
                        del(markers[i-1])
                if first_individual:
                    first_individual = False
                    if len(markers[0]) == 4: #2 digits per allele
                        marker_len = 2
                    else:
                        marker_len = 3
                    consumer.marker_len(marker_len)
                allele_list = []
                for marker in markers:
                    allele_list.append((
                        int(marker[0:marker_len]),
                        int(marker[marker_len:])
                        ))
                consumer.individual(indiv_name, allele_list)
            line = uhandle.readline().rstrip()
        consumer.end_record()

class _RecordConsumer(AbstractConsumer):
    """Consumer that converts a GenePop record to a Record object.

    Members:
    data    Record with GenePop data.

    """
    def __init__(self):
        self.data = None

    def start_record(self):
        self.data = Record()

    def end_record(self):
        pops = self.data.populations
        loci = self.data.loci_list
        for pop_i in range(len(pops)):
            for indiv_i in range(len(pops[pop_i])):
                for mk_i in range(len(loci)):
                    mk_orig = pops[pop_i][indiv_i][1][mk_i]
                    mk_real = []
                    for al in mk_orig:
                        if al == 0:
                            mk_real.append(None)
                        else:
                            mk_real.append(al)
                    pops[pop_i][indiv_i][1][mk_i] = tuple(mk_real)

    def comment(self, comment_line):
        self.data.comment_line = comment_line

    def loci_name(self, locus):
        self.data.loci_list.append(locus)

    def marker_len(self, marker_len):
        self.data.marker_len = marker_len

    def start_pop(self):
        self.current_pop = []
        self.data.populations.append(self.current_pop)

    def individual(self, indiv_name, allele_list):
        self.current_pop.append((indiv_name, allele_list))
    
