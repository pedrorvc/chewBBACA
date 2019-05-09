#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 17:07:23 2019

@author: pcerqueira
"""

from Bio import SeqIO
from Bio.Alphabet import generic_dna

import pytest

def is_fasta(filename):
    """ Checks if a file is a FASTA file."""
    
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        
        # returns True if FASTA file, False otherwise
        return any(fasta)
    

is_fasta("/home/pcerqueira/Lab_Software/testing/tests/invalid.fasta")

is_fasta("/home/pcerqueira/Lab_Software/testing/tests/empty.fasta")

is_fasta("/home/pcerqueira/Lab_Software/testing/tests/existing.fasta")
