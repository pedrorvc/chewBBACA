#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 10:16:33 2019

@author: pcerqueira
"""

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import os
import argparse
import time
import pickle
import shutil
import multiprocessing
import subprocess

def reverseComplement(strDNA):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    strDNArevC = ''
    for l in strDNA:
        strDNArevC += basecomplement[l]

    return strDNArevC[::-1]


def translateSeq(DNASeq):
    seq = DNASeq
    tableid = 11

    #look for ambiguous bases
    try:
        reverseComplement(seq)
    except:
        raise
    try:
        myseq = Seq(seq)
        protseq = Seq.translate(myseq, table=tableid, cds=True)
    except:
        try:
        	# include code to reverse and see if reverseComplement is a valid CDS...
            seq = reverseComplement(seq)
            myseq = Seq(seq)
            protseq = Seq.translate(myseq, table=tableid, cds=True)

        except:
            try:
                seq = seq[::-1]
                myseq = Seq(seq)
                protseq = Seq.translate(myseq, table=tableid, cds=True)
            except:
                try:
                    seq = seq[::-1]
                    seq = reverseComplement(seq)
                    myseq = Seq(seq)
                    protseq = Seq.translate(myseq, table=tableid, cds=True)
                except Exception as e:
                    #print "translation error"
                    #print e
                    raise

    return protseq, seq



#normal_fasta = "/home/pcerqueira/Lab_Software/testing/small_fasta_files/GCA_000007265.1_ASM726v1_genomic.fna"
#
#
#
#fasta_with_N = "teste_N.fna"
    
normal_fasta = "GTGAGATTAGACAAGTATCTAAAAGTGTCGCGTATCATTAAGCGTCGTCCCGTTGCTAAAGAAGTTGCAGACAAAGGACGTGTCAAGGTGAATGGCGTATTAGCAAAGTCATCAACTGACCTAAAGCTTAATGACCAAGTTGAAATTCGTTTTGGTAATAAATTGTTGACCGTTAAAGTATTAGAAATGAAAGATAGCACTAAGAAAGAAGATGCTATCAAAATGTATGAAATTATTAATGAAACAAGGATAGAAACAGATGAGCAAGCCTAA"

fasta_with_N = "GTGAXATTNGACAAGTATCTAAAAGTGTCGCGTATCATTAAGCGTCGTCCCGTTGCTAAAGAAGTTGCAGACAAAGGACGTGTCAAGGTGAATGGCGTATTAGCAAAGTCATCAACTGACCTAAAGCTTAATGACCAAGTTGAAATTCGTTTTGGTAATAAATTGTTGACCGTTAAAGTATTAGAAATGAAAGATAGCACTAAGAAAGAAGATGCTATCAAAATGTATGAAATTATTAATGAAACAAGGATAGAAACAGATGAGCAAGCCTAA"


rev_normal = reverseComplement(normal_fasta)

rec_N = reverseComplement(fasta_with_N)



normal, _ = translateSeq(normal_fasta)

#amb_char = 'RYWSMKHBVDN'

amb, _ = translateSeq(fasta_with_N)
#any(amb in translateSeq(normal_fasta)) == True

#def test_translate_normal():
#    assert  



#import string
tab = str.maketrans("ACTG", "TGAC")

def reverse_complement_table(seq):
    return seq.translate(tab)[::-1]



import random
import timeit

DNAlength = 50

#randomly generate 100k bases
int_to_basemap = {1: 'A', 2: 'C', 3: 'G', 4: 'T'}
num_strings = 500000
random.seed(775249)
DNAstrings = ["".join([int_to_basemap[random.randint(1,4)] for i in range(DNAlength)])
              for j in range(num_strings)]
#get an idea of what the DNAstrings look like
print(DNAstrings[0:5])


tic=timeit.default_timer()
rcs = [reverseComplement(seq) for seq in DNAstrings]
toc=timeit.default_timer()
baseline = toc - tic
namefunc = {"orignal implementation": reverseComplement,
            "translation table implementation" : reverse_complement_table}


for function_name in namefunc:
    func = namefunc[function_name]
    tic=timeit.default_timer()
    rcs = [func(seq) for seq in DNAstrings]
    toc=timeit.default_timer()
    walltime = toc-tic
    print("""{}
       {:.5f}s total,
       {:.1f} strings per second
       {:.1f}% increase over baseline""".format(
           function_name,
           walltime,
           num_strings/walltime,
           100- ((walltime/baseline)*100)   ))




