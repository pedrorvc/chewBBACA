#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 10:16:33 2019

@author: pcerqueira
"""

from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

from contextlib import nullcontext as does_not_raise
import pytest



def reverseComplement(strDNA):
    strDNA = strDNA.upper()
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




test_seqs = [
        ("atgggtcaataa", "MGQ"), # Test lower case sequence
        ("ATGGGTCAATAA", "MGQ"), # Test sense strand sequence
        ("TTATTGACCCAT", "MGQ"), # Test antisense strand sequence
        (pytest.param("AATAACTGGGTA", "MGQ", marks=pytest.mark.xfail(raises=TranslationError))), # Test reverse sense strand sequence
        ("TACCCAGTTATT", "MGQ"), # Test reverse antisense strand sequence
        ("ATGGGTCAATTCATGAAGTAA", "MGQFMK")] # Test start codon in the middle of the sequence


@pytest.mark.parametrize(
        "test_input, expectation",
        [("Aello World", pytest.raises(KeyError)), # Test string input
         ([], pytest.raises(AttributeError)), # Test empty list input
         (5, pytest.raises(AttributeError)), # Test integer input
         ("ATGGGTTAACAATAA", pytest.raises(TranslationError)), # Test extra stop codon in sequence 
         ("NTGAAACCCCTTTTGGCTTATCGTAN", pytest.raises(KeyError)), # Test ambiguous bases
         ("ATGGGTCAATTCATWSRGAAGTAA", pytest.raises(KeyError)), # Test IUPAC codons 
         ("ATGGGTCAATAA", does_not_raise()) # Test valid input to confirm that no exception is raised
         ])
def test_translate_invalid_inputs(test_input, expectation):
    """Tests the behaviour of the translateSeq function with unexpected inputs"""
    with expectation:
        translateSeq(test_input)

    

@pytest.mark.parametrize("seq_input, expected", test_seqs)
def test_translate_incorrect_sequences(seq_input, expected):
    """Tests the behaviour of the translateSeq function with different sequence inputs"""
    assert str(translateSeq(seq_input)[0]) == expected
 

    






















