#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  8 11:04:46 2019

@author: pcerqueira
"""

from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

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



def test_empty_translate():
    with pytest.raises(TypeError):
        translateSeq()
 
def test_translate_int():
    with pytest.raises(AttributeError):
        translateSeq(9)

def test_translate_string():
    with pytest.raises(KeyError):
        translateSeq("Aello World")
        
def test_translate_reverse_sense_strand():
    assert str(translateSeq("AATAACTGGGTA")[0]) == "MGQ"

def test_amb():
    with pytest.raises(TranslationError):
        translateSeq("ATGAAACCCCTTTTGGCTTATCGTAA")


def test_empty_list():
    with pytest.raises(AttributeError):
        translateSeq([])


def test_lower_case():
    assert str(translateSeq("atgggtcaataa")[0]) == "MGQ"


def test_translate_sense_strand():
    assert str(translateSeq("ATGGGTCAATAA")[0]) == "MGQ"


def test_translate_antisense_strand():
    assert str(translateSeq("TTATTGACCCAT")[0]) == "MGQ"


def test_translate_reverse_antisense_strand():
    assert str(translateSeq("TACCCAGTTATT")[0]) == "MGQ"


def test_translate_atg_mid():
    assert str(translateSeq("ATGGGTCAATTCATGAAGTAA")[0]) == "MGQFMK"


def test_translate_stop_mid():
    with pytest.raises(TranslationError):
        translateSeq("ATGGGTTAACAATAA")

def test_iupac_codons():
    with pytest.raises(KeyError):
        translateSeq("ATGGGTCAATTCATWSRGAAGTAA")

