#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  8 09:57:32 2019

@author: pcerqueira
"""


import pytest

def reverseComplement(strDNA):
    strDNA = strDNA.upper()
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    strDNArevC = ''
    for l in strDNA:
        strDNArevC += basecomplement[l]

    return strDNArevC[::-1]


def test_empty_reverse():
    with pytest.raises(TypeError):
        reverseComplement()


@pytest.mark.xfail(raises=TypeError)
def test_non_str_reverse():
    with pytest.raises(AttributeError):
        reverseComplement(9)
    
def test_amb():
    with pytest.raises(KeyError):
        reverseComplement("AAAAAACCCCTTTTGGCTTATCGX")
    
def test_empty_list():
    with pytest.raises(AttributeError):
        reverseComplement([])
    
def test_lower_case():
    assert reverseComplement("a") == "T"

def test_reverse_result():
    assert reverseComplement("AACTGGCATgctatgcat") == "ATGCATAGCATGCCAGTT"
    
def test_result_orientation():
    with pytest.raises(AssertionError):
        assert reverseComplement("AACTGGCATgctatgcat") == "TTGACCGTACGATACGTA"

