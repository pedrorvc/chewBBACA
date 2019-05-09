#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 10:16:33 2019

@author: pcerqueira
"""


from contextlib import nullcontext as does_not_raise
import pytest

def reverseComplement(strDNA):
    strDNA = strDNA.upper()
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    strDNArevC = ''
    for l in strDNA:
        strDNArevC += basecomplement[l]

    return strDNArevC[::-1]




@pytest.mark.parametrize(
        "test_input, expectation",
        [(5, pytest.raises(AttributeError)), # Tests integer input
         ([], pytest.raises(AttributeError)), # Tests list input
         ("Aello World", pytest.raises(KeyError)), # Tests string input
         ("AACTGGCATGCTATGCAT", does_not_raise()) # Tests valid input to confirm that no exception is raised
         ])
def test_reverseComplement_invalid_inputs(test_input, expectation):
    """Tests the behaviour of the reverseComplement function with unexpected inputs"""
    with expectation:
        reverseComplement(test_input)




@pytest.mark.parametrize(
        "test_input, expected",
        [("AACTGGCATgctatgcat", "ATGCATAGCATGCCAGTT"), # Test lower case sequence 
         (pytest.param("AAAAAACCCCTTTTGGCTTATCGX", "", marks=pytest.mark.xfail(raises=KeyError))), # Test ambiguous bases
         (pytest.param("AACTGGCATgctatgcat", "TTGACCGTACGATACGTA", marks=pytest.mark.xfail(raises=AssertionError))) # Test output orientation. Ensure the output sequence is reversed
         ])
def test_reverseComplement_invalid_sequences(test_input, expected):
    """Tests the behaviour of the reverseComplement function with different sequence inputs"""
    assert reverseComplement(test_input) == expected


