#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 15:56:48 2019

@author: pcerqueira
"""
import os

import chewBBACA_CreateSchema_nonPairwise_TEST as rfm
import pickle
from pickle import UnpicklingError
from Bio.Data.CodonTable import TranslationError

from contextlib import nullcontext as does_not_raise
import pytest


os.chdir("/home/pcerqueira/Lab_Software/testing/tests")

#### TEST genome_id ################################################

@pytest.mark.parametrize(
        "test_input, expectation",
        [(pytest.param([], "", marks=pytest.mark.xfail)),       # Test empty list and string input
         (pytest.param("", "", marks=pytest.mark.xfail)),       # Test empty string input
         (5, pytest.raises(AttributeError)),                    # Test integer input
         ("/home/pcerqueira/Lab_Software/testing/tests/GCA_000007265.1_ASM726v1_genomic.fna", does_not_raise()) # Test valid input to confirm that no exception is raised
         ])
def test_genome_id_invalid_inputs(test_input, expectation):
        """Tests the behaviour of the genome_id function with unexpected inputs"""
        with expectation:
                rfm.genome_id(test_input)


@pytest.mark.parametrize(
        "file_input, expected",
        [("unexisting.fasta", "unexisting"),    # Test with an unexisting fasta file input.
         ("existing.fasta", "existing"),        # Test with an existing fasta file input. The file has no content.
         ("invalid.fasta", "invalid"),          # Test with an invalid fasta file input. The file has an ID and an invalid DNA/protein sequence.
         ("empty.fasta", "empty"),              # Test with an empty fasta file input. The file has an ID but no sequence.
         ("/home/pcerqueira/Lab_Software/testing/tests/GCA_000007265.1_ASM726v1_genomic.fna", "GCA_000007265.1_ASM726v1_genomic")]) # Test with a valid input.
def test_genome_id_invalid_file_inputs(file_input, expected):
        """Tests the behaviour of the genome_id function with unexpected file inputs"""
        assert rfm.genome_id(file_input) == expected


##################################################################
    
    
#### TEST import_contigs ########################################
    
expected_output = {"GCA-000007265-protein1_1":"ATGGTACAATATAACAATAATTATCCACAAGACAATAAGGAAGAAGCTATGACGGAAAACGAACAACTATTTTGGAATAGAGTACTAGAGCTATCTCGTTCTCAAATAGCACCAGCAGCTTATGAATTTTTTGTTCTAGAGGCTAGACTCCTCAAAATTGAACATCAAACTGCAGTTATTACTTTAGATAACATTGAAATGAAAAAGCTATTCTGGGAACAAAATTTGGGGCCTGTTATCCTAACAGCTGGTTTTGAAATTTTCAATGCTGAAATTACAGCTAACTATGTCTCAAACGATTTACATTTACAAGAAACTAGTTTTTCTAACTACCAGCAATCTAGCAATGAAGTAAATACTTTACCAATTAGAAAAATCGACTCTAACCTTAAAGAGAAATATACTTTTGCTAATTTTGTTCAAGGAGATGAAAATAGATGGGCTGTTTCAGCATCAATTGCTGTAGCTGATAGTCCTGGCACGACTTATAATCCTCTATTTATCTGGGGGGGACCTGGTCTAGGAAAGACGCATCTACTAAATGCTATTGGAAATCAAGTCTTAAGAGATAATCCAAACGCGAGGGTTTTATACATCACTGCTGAGAATTTTATTAATGAATTTGTCAGTCATATTCGTTTAGATTCGATGGAAGAATTAAAAGAAAAGTTTCGCAACTTAGACTTACTCCTGATTGATGATATTCAGTCGCTTGCTAAGAAAACCTTAGGGGGGACCCAAGAGGAGTTCTTCAATACTTTCAATGCTTTACATACAAACGATAAACAAATCGTATTGACCAGTGACCGAAATCCAAATCAATTAAATGATCTAGAAGAACGTCTAGTCACGCGCTTTAGTTGGGGACTCCCAGTAAATATCACACCACCTGATTTTGAAACACGAGTTGCTATTTTAACCAATAAAATTCAAGAATATCCTTATGATTTTCCTCAAGATACCATTGAATACTTAGCAGGAGAATTTGATTCCAACGTACGTGAATTAGAAGGAGCCTTGAAAAATATTAGTCTAGTTGCTGACTTTAAGCATGCTAAAACTATTACAGTAGATATAGCTGCAGAAGCTATCAGAGCACGTAAAAATGATGGTCCTATTGTTACTGTCATTCCTATAGAAGAAATTCAAATACAAGTTGGTAAATTCTATGGCGTAACTGTAAAAGAGATAAAAGCAACTAAAAGAACACAAGATATTGTCCTTGCAAGACAGGTAGCCATGTACTTAGCTCGTGAGATGACAGATAACAGTCTCCCAAAAATAGGTAAAGAATTTGGGGGACGAGATCACTCAACTGTTCTCCACGCTTATAATAAAATAAAAAATATGGTTGCTCAAGATGACAACTTACGAATTGAGATAGAAACTATCAAAAATAAAATCAGATAG"}

@pytest.mark.parametrize(
        "test_input, expectation",
        [(pytest.param([], "", marks=pytest.mark.xfail)),       # Test empty list and string input
         (pytest.param("", "", marks=pytest.mark.xfail)),       # Test empty string input
         (5, pytest.raises(TypeError))                          # Test integer input
         ])
def test_import_contigs_invalid_inputs(test_input, expectation):
        """Tests the behaviour of the import_contigs function with unexpected inputs"""
        with expectation:
                rfm.import_contigs(test_input)    


@pytest.mark.parametrize(
        "file_input, expected",
        [(pytest.param("/home/pcerqueira/Lab_Software/testing/tests/unexisting.fasta", "", marks=pytest.mark.xfail(raises=FileNotFoundError))),  # Test with an unexisting fasta file input
         (pytest.param("/home/pcerqueira/Lab_Software/testing/tests/existing.txt", "", marks=pytest.mark.xfail)),                               # Test with an existing fasta file input. The file has no content.
         (pytest.param("invalid.fasta", {'ContigID': 'REALLIYINAVLIDFASTASEQUENCE'}, marks=pytest.mark.xfail)),                                 # Test with an invalid fasta file input. The file has an ID and an invalid DNA/protein sequence.
         (pytest.param("empty.fasta", {'ContigID': ''}, marks=pytest.mark.xfail)),                                                              # Test with an empty fasta file input. The file has an ID but no sequence.
         ("/home/pcerqueira/Lab_Software/testing/tests/GCA-000007265-protein1.fasta", expected_output)                                          # Test with a valid input.
         ])
def test_import_contigs_invalid_file_inputs(file_input, expected):
        """Tests the behaviour of the import_contigs function with unexpected file inputs"""
        assert rfm.import_contigs(file_input) == expected

##############################################################


#### TEST extract_coding_sequences ########################################

orf_file_path = "GCA_000007265.1_ASM726v1_genomic.fna_ORF.txt"

contigs = rfm.import_contigs("GCA_000007265.1_ASM726v1_genomic.fna")

valid_output = rfm.extract_coding_sequences(orf_file_path, contigs, 0, 0)

@pytest.mark.parametrize(
        "reading_frames, contigs, starting_id, genome_id, expectation",
        [([],[],[],[], pytest.raises(TypeError)),                       # Test empty list input
         ("","","","", pytest.raises(FileNotFoundError)),               # Test empty string input
#         (7,7,7,5, pytest.raises(TypeError)),
         (orf_file_path, "", "", "", pytest.raises(TypeError)),         # Test with vaild file and empty string input
         (orf_file_path, contigs, "", "", pytest.raises(TypeError)),    # Test with valid file, contig dict and empty string input
         ])
def test_extract_coding_sequences_invalid_inputs(reading_frames, contigs, starting_id, genome_id, expectation):
        """Tests the behaviour of the extract_coding_sequences function with unexpected inputs"""
        with expectation:
                rfm.extract_coding_sequences(reading_frames, contigs, starting_id, genome_id)   
        
def test_extract_coding_sequences_valid_input():
        """Tests the behaviour of the extract_coding_sequences function with a valid input"""
        assert rfm.extract_coding_sequences(orf_file_path, contigs, 0, 0) == valid_output


###########################################################################


#### TEST reverse_complement ########################################

@pytest.mark.parametrize(
        "test_input, expectation",
        [(5, pytest.raises(TypeError)),                                                 # Test integer input
         (pytest.param([], "", marks=pytest.mark.xfail(raises=AttributeError))),        # Test empty list input
         ("Aello World", pytest.raises(KeyError)),
         ("AAAAAACCCCTTTTGGCTTATCGX", pytest.raises(KeyError)),                         # Test ambiguous sequence input
         ("AACTGGCATGCTATGCAT", does_not_raise())                                       # Test valid input to confirm that no exception is raised
         ])
def test_reverse_complement_invalid_inputs(test_input, expectation):
        """Tests the behaviour of the reverse_complement function with unexpected inputs"""
        with expectation:
                rfm.reverse_complement(test_input)

@pytest.mark.parametrize(
        "test_seq, expected",
        [("AACTGGCATgctatgcat", "ATGCATAGCATGCCAGTT"),                                                                  # Test lower case sequence 
         ("AACTGGCATgctatgcan", "NTGCATAGCATGCCAGTT"),                                                                  # Test ambiguous base N
         (pytest.param("AACTGGCATgctatgcat", "TTGACCGTACGATACGTA", marks=pytest.mark.xfail(raises=AssertionError))),    # Test output orientation. Ensure the output sequence is reversed
         ])
def test_reverse_complement_invalid_sequences(test_seq, expected):
        """Tests the behaviour of the reverseComplement function with different sequence inputs"""
        assert rfm.reverse_complement(test_seq) == expected

###########################################################################


#### TEST translate #######################################################

test_seqs = [
        ("Aello World", "amb"),                 # Test with string input
        ("NTGAAACCCCTTTTGGCTTATCGTAN", "amb"),  # Test with ambiguous N base
        ("ATGGGTCAATTCATWSRGAAGTAA", "amb"),    # Test with other ambiguous bases
        ("atgggtcaataa", "MGQ"),                # Test with lower case sequence
        ("ATGGGTCAATAA", "MGQ"),                # Test sense strand sequence
        ("TTATTGACCCAT", "MGQ"),                # Test antisense strand sequence
        ("AATAACTGGGTA", "MGQ"),                # Test reverse sense strand sequence
        ("TACCCAGTTATT", "MGQ"),                # Test reverse antisense strand sequence
        ("ATGGGTCAATTCATGAAGTAA", "MGQFMK")]    # Test start codon in the middle of the sequence


@pytest.mark.parametrize(
        "test_input, expectation",
        [([], pytest.raises(AttributeError)),                   # Test with empty list input
         (5, pytest.raises(AttributeError)),                    # Test with integer input
         ("ATGGGTTAACAATAA", pytest.raises(TranslationError)),  # Test extra stop codon in sequence
         ("ATGGGTCAATAA", does_not_raise())                     # Test with valid input
         ])
def test_translate_invalid_inputs(test_input, expectation):
        """Tests the behaviour of the translate function with unexpected inputs"""
        with expectation:
                rfm.translate(test_input)

    
@pytest.mark.parametrize("seq_input, expected", test_seqs)
def test_translate_incorrect_sequences(seq_input, expected):
        """Tests the behaviour of the translate function with different sequence inputs"""
        assert str(rfm.translate(seq_input)[0]) == expected or str(rfm.translate(seq_input)) == 'amb' or str(rfm.translate(seq_input)) == 'untrans'

###########################################################################


#### TEST translate_coding_sequences #######################################################

test_seqs_tcs = [({1:"ATGGGTCAATAA"}, [{1: 'MGQ'}, {1: 'sense'}, {}, {}]),      # Test with valid sequence
             ({1:"ATGNGTCAATAA"}, [{}, {}, {1: 'ATGNGTCAATAA'}, {}]),           # Test with ambiguous sequence
             ({1:"AAAAAAAAAAAA"}, [{}, {}, {}, {1: 'AAAAAAAAAAAA'}])            # Test with untranslatable sequence
             ]

@pytest.mark.parametrize(
        "test_input, expectation",
        [([], pytest.raises(AttributeError)),   # Test with empty list input
         (5, pytest.raises(AttributeError)),    # Test with integer input
         ("",pytest.raises(AttributeError)),    # Test with empty string input
         ])
def test_translate_coding_sequences_invalid_args(test_input, expectation):
        """Tests the behaviour of the translate_coding_sequences function with unexpected inputs"""
        with expectation:
                rfm.translate_coding_sequences(test_input)
    

@pytest.mark.parametrize("seq_input, expected", test_seqs_tcs)
def test_translate_coding_sequences_incorrect_sequences(seq_input, expected):
        """Tests the behaviour of the translate_coding_sequences function with different sequence inputs"""
        assert rfm.translate_coding_sequences(seq_input) == expected

###########################################################################



#### TEST determine_repeated #######################################################

coding_sequences_repeated = {1:"ATGGTCGAACGGTCGTCAAAGTCAGAGTACCCCGGGTACCAACTTACGGAGGATATTGCTTGCAGCTGCTAA",
                             2:"ATGGTCGAACGGTCGTCAAAGTCAGAGTACCCCGGGTACCAACTTACGGAGGATATTGCTTGCAGCTGCTAA",
                             3:"ATGTTCAAAGGCGAAAGAACAAAGGCTTACTGTGCGCAGAGGAACGCCCATTTTGCGGCTGGCGTCTATTGA"}

coding_sequences_no_repeats = {1:"ATGGTCGAACGGTCGTCAAAGTCAGAGTACCCCGGGTACCAACTTACGGAGGATATTGCTTGCAGCTGCTAA",
                               2:"ATGTTCAAAGGCGAAAGAACAAAGGCTTACTGTGCGCAGAGGAACGCCCATTTTGCGGCTGGCGTCTATTGA"}

test_seqs_dr = [({1:"ATGGTCGAACGGTCGTCAAAGTCAGAGTACCCCGGGTACCAACTTACGGAGGATATTGCTTGCAGCTGCTAA",
                  2:"ATGGTCGAACGGTCGTCAAAGTCAGAGTACCCCGGGTACCAACTTACGGAGGATATTGCTTGCAGCTGCTAA",
                  3:"ATGTTCAAAGGCGAAAGAACAAAGGCTTACTGTGCGCAGAGGAACGCCCATTTTGCGGCTGGCGTCTATTGA"}, [[1,2]]),      # Test with valid input
                ({1:"ATGGTCGAACGGTCGTCAAAGTCAGAGTACCCCGGGTACCAACTTACGGAGGATATTGCTTGCAGCTGCTAA",
                 2:"ATGTTCAAAGGCGAAAGAACAAAGGCTTACTGTGCGCAGAGGAACGCCCATTTTGCGGCTGGCGTCTATTGA"}, [])]            # Test with sequences without repeats

@pytest.mark.parametrize(
        "test_input, expectation",
        [([], pytest.raises(AttributeError)),   # Test with empty list input
         (5, pytest.raises(AttributeError)),    # Test with integer input
         ("",pytest.raises(AttributeError)),    # Test with empty string input
         ])
def test_determine_repeated_sequences_invalid_args(test_input, expectation):
        """Tests the behaviour of the determine_repeated function with unexpected inputs"""
        with expectation:
                rfm.determine_repeated(test_input)


@pytest.mark.parametrize("seq_input, expected", test_seqs_dr)
def test_determine_repeated_incorrect_sequences(seq_input, expected):
        """Tests the behaviour of the translate_coding_sequences function with different sequence inputs"""
        assert rfm.determine_repeated(seq_input) == expected


############################################################################################

#### TEST remove_repeated #######################################################

coding_sequences_repeated = {1:"ATGGTCGAACGGTCGTCAAAGTCAGAGTACCCCGGGTACCAACTTACGGAGGATATTGCTTGCAGCTGCTAA",
                             2:"ATGGTCGAACGGTCGTCAAAGTCAGAGTACCCCGGGTACCAACTTACGGAGGATATTGCTTGCAGCTGCTAA",
                             3:"ATGTTCAAAGGCGAAAGAACAAAGGCTTACTGTGCGCAGAGGAACGCCCATTTTGCGGCTGGCGTCTATTGA"}

coding_sequences_no_repeats = {1:"ATGGTCGAACGGTCGTCAAAGTCAGAGTACCCCGGGTACCAACTTACGGAGGATATTGCTTGCAGCTGCTAA",
                               2:"ATGTTCAAAGGCGAAAGAACAAAGGCTTACTGTGCGCAGAGGAACGCCCATTTTGCGGCTGGCGTCTATTGA"}

proteins_repeated = rfm.translate_coding_sequences(coding_sequences_repeated)[0]

proteins_repeated2 = rfm.translate_coding_sequences(coding_sequences_repeated)[0]

proteins_no_repeat = rfm.translate_coding_sequences(coding_sequences_no_repeats)[0]

repeated_list = rfm.determine_repeated(coding_sequences_repeated)




@pytest.mark.parametrize(
        "proteins, repeated_ids, expectation",
        [(5,5, pytest.raises(TypeError))
         ])
def test_remove_repeated_sequences_invalid_args(proteins, repeated_ids, expectation):
        """Tests the behaviour of the remove_repeated function with unexpected inputs"""
        with expectation:
                rfm.remove_repeated(proteins, repeated_ids)
        
#@pytest.mark.parametrize(
#        "proteins_incorrect, repeated_ids_incorrect, expected"
#        [([],[], []),
#         ("","", ""),
#         ({},[], {}),
#         ])
#def test_remove_repeated_incorrect_inputs(proteins_incorrect, repeated_ids_incorrect, expected):
#    assert rfm.remove_repeated(proteins_incorrect, repeated_ids_incorrect) == expected

def test_remove_repeated_valid_input():
        assert rfm.remove_repeated(proteins_repeated, repeated_list) == {1:"MVERSSKSEYPGYQLTEDIACSC", 3:"MFKGERTKAYCAQRNAHFAAGVY"}

def test_dremove_repeated_no_repeat_list_input():
        assert rfm.remove_repeated(proteins_repeated2,[]) == proteins_repeated2



############################################################################################


#### TEST determine_contained #######################################################

contained_dna_seq = "ATGTACTACGCATGCGGGCAGATTATGCAGGTCGAGAGATGCGGGAGAAGTTCTCGACCTTCCCGTGGGACGTGTACCTATCCCCTTATCGAGCATTCCGTTTAA"
contained_protein_seq = rfm.translate(contained_dna_seq)


start_protein_seq = "MYYACGQIM"
end_protein_seq= "TYPLIEHSV"
mid_protein_seq = "MQVERCGRSSRPSRGTCT"
mid2 = "CGQIMQVERCGRSSR"
mid3 = "IMQVERCGRSSRPSR"
mid4 = "RPSRGTCTY"
equal = contained_protein_seq[0]


proteins = {1:contained_protein_seq[0],
            2:start_protein_seq,
            3:end_protein_seq,
            4:mid_protein_seq,
            5:mid2,
            6:mid3,
            7:mid4}

proteins_eq = {1:contained_protein_seq[0],
            2:start_protein_seq,
            3:end_protein_seq,
            4:mid_protein_seq,
            5:mid2,
            6:mid3,
            7:mid4,
            8:equal}

proteins_inv = rfm.invert_dictionary(proteins)
proteins_eq_inv = rfm.invert_dictionary(proteins_eq)

test_seqs_dc = [(proteins_inv, {'MYYACGQIM':2, 'TYPLIEHSV':3, 'RPSRGTCTY':7, 'CGQIMQVERCGRSSR':5, 'IMQVERCGRSSRPSR':6, 'MQVERCGRSSRPSRGTCT':4}),
                (proteins_eq_inv, {'MYYACGQIM':2, 'TYPLIEHSV':3, 'RPSRGTCTY':7, 'CGQIMQVERCGRSSR':5, 'IMQVERCGRSSRPSR':6, 'MQVERCGRSSRPSRGTCT':4})]

@pytest.mark.parametrize(
        "test_input, expectation",
        [([], pytest.raises(AttributeError)),   # Test empty list input 
         (5, pytest.raises(AttributeError)),    # Test integer input
         ("",pytest.raises(AttributeError)),    # Test empty string input
         ])
def test_determine_contained_invalid_args(test_input, expectation):
        """Tests the behaviour of the determine_contained function with unexpected inputs"""
        with expectation:
                rfm.determine_contained(test_input)


@pytest.mark.parametrize("seq_input_dc, expected_dc", test_seqs_dc)
def test_determine_contained_sequences(seq_input_dc, expected_dc):
        """Tests the behaviour of the determine_contained function with different sequence inputs"""
        assert rfm.determine_contained(seq_input_dc) == expected_dc

############################################################################################



