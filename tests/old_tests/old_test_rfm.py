#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 08 13:53:42 2019

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

xfail = pytest.mark.xfail

#### TEST genome_id ################################################

def test_genome_id_no_args():
   with pytest.raises(TypeError):
       rfm.genome_id()

@xfail
def test_genome_id_empty_list_input():
   rfm.genome_id([])


# Add check for file existence
#@xfail
def test_genome_id_unexisting_fasta_input():
   assert rfm.genome_id("unexisting.fasta") == "unexisting"

# Add check for file contents ??
#@xfail    
def test_genome_id_existing_fasta_input():
   assert rfm.genome_id("existing.fasta") == "existing"
   
#@xfail    
def test_genome_id_invalid_fasta_input():
   assert rfm.genome_id("invalid.fasta") == "invalid"
   
#@xfail    
def test_genome_id_empty_fasta_input():
   assert rfm.genome_id("empty.fasta") == "empty"

def test_genome_id_valid_input():
   
   expected_output = "GCA_000007265.1_ASM726v1_genomic"
   
   output = rfm.genome_id("GCA_000007265.1_ASM726v1_genomic.fna")
   
   assert output == expected_output
    

##################################################################
    
    
#### TEST import_contigs ########################################
  
def test_import_contigs_no_args():
   with pytest.raises(TypeError):
       rfm.import_contigs()
      
expected_output = {"GCA-000007265-protein1_1":"ATGGTACAATATAACAATAATTATCCACAAGACAATAAGGAAGAAGCTATGACGGAAAACGAACAACTATTTTGGAATAGAGTACTAGAGCTATCTCGTTCTCAAATAGCACCAGCAGCTTATGAATTTTTTGTTCTAGAGGCTAGACTCCTCAAAATTGAACATCAAACTGCAGTTATTACTTTAGATAACATTGAAATGAAAAAGCTATTCTGGGAACAAAATTTGGGGCCTGTTATCCTAACAGCTGGTTTTGAAATTTTCAATGCTGAAATTACAGCTAACTATGTCTCAAACGATTTACATTTACAAGAAACTAGTTTTTCTAACTACCAGCAATCTAGCAATGAAGTAAATACTTTACCAATTAGAAAAATCGACTCTAACCTTAAAGAGAAATATACTTTTGCTAATTTTGTTCAAGGAGATGAAAATAGATGGGCTGTTTCAGCATCAATTGCTGTAGCTGATAGTCCTGGCACGACTTATAATCCTCTATTTATCTGGGGGGGACCTGGTCTAGGAAAGACGCATCTACTAAATGCTATTGGAAATCAAGTCTTAAGAGATAATCCAAACGCGAGGGTTTTATACATCACTGCTGAGAATTTTATTAATGAATTTGTCAGTCATATTCGTTTAGATTCGATGGAAGAATTAAAAGAAAAGTTTCGCAACTTAGACTTACTCCTGATTGATGATATTCAGTCGCTTGCTAAGAAAACCTTAGGGGGGACCCAAGAGGAGTTCTTCAATACTTTCAATGCTTTACATACAAACGATAAACAAATCGTATTGACCAGTGACCGAAATCCAAATCAATTAAATGATCTAGAAGAACGTCTAGTCACGCGCTTTAGTTGGGGACTCCCAGTAAATATCACACCACCTGATTTTGAAACACGAGTTGCTATTTTAACCAATAAAATTCAAGAATATCCTTATGATTTTCCTCAAGATACCATTGAATACTTAGCAGGAGAATTTGATTCCAACGTACGTGAATTAGAAGGAGCCTTGAAAAATATTAGTCTAGTTGCTGACTTTAAGCATGCTAAAACTATTACAGTAGATATAGCTGCAGAAGCTATCAGAGCACGTAAAAATGATGGTCCTATTGTTACTGTCATTCCTATAGAAGAAATTCAAATACAAGTTGGTAAATTCTATGGCGTAACTGTAAAAGAGATAAAAGCAACTAAAAGAACACAAGATATTGTCCTTGCAAGACAGGTAGCCATGTACTTAGCTCGTGAGATGACAGATAACAGTCTCCCAAAAATAGGTAAAGAATTTGGGGGACGAGATCACTCAACTGTTCTCCACGCTTATAATAAAATAAAAAATATGGTTGCTCAAGATGACAACTTACGAATTGAGATAGAAACTATCAAAAATAAAATCAGATAG"}


@xfail
def test_import_contigs_empty_list_input():
   rfm.import_contigs([])
   
def test_import_contigs_int_input():
   with pytest.raises(TypeError):
       rfm.import_contigs(9)

# Add check for file existence
#@xfail
def test_import_contigs_unexisting_txt_input():
   with pytest.raises(FileNotFoundError):
       rfm.import_contigs("unexisting.txt")

Add check for file contents ??
@xfail    
def test_import_contigs_existing_empty_txt_input():
   rfm.import_contigs("existing.txt")
   
        

    
def test_import_contigs_valid_input():
   
   expected_output = {"GCA-000007265-protein1_1":"ATGGTACAATATAACAATAATTATCCACAAGACAATAAGGAAGAAGCTATGACGGAAAACGAACAACTATTTTGGAATAGAGTACTAGAGCTATCTCGTTCTCAAATAGCACCAGCAGCTTATGAATTTTTTGTTCTAGAGGCTAGACTCCTCAAAATTGAACATCAAACTGCAGTTATTACTTTAGATAACATTGAAATGAAAAAGCTATTCTGGGAACAAAATTTGGGGCCTGTTATCCTAACAGCTGGTTTTGAAATTTTCAATGCTGAAATTACAGCTAACTATGTCTCAAACGATTTACATTTACAAGAAACTAGTTTTTCTAACTACCAGCAATCTAGCAATGAAGTAAATACTTTACCAATTAGAAAAATCGACTCTAACCTTAAAGAGAAATATACTTTTGCTAATTTTGTTCAAGGAGATGAAAATAGATGGGCTGTTTCAGCATCAATTGCTGTAGCTGATAGTCCTGGCACGACTTATAATCCTCTATTTATCTGGGGGGGACCTGGTCTAGGAAAGACGCATCTACTAAATGCTATTGGAAATCAAGTCTTAAGAGATAATCCAAACGCGAGGGTTTTATACATCACTGCTGAGAATTTTATTAATGAATTTGTCAGTCATATTCGTTTAGATTCGATGGAAGAATTAAAAGAAAAGTTTCGCAACTTAGACTTACTCCTGATTGATGATATTCAGTCGCTTGCTAAGAAAACCTTAGGGGGGACCCAAGAGGAGTTCTTCAATACTTTCAATGCTTTACATACAAACGATAAACAAATCGTATTGACCAGTGACCGAAATCCAAATCAATTAAATGATCTAGAAGAACGTCTAGTCACGCGCTTTAGTTGGGGACTCCCAGTAAATATCACACCACCTGATTTTGAAACACGAGTTGCTATTTTAACCAATAAAATTCAAGAATATCCTTATGATTTTCCTCAAGATACCATTGAATACTTAGCAGGAGAATTTGATTCCAACGTACGTGAATTAGAAGGAGCCTTGAAAAATATTAGTCTAGTTGCTGACTTTAAGCATGCTAAAACTATTACAGTAGATATAGCTGCAGAAGCTATCAGAGCACGTAAAAATGATGGTCCTATTGTTACTGTCATTCCTATAGAAGAAATTCAAATACAAGTTGGTAAATTCTATGGCGTAACTGTAAAAGAGATAAAAGCAACTAAAAGAACACAAGATATTGTCCTTGCAAGACAGGTAGCCATGTACTTAGCTCGTGAGATGACAGATAACAGTCTCCCAAAAATAGGTAAAGAATTTGGGGGACGAGATCACTCAACTGTTCTCCACGCTTATAATAAAATAAAAAATATGGTTGCTCAAGATGACAACTTACGAATTGAGATAGAAACTATCAAAAATAAAATCAGATAG"}
   
   output = rfm.import_contigs("GCA-000007265-protein1.fasta")
   
   assert output == expected_output

##############################################################


#### TEST extract_coding_sequences ########################################

orf_file_path = "GCA_000007265.1_ASM726v1_genomic.fna_ORF.txt"

#genome_id = rfm.genome_id()

contigs = rfm.import_contigs("GCA_000007265.1_ASM726v1_genomic.fna")

valid_output = rfm.extract_coding_sequences(orf_file_path, contigs, 0, 0)

#with open(orf_file_path, 'rb') as orf_file:
#    start_stop_codons = pickle.load(orf_file)
#    
#invalid_dict = {"INVALIDID": [[]]}

def test_extract_coding_sequences_no_args():
   with pytest.raises(TypeError):
       rfm.extract_coding_sequences()
        
def test_extract_coding_sequences_empty_list_input():
   with pytest.raises(TypeError):
       rfm.extract_coding_sequences([],[],[],[])
       
def test_extract_coding_sequences_invalid_orf_file():
   with pytest.raises(UnpicklingError):
       rfm.extract_coding_sequences("existing.txt", "","","")
       
def test_extract_coding_sequences_int_input():
   with pytest.raises(TypeError):
       rfm.extract_coding_sequences(7,7,7,5)
    
def test_extract_coding_sequences_empty_str_input():
   with pytest.raises(FileNotFoundError):
       rfm.extract_coding_sequences("","","","")
   
def test_extract_coding_sequences_no_dict_input():
   with pytest.raises(TypeError):
       rfm.extract_coding_sequences(orf_file_path, "", "", "")
       
def test_extract_coding_sequences_no_int_input():
   with pytest.raises(TypeError):
       rfm.extract_coding_sequences(orf_file_path, contigs, "", "")

def test_extract_coding_sequences_valid_input():
    assert rfm.extract_coding_sequences(orf_file_path, contigs, 0, 0) == valid_output


###########################################################################


#### TEST reverse_complement ########################################

def test_empty_reverse():
   with pytest.raises(TypeError):
       rfm.reverse_complement()


#@pytest.mark.xfail(raises=TypeError)
def test_non_str_reverse():
   with pytest.raises(TypeError):
       rfm.reverse_complement(9)
   
def test_amb():
   with pytest.raises(KeyError):
       rfm.reverse_complement("AAAAAACCCCTTTTGGCTTATCGX")

@xfail  
def test_empty_list():
   rfm.reverse_complement([])
   
def test_lower_case():
   rfm.reverse_complement("a")

def test_reverse_result():
   assert rfm.reverse_complement("AACTGGCATgctatgcat") == "ATGCATAGCATGCCAGTT"
   
def test_result_orientation():
   with pytest.raises(AssertionError):
       assert rfm.reverse_complement("AACTGGCATgctatgcat") == "TTGACCGTACGATACGTA"


#def test_result_orientation2():
#    assert rfm.reverse_complement("TTGACCGTACGATACGTA") == "TACGTATCGTACGGTCAA"
#


###########################################################################


#### TEST translate ######################################################

def test_translate_string():
   assert rfm.translate("Aello World") == "amb"

def test_seq_length():
   with pytest.raises(TranslationError):
       rfm.translate("ATGAAACCCCTTTTGGCTTATCGTAA")
   
def test_empty_list():
   with pytest.raises(AttributeError):
       rfm.translate([])
   
def test_lower_case():
   assert str(rfm.translate("atgggtcaataa")[0]) == "MGQ"

def test_translate_sense_strand():
   assert str(rfm.translate("ATGGGTCAATAA")[0]) == "MGQ"
   
def test_translate_antisense_strand():
   assert str(rfm.translate("TTATTGACCCAT")[0]) == "MGQ"

def test_translate_reverse_sense_strand():
   assert str(rfm.translate("AATAACTGGGTA")[0]) == "MGQ"
   
def test_translate_reverse_antisense_strand():
   assert str(rfm.translate("TACCCAGTTATT")[0]) == "MGQ"
   
def test_translate_atg_mid():
   assert str(rfm.translate("ATGGGTCAATTCATGAAGTAA")[0]) == "MGQFMK"
   
def test_translate_stop_mid():
   with pytest.raises(TranslationError):
       rfm.translate("ATGGGTTAACAATAA")
        
def test_translate_N():
   rfm.translate("ATGGTNCAATATAACAATAATTATCCACAAGACAATAAGGAAGAAGCTATGACGGAAAACGAACAACTATTTTGGAATAGAGTACTAGAGCTATCTCGTTCTCAAATAGCACCAGCAGCTTATGAATTTTTTGTTCTAGAGGCTAGACTCCTCAAAATTGAACATCAAACTGCAGTTATTACTTTAGATAACATTGAAATGAAAAAGCTATTCTGGGAACAAAATTTGGGGCCTGTTATCCTAACAGCTGGTTTTGAAATTTTCAATGCTGAAATTACAGCTAACTATGTCTCAAACGATTTACATTTACAAGAAACTAGTTTTTCTAACTACCAGCAATCTAGCAATGAAGTAAATACTTTACCAATTAGAAAAATCGACTCTAACCTTAAAGAGAAATATACTTTTGCTAATTTTGTTCAAGGAGATGAAAATAGATGGGCTGTTTCAGCATCAATTGCTGTAGCTGATAGTCCTGGCACGACTTATAATCCTCTATTTATCTGGGGGGGACCTGGTCTAGGAAAGACGCATCTACTAAATGCTATTGGAAATCAAGTCTTAAGAGATAATCCAAACGCGAGGGTTTTATACATCACTGCTGAGAATTTTATTAATGAATTTGTCAGTCATATTCGTTTAGATTCGATGGAAGAATTAAAAGAAAAGTTTCGCAACTTAGACTTACTCCTGATTGATGATATTCAGTCGCTTGCTAAGAAAACCTTAGGGGGGACCCAAGAGGAGTTCTTCAATACTTTCAATGCTTTACATACAAACGATAAACAAATCGTATTGACCAGTGACCGAAATCCAAATCAATTAAATGATCTAGAAGAACGTCTAGTCACGCGCTTTAGTTGGGGACTCCCAGTAAATATCACACCACCTGATTTTGAAACACGAGTTGCTATTTTAACCAATAAAATTCAAGAATATCCTTATGATTTTCCTCAAGATACCATTGAATACTTAGCAGGAGAATTTGATTCCAACGTACGTGAATTAGAAGGAGCCTTGAAAAATATTAGTCTAGTTGCTGACTTTAAGCATGCTAAAACTATTACAGTAGATATAGCTGCAGAAGCTATCAGAGCACGTAAAAATGATGGTCCTATTGTTACTGTCATTCCTATAGAAGAAATTCAAATACAAGTTGGTAAATTCTATGGCGTAACTGTAAAAGAGATAAAAGCAACTAAAAGAACACAAGATATTGTCCTTGCAAGACAGGTAGCCATGTACTTAGCTCGTGAGATGACAGATAACAGTCTCCCAAAAATAGGTAAAGAATTTGGGGGACGAGATCACTCAACTGTTCTCCACGCTTATAATAAAATAAAAAATATGGTTGCTCAAGATGACAACTTACGAATTGAGATAGAAACTATCAAAAATAAAATCAGATAG")


###########################################################################


#### TEST translate_coding_sequences #######################################################

def test_translate_coding_sequences_no_args():
   with pytest.raises(TypeError):
       rfm.translate_coding_sequences()
    
def test_translate_coding_sequences_empty_list():
   with pytest.raises(AttributeError):
       rfm.translate_coding_sequences([])

#@xfail      
#def test_translate_coding_sequences_empty_dict():
#    assert rfm.translate_coding_sequences({}) == []

def test_translate_coding_sequences_int_input():
   with pytest.raises(AttributeError):
       rfm.translate_coding_sequences(5)

def test_translate_coding_sequences_empty_str():
   with pytest.raises(AttributeError):
       rfm.translate_coding_sequences("")


def test_translate_coding_sequences_valid_input():
   assert rfm.translate_coding_sequences({1:"ATGGGTCAATAA"}) == [{1: 'MGQ'}, {1: 'sense'}, {}, {}]
   

def test_translate_coding_sequences_ambiguous_input():
   assert rfm.translate_coding_sequences({1:"ATGNGTCAATAA"}) == [{}, {}, {1: 'ATGNGTCAATAA'}, {}]


def test_translate_coding_sequences_untranslatable_input():
   assert rfm.translate_coding_sequences({1:"AAAAAAAAAAAA"}) == [{}, {}, {}, {1: 'AAAAAAAAAAAA'}]


###########################################################################



#### TEST determine_repeated #######################################################

coding_sequences_repeated = {1:"ATGGTCGAACGGTCGTCAAAGTCAGAGTACCCCGGGTACCAACTTACGGAGGATATTGCTTGCAGCTGCTAA",
                    2:"ATGGTCGAACGGTCGTCAAAGTCAGAGTACCCCGGGTACCAACTTACGGAGGATATTGCTTGCAGCTGCTAA",
                    3:"ATGTTCAAAGGCGAAAGAACAAAGGCTTACTGTGCGCAGAGGAACGCCCATTTTGCGGCTGGCGTCTATTGA"}

coding_sequences_no_repeats = {1:"ATGGTCGAACGGTCGTCAAAGTCAGAGTACCCCGGGTACCAACTTACGGAGGATATTGCTTGCAGCTGCTAA",
                               2:"ATGTTCAAAGGCGAAAGAACAAAGGCTTACTGTGCGCAGAGGAACGCCCATTTTGCGGCTGGCGTCTATTGA"}

proteins_repeated = rfm.translate_coding_sequences(coding_sequences_repeated)[0]

proteins_no_repeat = rfm.translate_coding_sequences(coding_sequences_no_repeats)[0]

def test_determine_repeated_no_args():
   with pytest.raises(TypeError):
       rfm.determine_repeated()

def test_determine_repeated_empty_list():
   with pytest.raises(AttributeError):
       rfm.determine_repeated([])

def test_determine_repeated_empty_str():
   with pytest.raises(AttributeError):
       rfm.determine_repeated("")

def test_determine_repeated_int_input():
   with pytest.raises(AttributeError):
       rfm.determine_repeated(5)

def test_determine_repeated_empty_dict():
   assert rfm.determine_repeated({}) == []

def test_determine_repeated_valid_input():
   assert rfm.determine_repeated(coding_sequences_repeated) == [[1,2]]
   
def test_determine_repeated_no_repeat_input():
   assert rfm.determine_repeated(coding_sequences_no_repeats) == []



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


def test_remove_repeated_no_args():
   with pytest.raises(TypeError):
       rfm.remove_repeated()
        
def test_remove_repeated_empty_list():
   assert rfm.remove_repeated([],[]) == []

def test_remove_repeated_empty_str():
   assert rfm.remove_repeated("","") == ""

#def test_remove_repeated_int_input():
#    with pytest.raises(TypeError):
#        rfm.remove_repeated(5,5)

def test_remove_repeated_empty_dict():
   assert rfm.remove_repeated({},[]) == {}

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

#rfm.determine_contained(proteins_inv)

def test_determine_contained_no_args():
   with pytest.raises(TypeError):
       rfm.determine_contained()

def test_determine_contained_empty_list():
   with pytest.raises(AttributeError):
       rfm.determine_contained([])
       
def test_determine_contained_empty_str():
   with pytest.raises(AttributeError):
       rfm.determine_contained("")
       
def test_determine_contained_empty_dict():
   assert rfm.determine_contained({}) == {}


def test_determine_contained_int_input():
   with pytest.raises(AttributeError):
       rfm.determine_contained(5)

def test_determine_contained_valid_input():
   assert rfm.determine_contained(proteins_inv) == {'MYYACGQIM':2,
                                                    'TYPLIEHSV':3,
                                                    'RPSRGTCTY':7,
                                                    'CGQIMQVERCGRSSR':5,
                                                    'IMQVERCGRSSRPSR':6,
                                                    'MQVERCGRSSRPSRGTCT':4}
   

def test_determine_contained_input_with_equal_proteins():
   assert rfm.determine_contained(proteins_eq_inv) == {'MYYACGQIM':2,
                                                'TYPLIEHSV':3,
                                                'RPSRGTCTY':7,
                                                'CGQIMQVERCGRSSR':5,
                                                'IMQVERCGRSSRPSR':6,
                                                'MQVERCGRSSRPSRGTCT':4}

        

############################################################################################
