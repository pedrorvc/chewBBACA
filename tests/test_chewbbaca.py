#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 11:07:50 2019

@author: pcerqueira
"""

from unittest.mock import patch
import sys
import os
import argparse
import pytest
import chewBBACA
import pickle

from Bio import SeqIO

def is_fasta(filename):
    """ Checks if a file is a FASTA file."""
    
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        
        # returns True if FASTA file, False otherwise
        return any(fasta)

#import pytest_mock 

#chewBBACA.create_schema()


#@mock.patch('argparse.ArgumentParser.parse_args',
#            return_value=argparse.Namespace("CreateSchema",
#                                            "/home/pcerqueira/Lab_Software/testing/small_fasta_files",
#                                            "/home/pcerqueira/Lab_Software/testing/chewbbaca_test",
#                                            6,
#                                            "/home/pcerqueira/SW/anaconda3/bin/blastp",
#                                            False,
#                                            0.6,
#                                            "/home/pcerqueira/Lab_Software/testing/Streptococcus_agalactiae.trn",
#                                            False,
#                                            201))




#with patch.object(sys, "argv", testargs):
#    chewBBACA.main()



@pytest.fixture
def gen_chewie_args():
    
    testargs = ["chewBBACA.py",
            "CreateSchema",
            "-i",
            "/home/pcerqueira/Lab_Software/testing/small_fasta_files/",
            "-o",
            "/home/pcerqueira/Lab_Software/testing/schema_seed_mock",
            "--cpu",
            "6",
            "--ptf",
            "/home/pcerqueira/Lab_Software/testing/Streptococcus_agalactiae.trn"]    
    return testargs


@pytest.fixture
def gen_chewie_no_args():
    
    testargs = ["chewBBACA.py",
            "CreateSchema"]    
    return testargs


@pytest.fixture
def gen_chewie_args_empty_fasta_input():
    
    testargs = ["chewBBACA.py",
            "CreateSchema",
            "-i",
            "/home/pcerqueira/Lab_Software/testing/tests/empty.fasta",
            "-o",
            "/home/pcerqueira/Lab_Software/testing/schema_seed_mock",
            "--cpu",
            "6",
            "--ptf",
            "/home/pcerqueira/Lab_Software/testing/Streptococcus_agalactiae.trn"]    
    return testargs

@pytest.fixture
def gen_chewie_args_invalid_fasta_input():
    
    testargs = ["chewBBACA.py",
            "CreateSchema",
            "-i",
            "/home/pcerqueira/Lab_Software/testing/tests/invalid_fastas/",
            "-o",
            "/home/pcerqueira/Lab_Software/testing/schema_seed_mock",
            "--cpu",
            "6",
            "--ptf",
            "/home/pcerqueira/Lab_Software/testing/Streptococcus_agalactiae.trn"]    
    return testargs





# @pytest.mark.parametrize()
# def test_chewie():
#     """
#     """


# @pytest.mark.skip
def test_chewbbaca_create_schema(gen_chewie_args):
    
    with patch.object(sys, "argv", gen_chewie_args):
        chewBBACA.main()
    
    short_path = "/home/pcerqueira/Lab_Software/testing/schema_seed_mock/short"

    for fasta in os.listdir(short_path):
        with open(os.path.join(short_path, fasta)) as f:
            assert len(f.read()) != 0

        assert is_fasta(os.path.join(short_path, fasta)) == True

        fasta_file = SeqIO.read(os.path.join(short_path, fasta), "fasta")
        assert str(fasta_file) != ""


    assert len(os.listdir("/home/pcerqueira/Lab_Software/testing/schema_seed_mock")) == 2163
    assert len(os.listdir("/home/pcerqueira/Lab_Software/testing/schema_seed_mock/short")) == 2162





@pytest.mark.skip
# @pytest.mark.xfail
def test_chewbbaca_create_schema_no_args(gen_chewie_no_args, capfd):
   
   with patch.object(sys, "argv", gen_chewie_no_args):
       chewBBACA.main()
       out, err = capfd.readouterr()

   assert "chewBBACA.py: error: the following arguments are required: -i, -o, --cpu" in err
       
       
   
@pytest.mark.skip
def test_chewbbaca_create_schema_empty_fasta_input(gen_chewie_args_empty_fasta_input, capfd):
   
   with patch.object(sys, "argv", gen_chewie_args_empty_fasta_input):
       chewBBACA.main()
       out, err = capfd.readouterr()


   assert "Error: can't open input file" in err


@pytest.mark.skip
def test_chewbbaca_create_schema_invalid_fasta_input(gen_chewie_args_invalid_fasta_input, capfd):
    
    prodigal_output_file_1 = "/home/pcerqueira/Lab_Software/testing/temp/invalid.fasta_ORF.txt" 
    prodigal_output_file_2 = "/home/pcerqueira/Lab_Software/testing/temp/invalid2.fasta_ORF.txt"

    expected_prodigal_output_1 = {}
    expected_prodigal_output_2 = {}
    
    with patch.object(sys, "argv", gen_chewie_args_invalid_fasta_input):
        chewBBACA.main()
        out, err = capfd.readouterr()
    
    with open(prodigal_output_file_1, "rb") as out_file1:
        prodigal_output_1 = pickle.load(out_file1)
        
    with open(prodigal_output_file_2, "rb") as out_file2:
        prodigal_output_2 = pickle.load(out_file2)
        

    assert expected_prodigal_output_1 == prodigal_output_1
    assert expected_prodigal_output_2 == prodigal_output_2
    
    assert "BLAST Database error" in err










