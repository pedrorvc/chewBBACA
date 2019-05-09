#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 08 13:50:17 2019

@author: pcerqueira
"""


import os
import pickle
import subprocess

import pytest

def main(input_file,tempPath,choosenTaxon):

    contigsFasta = input_file
    basepath = tempPath


    # ------------ #
    # RUN PRODIGAL #
    # ------------ #


    if choosenTaxon == "False":
    
        proc = subprocess.Popen(
            ['prodigal', '-i', contigsFasta, '-c', '-m', '-g', '11', '-p', 'single', '-f', 'sco', '-q'],
            stdout=subprocess.PIPE)
    else:
        proc = subprocess.Popen(
            ['prodigal', '-i', contigsFasta, '-c', '-m', '-g', '11', '-p', 'single', '-f', 'sco', '-q', '-t',
                choosenTaxon], stdout=subprocess.PIPE)
    
    cdsDict = {}
    tempList = []
    # Reads the stdout from Prodigal
    prodigal_out = proc.stdout.readlines()

    # Parse 'Prodigal's output 
    for line in prodigal_out:
        line_decoded = line.decode("utf-8")
        
        if "seqhdr" in line_decoded:
            seqid = line_decoded.split('"')[1].split()[0]

        # Obtain the start and end positions of the CDSs
        elif ">" in line_decoded:
            cdsL = line_decoded.split("_")
            
            # Start index correction needed because Prodigal indexes start in 1 instead of 0
            start_position = int(cdsL[1]) - 1
            end_position = int(cdsL[2])
            
            # Strand of the CDS. 1 if sense (+), else 0 (-) (antisense)
            strand = 1 if cdsL[-1].strip() == "+" else 0
            tempList = [start_position, end_position, strand]
        
            if len(tempList) > 0:
                if seqid not in cdsDict:
                    # Add the sequence ID as the key and the list of lists of the CDSs' positions as the value
                    cdsDict[seqid] = [tempList]
                else:
                    cdsDict[seqid] += [tempList]
                    tempList = []
                
#    return cdsDict
    
#    Add the cdsDict to a file
    filepath = os.path.join(basepath, str(os.path.basename(contigsFasta)) + "_ORF.txt")
    #print(filepath)
    with open(filepath, 'wb') as f:
        var = cdsDict
        pickle.dump(var, f)
    
    print("done prodigal run on: " + str(os.path.basename(contigsFasta)))

    return True


#
#
#files = ["/home/pcerqueira/Lab_Software/testing/small_fasta_files/GCA_000007265.1_ASM726v1_genomic.fna",
#         "/home/pcerqueira/Lab_Software/testing/inv_small_fasta/inv_GCA_000007265.1_ASM726v1_genomic.fna",
#         "/home/pcerqueira/Lab_Software/testing/tests/invalid.fasta",
#         "/home/pcerqueira/Lab_Software/testing/tests/empty.fasta"]
#
#tempPath = ["/home/pcerqueira/Lab_Software/testing/Temp_meu",
#            "/home/pcerqueira/Lab_Software/testing/Temp_inv",
#            "/home/pcerqueira/Lab_Software/testing/Temp_test"]
#
#training_file = "/home/pcerqueira/Lab_Software/testing/Streptococcus_agalactiae.trn"
#
#expected_output_files = ["/home/pcerqueira/Lab_Software/testing/Temp_test/GCA_000007265.1_ASM726v1_genomic.fna_ORF.txt",
#                "/home/pcerqueira/Lab_Software/testing/Temp_test/inv_GCA_000007265.1_ASM726v1_genomic.fna_ORF.txt"]
#
#output_files = [tempPath[0] + "/GCA_000007265.1_ASM726v1_genomic.fna_ORF.txt",
#                tempPath[1] + "/inv_GCA_000007265.1_ASM726v1_genomic.fna_ORF.txt",
#                tempPath[2] + "/invalid.fasta_ORF.txt",
#                tempPath[3] + "/empty.fasta_ORF.txt"]
#
#valid_input = (files[0], tempPath[0], training_file, expected_output_files[0], )

@pytest.mark.parametrize(
        "input_file, tempPath, choosenTaxon, expectation",
        [([],[],[], pytest.raises(TypeError)),
         (7,7,7, pytest.raises(TypeError)),
         (pytest.param("","","", "", marks=pytest.mark.xfail))
         ])
def test_prodigal_invalid_inputs(input_file, tempPath, choosenTaxon, expectation):
    with expectation:
        main(input_file, tempPath, choosenTaxon, expectation)





#@pytest.mark.skip
def test_output_prodigal():
#    p = tmpdir.join("output_file.txt")
#    p.write(expected_output)
    
    training = "/home/pcerqueira/Lab_Software/testing/Streptococcus_agalactiae.trn"
    tempPath = "/home/pcerqueira/Lab_Software/testing/Temp_meu"
    file = "/home/pcerqueira/Lab_Software/testing/small_fasta_files/GCA_000007265.1_ASM726v1_genomic.fna"


    expected_output_file = "/home/pcerqueira/Lab_Software/testing/Temp_test/GCA_000007265.1_ASM726v1_genomic.fna_ORF.txt"
    
    with open(expected_output_file, "rb") as ex_out:
        expected_output = pickle.load(ex_out)
    
    main(file, tempPath, training)
    
    
    output_file = tempPath + "/GCA_000007265.1_ASM726v1_genomic.fna_ORF.txt"
    
    assert os.path.isfile(output_file)
    
    with open(output_file, "rb") as out:
        output = pickle.load(out)
        
    assert expected_output == output

#@pytest.mark.skip
def test_prodigal_inverted_input():
    training = "/home/pcerqueira/Lab_Software/testing/Streptococcus_agalactiae.trn"
    tempPath = "/home/pcerqueira/Lab_Software/testing/Temp_inv"
    file = "/home/pcerqueira/Lab_Software/testing/inv_small_fasta/inv_GCA_000007265.1_ASM726v1_genomic.fna"


    expected_output_file = "/home/pcerqueira/Lab_Software/testing/Temp_test/inv_GCA_000007265.1_ASM726v1_genomic.fna_ORF.txt"
    
    with open(expected_output_file, "rb") as ex_out:
        expected_output = pickle.load(ex_out)
    
    main(file, tempPath, training)
    
    output_file = tempPath + "/inv_GCA_000007265.1_ASM726v1_genomic.fna_ORF.txt"
    
    assert os.path.isfile(output_file)
    
    with open(output_file, "rb") as out:
        output = pickle.load(out)
        
    assert expected_output == output
    
#@pytest.mark.skip
def test_prodidgal_invalid_fasta_input():
    training = "/home/pcerqueira/Lab_Software/testing/Streptococcus_agalactiae.trn"
    tempPath = "/home/pcerqueira/Lab_Software/testing/Temp_test"
    file = "/home/pcerqueira/Lab_Software/testing/tests/invalid.fasta"
    
    expected_output = {}
    
    main(file, tempPath, training)
    
    output_file = tempPath + "/invalid.fasta_ORF.txt"
    
    assert os.path.isfile(output_file)
    
    with open(output_file, "rb") as out:
        output = pickle.load(out)
        
    assert expected_output == output

#@pytest.mark.skip    
def test_prodigal_empty_fasta_input(capfd):
    training = "/home/pcerqueira/Lab_Software/testing/Streptococcus_agalactiae.trn"
    tempPath = "/home/pcerqueira/Lab_Software/testing/Temp_test"
    file = "/home/pcerqueira/Lab_Software/testing/tests/empty.fasta"
    
    expected_output = {}
    
    main(file, tempPath, training)
    
    out, err = capfd.readouterr()
    
    output_file = tempPath + "/empty.fasta_ORF.txt"
    
    assert os.path.isfile(output_file)
    
    with open(output_file, "rb") as out_file:
        output = pickle.load(out_file)
    
    assert "Error:  no input sequences to analyze." in err    
    assert expected_output == output    