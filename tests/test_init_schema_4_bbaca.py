#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 16:07:00 2019

@author: pcerqueira
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna, _verify_alphabet
import os
import argparse

from contextlib import nullcontext as does_not_raise
import pytest


def get_Short(genesList):
    """Creates a short version of each fasta file with only the 1st allele"""
    
    if not genesList:
        #print("An empty list was provided. Stopping execution...")
        #return None
        raise Exception("An empty list was provided. Stopping execution...")
    else:
        
        for gene in genesList:
            # gene = gene.rstrip('\n')
            pathtoDir = os.path.join(os.path.dirname(gene), "short")
            if not os.path.exists(pathtoDir):
                os.makedirs(pathtoDir)
            shortgene = os.path.join(os.path.dirname(gene), "short", os.path.basename(gene))
            shortgene = shortgene.replace(".fasta", "_short.fasta")
            
            
            first_allele = next(SeqIO.parse(gene, "fasta", IUPAC.unambiguous_dna))
            
            if not _verify_alphabet(first_allele.seq.upper()):
                print("The DNA sequence has invalid nucleotides. Execution will not be stopped.")
                raise Exception("The DNA sequence has invalid nucleotides. Execution will not be stopped.")
                
            else:
                with open(shortgene, "w") as fG:
                    fG.write('>' + str(first_allele.id) + '\n' + str(first_allele.seq.upper()) + '\n')
            

    
#            with open(shortgene, "w") as fG:
#                first_allele = next(SeqIO.parse(gene, "fasta", generic_dna))
#                fG.write('>' + str(first_allele.id) + '\n' + str(first_allele.seq.upper()) + '\n')
    
            #gene_fp2 = HTSeq.FastaReader(gene)
            # for allele in SeqIO.parse(gene, "fasta", generic_dna):
            #     fG = open(shortgene, 'w')
            #     fG.write('>' + str(allele.id) + '\n' + str(allele.seq.upper()) + '\n')
            #     fG.close()
            #     break

    return True


def main(geneFiles):
    #~ parser = argparse.ArgumentParser(
        #~ description="This program prepares a schema for a chewBBACA allele call, creating a short version of each fast with only the 1st allele")
    #~ parser.add_argument('-i', nargs='?', type=str, help='List of genes files (list of fasta files)', required=True)
#~ 
    #~ args = parser.parse_args()

    geneFiles = args.i
    listGenes = []
    gene_fp = open(geneFiles, 'r')
    for gene in gene_fp:
        gene = gene.rstrip('\n')
        listGenes.append(gene)
    gene_fp.close()

    get_Short(listGenes)


#if __name__ == "__main__":
#    main()


expected_content = ">GCA-000007265-protein1_1\nATGGTACAATATAACAATAATTATCCACAAGACAATAAGGAAGAAGCTATGACGGAAAACGAACAACTATTTTGGAATAGAGTACTAGAGCTATCTCGTTCTCAAATAGCACCAGCAGCTTATGAATTTTTTGTTCTAGAGGCTAGACTCCTCAAAATTGAACATCAAACTGCAGTTATTACTTTAGATAACATTGAAATGAAAAAGCTATTCTGGGAACAAAATTTGGGGCCTGTTATCCTAACAGCTGGTTTTGAAATTTTCAATGCTGAAATTACAGCTAACTATGTCTCAAACGATTTACATTTACAAGAAACTAGTTTTTCTAACTACCAGCAATCTAGCAATGAAGTAAATACTTTACCAATTAGAAAAATCGACTCTAACCTTAAAGAGAAATATACTTTTGCTAATTTTGTTCAAGGAGATGAAAATAGATGGGCTGTTTCAGCATCAATTGCTGTAGCTGATAGTCCTGGCACGACTTATAATCCTCTATTTATCTGGGGGGGACCTGGTCTAGGAAAGACGCATCTACTAAATGCTATTGGAAATCAAGTCTTAAGAGATAATCCAAACGCGAGGGTTTTATACATCACTGCTGAGAATTTTATTAATGAATTTGTCAGTCATATTCGTTTAGATTCGATGGAAGAATTAAAAGAAAAGTTTCGCAACTTAGACTTACTCCTGATTGATGATATTCAGTCGCTTGCTAAGAAAACCTTAGGGGGGACCCAAGAGGAGTTCTTCAATACTTTCAATGCTTTACATACAAACGATAAACAAATCGTATTGACCAGTGACCGAAATCCAAATCAATTAAATGATCTAGAAGAACGTCTAGTCACGCGCTTTAGTTGGGGACTCCCAGTAAATATCACACCACCTGATTTTGAAACACGAGTTGCTATTTTAACCAATAAAATTCAAGAATATCCTTATGATTTTCCTCAAGATACCATTGAATACTTAGCAGGAGAATTTGATTCCAACGTACGTGAATTAGAAGGAGCCTTGAAAAATATTAGTCTAGTTGCTGACTTTAAGCATGCTAAAACTATTACAGTAGATATAGCTGCAGAAGCTATCAGAGCACGTAAAAATGATGGTCCTATTGTTACTGTCATTCCTATAGAAGAAATTCAAATACAAGTTGGTAAATTCTATGGCGTAACTGTAAAAGAGATAAAAGCAACTAAAAGAACACAAGATATTGTCCTTGCAAGACAGGTAGCCATGTACTTAGCTCGTGAGATGACAGATAACAGTCTCCCAAAAATAGGTAAAGAATTTGGGGGACGAGATCACTCAACTGTTCTCCACGCTTATAATAAAATAAAAAATATGGTTGCTCAAGATGACAACTTACGAATTGAGATAGAAACTATCAAAAATAAAATCAGATAG\n"

def test_get_Short_valid_input(tmpdir):
    """Tests the output of the get_Short function"""
    p = tmpdir.join("expected_output.fasta")
    p.write(expected_content)
    
    get_Short(["/home/pcerqueira/Lab_Software/testing/tests/GCA-000007265-protein1.fasta"])
    output_file = "/home/pcerqueira/Lab_Software/testing/tests/short/GCA-000007265-protein1_short.fasta"
    
    with open(output_file, "r") as out:
        assert p.read() == out.read()
    

@pytest.mark.parametrize(
        "test_input, expectation",
        [([], pytest.raises(Exception)), # Tests empty list input
         ([5], pytest.raises(TypeError)), # Tests integer input
         (["ola"], pytest.raises(FileNotFoundError)), # Tests string input
         (["/home/pcerqueira/Lab_Software/testing/tests/unexisting.fasta"], pytest.raises(FileNotFoundError)),  # Tests unexisting fasta file input
         (["/home/pcerqueira/Lab_Software/testing/tests/invalid.fasta"], pytest.raises(Exception)),             # Tests invalid fasta file input
         (["/home/pcerqueira/Lab_Software/testing/tests/GCA-000007265-protein1.fasta"], does_not_raise())       # Tests valid input to confirm that no exception is raised
         ])
def test_get_Short_invalid_inputs(test_input, expectation):
    """Tests the behaviour of the get_Short function with unexpected inputs"""
    with expectation:
        get_Short(test_input)


 


