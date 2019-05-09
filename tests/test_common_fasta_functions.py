#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 15:35:46 2019

@author: pcerqueira
"""

from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
import os
import shutil

from contextlib import nullcontext as does_not_raise
import pytest


def ensure_dir(f):
    if not os.path.isdir(f):
        # print "No BLAST db dir was found. Creating it..."
        os.makedirs(f)


# overwrite has to be 1 or 0 (True or False)
def Create_Blastdb(questionDB, overwrite, dbtypeProt):
    """questionDB is a protein fasta file
       overwrite is 1 or 0
       dbtypeProt is a bool"""
    base = os.path.basename(questionDB)
    dirname = os.path.dirname(questionDB)
    isProt = dbtypeProt

    if len(dirname) == 0:
        dirname = '.'
    basename = os.path.splitext(base)[0]
    ensure_dir(dirname + "/blastdbs")
    name = dirname + "/blastdbs/" + basename + "_db"

    if not os.path.isfile(name + ".nin") and not os.path.isfile(name + ".nhr") and not os.path.isfile(name + ".nsq"):

        if not isProt:
            os.system(
                "makeblastdb -in " + questionDB + " -out " + name + " -dbtype nucl -logfile " + name + "_blast.log")
        else:
            os.system(
                "makeblastdb -in " + questionDB + " -out " + name + " -dbtype prot -logfile " + name + "_blast.log")

    elif overwrite:
        if not isProt:
            os.system(
                "makeblastdb -in " + questionDB + " -out " + name + " -dbtype nucl -logfile " + name + "_blast.log")
        else:
            os.system(
                "makeblastdb -in " + questionDB + " -out " + name + " -dbtype prot -logfile " + name + "_blast.log")

    else:
        print("BLAST DB files found. Using existing DBs..")
    return name


#def runBlast(cline, bOutFile, locus_sbjct):
#    os.system(str(cline))
#    rec = open(bOutFile)
#    blast_record = NCBIXML.read(rec)
#
#    if os.path.isfile(locus_sbjct):
#        os.remove(locus_sbjct)
#    os.remove(bOutFile)
#
#    return blast_record


def runBlastParser(cline, bOutFile):
    """Ensure cline is in fact the blastp command"""
    if str(shutil.which("blastp")) in str(cline):
        os.system(str(cline)) 
#        print("opening xml")
        rec = open(bOutFile)
#        print("parsing...")
        blast_records = NCBIXML.parse(rec)
    else:
        raise Exception("Blastp path not found in command line argument")

    #	if os.path.isfile(locus_sbjct):
    #		os.remove(locus_sbjct)

    # os.remove(bOutFile)

    return blast_records



#expected_output = SeqIO.parse(open("/home/pcerqueira/Lab_Software/testing/temp_copy/protogenome0/protogenome0_proteins.fasta"),'fasta')


os.chdir("/home/pcerqueira/Lab_Software/testing/tests")

@pytest.mark.parametrize(
        "questionDB, overwrite, dbtypeProt, expectation",
        [([],[],[], pytest.raises(TypeError)),                                  # Test empty list input
         (pytest.param("","","","", marks=pytest.mark.xfail)),                  # Test empty string input
         (pytest.param("unexisting.fasta","","","", marks=pytest.mark.xfail)),  # Test unexisting fasta file input
         (pytest.param("invalid.fasta","","","", marks=pytest.mark.xfail)),     # Test invalid fasta file. The file contains a valid id and an invalid DNA/protein sequence 
         (pytest.param("empty.fasta","","","", marks=pytest.mark.xfail)),       # Test empty fasta file
         ("/home/pcerqueira/Lab_Software/testing/temp_copy/protogenome0/protogenome0_proteins.fasta", 0, True, does_not_raise()) # Test valid input to confirm that no exception is raised
         ])
def test_Create_Blastdb_invalid_inputs(questionDB, overwrite, dbtypeProt, expectation):
    """Tests the behaviour of the Create_Blastdb function with unexpected inputs"""
    with expectation:
        Create_Blastdb(questionDB, overwrite, dbtypeProt)



# os.chdir("/home/pcerqueira/Lab_Software/testing/temp_copy/protogenome0")


# cmd_line = NcbiblastpCommandline(cmd=shutil.which("blastp"), query="/home/pcerqueira/Lab_Software/testing/temp_copy/protogenome0/protogenome0_proteins.fasta"
#                               , db="/home/pcerqueira/Lab_Software/testing/temp_copy/protogenome0/blastdbs/protogenome0_proteins_db",
#                               evalue=0.001,
#                               out="/home/pcerqueira/Lab_Software/testing/temp_copy/protogenome0/output.xml", outfmt=5, num_threads=1)

# @pytest.mark.parametrize(
#         "cline, bOutFile, expectation",
#         [([],[], pytest.raises(Exception)), # Test empty list input
#          ("hello world", "oops", pytest.raises(Exception)), # Test invalid string input
#          (cmd_line, "oops", pytest.raises(FileNotFoundError)), # Test command line and invalid string input
#          (pytest.param(cmd_line,"/home/pcerqueira/Lab_Software/testing/temp_copy/protogenome0/output_empty.xml", "", marks=pytest.mark.xfail)), # Test command line and empty xml file input
#          (cmd_line, "/home/pcerqueira/Lab_Software/testing/temp_copy/protogenome0/output.xml", does_not_raise()) # Test valid input to confirm that no exception is raised
#          ])
# def test_runBlastParser_invalid_inputs(cline, bOutFile, expectation):
#     """Tests the behaviour of the runBlastParser function with unexpected inputs"""
#     with expectation:
#         runBlastParser(cline, bOutFile)









