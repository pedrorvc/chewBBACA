#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  8 12:14:00 2019

@author: pcerqueira
"""

from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
import os
import shutil

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

def test_Create_Blastdb_no_args():
    with pytest.raises(TypeError):
        Create_Blastdb()

def test_Create_Blastdb_empty_args():
    with pytest.raises(TypeError):
        Create_Blastdb([],[],[])
        
@pytest.mark.xfail
def test_Create_Blastdb_empty_str_args():    
    Create_Blastdb("","","")

@pytest.mark.xfail    
def test_Create_Blastdb_unexisting_fasta():
    Create_Blastdb("ola.fasta", "", "")

@pytest.mark.xfail        
def test_Create_Blastdb_invalid_fasta():
    Create_Blastdb("invalid.fasta", "", "")
    
@pytest.mark.xfail
def test_Create_Blastdb_empty_fasta():
    Create_Blastdb("empty.fasta", "", "")
   
def test_Create_Blastdb_valid_fasta():
    Create_Blastdb("/home/pcerqueira/Lab_Software/testing/temp_copy/protogenome0/protogenome0_proteins.fasta",
                   0,
                   True)




#
#@pytest.fixture()
#def blastp_command():
#    return NcbiblastpCommandline(cmd=shutil.which("blastp"), query="/home/pcerqueira/Lab_Software/testing/temp_copy/protogenome0/protogenome0_proteins.fasta",
#                                 db="/home/pcerqueira/Lab_Software/testing/temp_copy/protogenome0/blastdbs/protogenome0_proteins_db",
#                                 evalue=0.001,
#                                 out="/home/pcerqueira/Lab_Software/testing/temp_copy/protogenome0/output.xml", 
#                                 outfmt=5, num_threads=1)


os.chdir("/home/pcerqueira/Lab_Software/testing/temp_copy/protogenome0")


cmd_line = NcbiblastpCommandline(cmd=shutil.which("blastp"), query="/home/pcerqueira/Lab_Software/testing/temp_copy/protogenome0/protogenome0_proteins.fasta"
                              , db="/home/pcerqueira/Lab_Software/testing/temp_copy/protogenome0/blastdbs/protogenome0_proteins_db",
                              evalue=0.001,
                              out="/home/pcerqueira/Lab_Software/testing/temp_copy/protogenome0/output.xml", outfmt=5, num_threads=1)


def test_runBlastParser_no_args():
    with pytest.raises(TypeError):
        runBlastParser()
     
def test_runBlastParser_diff_args():
    with pytest.raises(Exception):
        runBlastParser("hello world", "oops")

def test_runBlastParser_no_blastOut():
    with pytest.raises(FileNotFoundError):
        runBlastParser(cmd_line, "oops")

@pytest.mark.xfail       
def test_runBlastParser_empty_blastOut_file():
    runBlastParser(cmd_line,
                   "/home/pcerqueira/Lab_Software/testing/temp_copy/protogenome0/output_empty.xml")


def test_runBlastParser_normal_run():
        runBlastParser(cmd_line,
                   "/home/pcerqueira/Lab_Software/testing/temp_copy/protogenome0/output.xml")
