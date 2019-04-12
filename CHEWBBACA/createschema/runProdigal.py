#!/usr/bin/env python3
import sys
import os
import subprocess
import pickle
#import shutil

def main(input_file,tempPath,choosenTaxon):


    # choosenTaxon = "/home/pcerqueira/Lab_Software/testing/Streptococcus_agalactiae.trn"
    # basepath = "/home/pcerqueira/Lab_Software/testing/temp"

    contigsFasta = input_file
    #genomes = os.listdir("/home/pcerqueira/Lab_Software/testing/fasta_files")
    # for genome in genomes:
    #     contigsFasta = "/home/pcerqueira/Lab_Software/testing/fasta_files/" + genome

    #contigsFasta = "/home/pcerqueira/Lab_Software/testing/fasta_files/GCA_000007265.1_ASM726v1_genomic.fna"

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
    prodigal_out = proc.stdout.readlines()

    # Parse Prodigal's output 
    for line in prodigal_out:
        line_decoded = line.decode("utf-8")

        # The first line of the prodigal output contains the sequence ID
        if line_decoded.startswith("#") and "seqhdr" in line_decoded:
            seqid = line_decoded.split('"')[1].split()[0]

        # The second line of the prodigal output is not used 
        elif line_decoded.startswith("#") and "model" in line_decoded:
            continue
        
        # Obtain the start and end positions of the CDSs
        else:
            cdsL = line_decoded.split("_")
            start_position = int(cdsL[1]) - 1
            end_position = int(cdsL[2])
            tempList.append([start_position, end_position])
    
    # Add the sequence ID as the key and the list of lists of the CDSs' positions as the value
    cdsDict[seqid] = tempList
    
    # Write a file with the cdsDict
    filepath = os.path.join(basepath, str(os.path.basename(contigsFasta)) + "_ORF.txt")
    with open(filepath, 'wb') as f:
        var = cdsDict
        pickle.dump(var, f)
    
    print("done prodigal run on: " + str(os.path.basename(contigsFasta)))

    return True


if __name__ == "__main__":
    main()
