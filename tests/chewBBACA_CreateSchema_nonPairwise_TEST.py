#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 20:27:53 2019

@author: rfm
"""


###############################################################################
#                     CreateSchema - nonPairwise TESTS                        #
###############################################################################

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import os
import argparse
import time
import pickle
import shutil
import multiprocessing
import subprocess

#from get_varSize_deep import *


def genome_id(fasta_path):
    """ Get genome/assembly identifier from full path.
    
        Args: 
            fasta_path (str): the full path to the FASTA file.
        
        Returns: 
            file (str): the last string resulting from the splitting of the path by 
            "/" and stripped of the file extension suffix.
        
        Example:
            
            >>> get_id("/home/user/wd/file.fasta")
            file
        
    """
    
    # possible FASTA file suffixes
    extension_format = ['.fasta','.fna','.ffn']
    
    file = fasta_path.split('/')[-1]
    
    for extension in extension_format:
        if extension in file:
            file = file.replace(extension, '')
            
    return file
    

# Import contigs for each genome
def import_contigs(fasta_path):
    """ Imports contigs from a FASTA file.
    
        Args: 
            fasta_path (str): full path to the FASTA file.
        
        Returns: 
            dictionary that has contigs ids as keys and contigs DNA 
            sequences as values.
            
        Example:
            
            >>> import_contigs("/home/user/wd/file.fasta")
            {contig_1:'TCGAACCACCGACCTCACGCTTATCAGG...',
            contig_2:'ATAAATGGCGCGAGACGGAATCGAACCGCCGA...'
            ...}
            
    """
    
    contigs_dict = {}
    # use BioPython to read FASTA file and get each contig sequence
    for contig in SeqIO.parse(fasta_path, 'fasta', generic_dna):
        # seq object has to be converted to string
        sequence = str(contig.seq.upper())
        contig_id = contig.id
        
        # add contig id as key and DNA sequence as value
        contigs_dict[contig_id] = sequence
    
    return contigs_dict


# Extract CDSs from contigs
def extract_coding_sequences(reading_frames, contigs, starting_id, genome_id):
    """ Extracts CDSs from contigs based on the start codon and stop codon 
        positions determined by Prodigal.
        
        Args:
            orf_file_path (str): full path to the ORF file derived from Prodigal.
            contigs (dict): a dictionary with contigs ids as keys and contigs 
            sequences as values.
            starting_protid (int): integer identifier to give to the first CDS
            extracted and that will be incremented to serve as identifier for
            subsequent CDSs.
            genome_id (str): id of the genome or assembly.
            
        Returns:
            cds_dict (dict): dictionary with CDSs ids as keys and CDSs DNA 
            sequences as values.
            cdss_lines (list): list of lists where each sublist has information 
            about the CDS.
            cdss_contigs (dict): dictionary with CDSs ids as keys and contigs 
            ids as values.
            protid (int): last extracted CDS id + 1. integer value that might 
            serve as starting id if the function will be used to extract CDSs 
            from more genomes/assemblies.
        
        Example:
            >>> extract_cdss("/home/user/wd/file_ORF.txt", contigs, 1, "genome1")
            {1: 'GTGACACCAAAACCAGTACCGGATAAA...',
            2: 'TTAGTCTAATTCTATTTGGAGAAATTTA...',
            ...},
            [['genome1','contig_1','102','528','1'],
            ['genome1','contig_1','1861','2110','2'],
            ...],
            {1: 'contig_1',
            2: 'contig_1,
            ...},
            1914
            
    """
    
    # load binary file with a list of lists for each contig
    # each sublist has a start codon and stop codon positions in the contig
    with open(reading_frames, 'rb') as orf_file:
        start_stop_codons = pickle.load(orf_file)

    protid = starting_id
    coding_sequences = {}
    coding_sequences_info = []
    coding_sequences_contig = {}
    coding_sequences_genome = {}
    # for each contig
    for contig_id, frames in start_stop_codons.items():
        # for each start and stop codon in that contig
        for coding_sequence in frames:
            start_codon = coding_sequence[0]
            stop_codon = coding_sequence[1]
            # extract CDS sequence
            cds_sequence = contigs[contig_id][start_codon:stop_codon].upper()
            
            # store CDS with unique id
            coding_sequences[protid] = cds_sequence
            
            # store CDS information
            coding_sequences_info.append([genome_id, contig_id, str(start_codon), str(stop_codon), str(protid)])
            
            # store the contig id from where each CDS was extracted
            coding_sequences_contig[protid] = contig_id
            
            # store the genome id from where each CDS was extracted
            coding_sequences_genome[protid] = genome_id
            
            # increment the CDS id by 1 so that it can be an unique identifier
            # for the next CDS
            protid += 1
            
    return [coding_sequences, coding_sequences_info, coding_sequences_contig, coding_sequences_genome, protid]


# translate CDSs
def reverse_complement(dna_sequence):
    """ Determines the reverse complement of given DNA strand.
    
        Args:
            strDNA (str): string representing a DNA sequence.
        
        Returns:
            revC_dna (str): the reverse complement of the DNA sequence, without
            lowercase letters.
        
        Example:
            >>> reverse_complement('ATCGgcaNn')
            'NNTGCCGAT'
            
    """
    
    base_complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A',
                      'a':'T', 'c':'G', 'g':'C', 't':'A',
                      'n':'N', 'N':'N'}
    
    # convert string into list with each character as a separate element
    bases = list(dna_sequence)
    
    # determine complement strand
    bases = [base_complement[base] for base in bases]
    
    complement_strand = ''.join(bases)
    
    # reverse strand
    reverse_complement_strand = complement_strand[::-1]
    
    return reverse_complement_strand


def translate(dna_sequence):
    """ Converts a given DNA sequence into a protein sequence.
    
    Args:
        DNASeq (str): a string representing a DNA sequence.
    
    Returns:
        protseq (str): a string representing the protein sequence obtained from 
        the translation of the input DNA sequence.
        trans_state (str): indicates the DNA strand that coded for the protein
        sequence (sense or antisense).
        
    Example:
        >>> translate_dna("GTGACACCAAAACCATGA")
        ['MTPKP', 'sense']
    
    Note:
        The sequence is scanned for ambiguous nucleotides according to IUPAC definitions.
        The translation table used is the 'The Bacterial, Archaeal and Plant Plastid Code'
        from NCBI.
        For more information about the translation table used visit:
            https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    
    """
    
    sequence = dna_sequence.upper()
    tableid = 11
    start_codons = ('ATG','GTG','TTG','ATT','CTG')
    stop_codons = ('TAA','TAG','TGA')
    translated_strand = ''

    # sequences with ambiguous nucleotides are not translated
    if any((nuc in 'RYWSMKHBVDN') for nuc in sequence):
        ambiguous = 'amb'
        
        return ambiguous
    
    # if the sequence has no ambiguous nucleotides
    else:
        # check for start and stop codons
        start = sequence.startswith(start_codons)
        stop = sequence.endswith(stop_codons)
        
        # translate sequence if it is a CDS
        if start and stop:
            #print('sense')
            myseq = Seq(sequence)
            protseq = Seq.translate(myseq, table=tableid, cds=True)
            translated_strand = 'sense'
        
        # reverse complement the strand to check for CDS if sense strand had no CDS
        else:
            #print('going to antisense')
            revC_seq = reverse_complement(sequence)
            start = revC_seq.startswith(start_codons)
            stop = revC_seq.endswith(stop_codons)

            if start and stop:
                #print('antisense')
                myseq = Seq(revC_seq)
                protseq = Seq.translate(myseq, table=tableid, cds=True)
                translated_strand = 'antisense'

    # if previous translation attempts failed, the sense and antisense
    # sequences are reversed to check for the possibility that sequences 
    # provided by users had reverse orientation
    
    if translated_strand == '':
        rev_seq = sequence[::-1]
        start = rev_seq.startswith(start_codons)
        stop = rev_seq.endswith(stop_codons)
        
        if start and stop:
            #print("r-sense")
            myseq = Seq(rev_seq)
            protseq = Seq.translate(myseq, table=tableid, cds=True)
            translated_strand = 'r-sense'
        
        else:
            #print("going to r-antisense")
            revrevC_seq = reverse_complement(rev_seq)
            start = revrevC_seq.startswith(start_codons)
            stop = revrevC_seq.endswith(stop_codons)
            
            if start and stop:
                #print("r-antisense")
                myseq = Seq(revrevC_seq)
                protseq = Seq.translate(myseq, table=tableid, cds=True)
                translated_strand = 'r-antisense'
                
            else:
                non_translatable = 'untrans'
                
                return non_translatable

    return [str(protseq), translated_strand]


# CDSs determined by Prodigal are translatable at the first try or by reverse complementing.
def translate_coding_sequences(coding_sequences):
    """ Translates CDSs into protein sequences.
    
        Args:
            cdss (dict): a dictionary with CDSs ids as keys and CDSs DNA 
            sequences as values.
            
        Returns:
            prots (dict): a dictionary with CDSs/proteins ids as keys and protein
            sequences as values.
            trans_state (dict): a dictionary with the CDSs/proteins ids as keys 
            and the DNA strand that coded for those proteins.
            ambiguous (dict): a dictionary with CDSs/proteins ids as keys and 
            CDSs DNA sequences that had ambiguous bases as values.
            untranslatable (dict): a dictionary with CDSs/proteins ids as keys 
            and CDSs DNA sequences that could not be translated.
            
        Example:
            
            >>> translate_cds(coding_sequences)
            {1: 'MTPKPVPDKDKYDPTG',
            2: 'MSPLGMIKDEGLFNTELD',
            ...},
            {1: 'sense',
            2: 'antisense',
            ...},
            {3: 'AAATTTNTGCATGA',
            4: 'ATGATTTTTGCBTGA',
            ...}
            
    """
    
    proteins = {}
    translated_strand = {}
    ambiguous = {}
    untranslatable = {}
    # for each CDS id and sequence
    for protid, sequence in coding_sequences.items():
        #print(protid)
        # attempt to translate the DNA sequence into protein
        translation = translate(sequence)
        # store the protein with the same id as the CDS and the information
        # about the coding strand
        if type(translation) == list:
            protein = translation[0]
            translated_strand[protid] = translation[1]
            proteins[protid] = protein
        
        # store the CDSs sequences with ambiguous nucleotides that are not translated
        elif translation == 'amb':
            ambiguous[protid] = sequence
        
        # also store the DNA sequences that cannot be translated because no 
        # valid CDS features can be found
        elif translation == 'untrans':
            untranslatable[protid] = sequence
    
    return [proteins, translated_strand, ambiguous, untranslatable]


# find and remove equal proteins in same set
def determine_repeated(proteins):
    """
    """
    
    repeated_proteins = {}
    
    for key, value in proteins.items():
        if value not in repeated_proteins:
            repeated_proteins[value] = [key]
        else:
            repeated_proteins[value].append(key)
    
    repeated_ids = [id_list for id_list in list(repeated_proteins.values()) if len(id_list) > 1]
    
    return repeated_ids


def remove_repeated(proteins, repeated_ids):
    """ Finds repeated proteins and removes them from a dictionary.
    
        Args:
            prots (dict):
            repeated_ids (list)
        Returns:
            prots (dict):
            
        Example:
            
    """
    
    for id_list in repeated_ids:
        to_remove = id_list[1:]
        for protein in to_remove:
            proteins.pop(protein, None)
            
    return proteins


def invert_dictionary(dictionary):
    """
    """
    
    inverted = {value:key for key, value in dictionary.items()}
    
    return inverted


# alter to get DNA sequence instead of protein sequence???
def determine_small(proteins, minimum_length):
    """ Find protein sequences that are shorter than desired length.
    
        Args:
            prots (dict): a dictionary with protein ids as keys and protein 
            sequences as values.
            min_len (int): Proteins with a number of amino acids lower than 
            this value are considered small.
            
        Returns:
            small_proteins (dict): a dictionary with the ids of small proteins 
            as keys and their amino acid sequence as values.
                
        Example:
            
            >>> small_prots(prots, 67)
            {5: 'MTPKPVDKDKYD',
            19: 'MTLNEMVGYVISAHHGMYDFCYCSDDAE',
            ...}
            
    """
    
    small_proteins = {}
    for protid, sequence in proteins.items():
        if len(sequence) < minimum_length:
            small_proteins[protid] = sequence
            
    return small_proteins


# remove small proteins
def remove_small(proteins, small_proteins):
    """ Removes small proteins from the dictionary with all proteins.
    
        Args:
            proteins (dict): a dictionary with protein ids as keys and protein 
            sequences as values.
            small_prots (dict): a dictionary with small proteins ids as keys and 
            proteins sequences for those small proteins as values.
            
        Returns:
            proteins (dict): a dictionary with all proteins except the small 
            proteins.
            
    """

    for protid in small_proteins:
        if protid in proteins:
            del proteins[protid]

    return proteins


def determine_common(protein_group1, protein_group2):
    """
    """
    
    # invert proteins dictionaries
    # proteins become keys and ids become values
    inverted_group1 = invert_dictionary(protein_group1)
    inverted_group2 = invert_dictionary(protein_group2)
    
    # find the intersection of both protein sets
    common_proteins = [protein for protein in inverted_group1 if protein in inverted_group2]
    
    return [inverted_group1, inverted_group2, common_proteins]


def protein_groups_union(protein_group1, protein_group2, common_proteins):
    """ Determines the proteins common to two protein sets and returns the 
        union of both sets without duplicates.
        
        Args:
            prots1, prots2 (dict): dictionaries with proteins ids as keys and 
            proteins sequences as values.
        
        Returns:
            prot_merge (dict): a dictionary with protein sequences as keys and
            protein ids as values. This dictionary results from merging both 
            dictionaries given as input, excluding duplicates.
            
        Example:
            
            
    """
    
    # remove common proteins from one of the sets
    for protein in common_proteins:
        protein_group2.pop(protein, None)
    
    # merge both sets to obtain the unique set of proteins
    merged_groups = {**protein_group1, **protein_group2}
    
    return merged_groups


# check if any protein is contained in another
def determine_contained(proteins):
    """
    """
    
    # get all prots
    protein_sequences = list(proteins.keys())
    # sort by ascending length to facilitate search
    protein_sequences = sorted(protein_sequences, key=len)
    
    # verify if each prot is contained in another, starting with smaller ones
    contained = {}
    for protein in range(len(protein_sequences)):
        current_protein = protein_sequences[protein]
        
        # the set of searched proteins will decrease with each iteration
        # due to initial sorting
        for protein2 in range(protein+1,len(protein_sequences)):
            if current_protein in protein_sequences[protein2]:
                contained[current_protein] = proteins[current_protein]
    
    return contained
    

def remove_contained(proteins, contained):
    """ Finds proteins that are contained in other proteins and removes those 
        proteins from the complete set of proteins.
        
        Args:
            prots (dict): a dictionary with proteins as keys and proteins 
            ids as values.
        
        Returns:
            final_proteins (dict): equivalent to the input dictionary without 
            proteins contained in other proteins.
        
        Example:    
            
    """
                
    # remove proteins contained in others
    for protein in contained:
        proteins.pop(protein, None)
    
    # invert dictionary to get protein ids as keys and protein sequences as values
    final_proteins = invert_dictionary(proteins)
    
    return final_proteins


def create_fasta_lines(sequences, genome_ids):
    """
    
        Args:
            
        Returns:
            
        Example:
    """
    
    lines = []
    for seqid in sequences:
        header = '>' + genome_ids[seqid] + '|protein' + str(seqid)
        sequence = sequences[seqid]
        
        lines.append(header)
        lines.append(sequence)

    return lines


def write_fasta(fasta_lines, output_file):
    """
    """
    
    joined_lines = '\n'.join(fasta_lines)
    
    with open(output_file, 'a') as file:
        file.write(joined_lines)
    
    return 'Wrote FASTA sequences to ' + output_file


def write_protein_table(file_name, *argv):
    """
    """
    
    total_proteins = 0
    with open(file_name, 'a') as file:
        file.write('Genome\tcontig\tStart\tStop\tprotID\n')
        
        for protein_list in argv:
            lines = ['\t'.join(protein)+'\n' for protein in protein_list]
            file.writelines(lines)
            
            total_proteins += len(lines)
    
    return 'Wrote information about ' + str(total_proteins) + ' proteins to ' + file_name


######################################################
# USING FUNCTIONS TO GET TO SAME RESULT AS CHEWBBACA #
######################################################

#genomes_path = '/home/pcerqueira/DATA/chewbbaca/GBS_Aug2016'
##genomes_path = '/home/pcerqueira/Lab_Software/testing/one_fasta_files'
#all_fasta_files = os.listdir(genomes_path)
#for f in range(len(all_fasta_files)):
#    all_fasta_files[f] = os.path.join(genomes_path, all_fasta_files[f])
#
#protogenome_name = 'protogenome0'
##temp_path = '/home/pcerqueira/Lab_Software/testing/temp_gbs'
#
#temp_path = '/home/pcerqueira/Lab_Software/testing/temp'
#
#protogenome_path = os.path.join(temp_path, protogenome_name)
#cpu = 6
## blastp = 
## verbose = 
#bsr = 0.6
#
## get genome/assemblies ids
#genomes_identifiers = []
#for genome in all_fasta_files:
#    genomes_identifiers.append(genome_id(genome))
#
## import the contigs of each assembly
#genomes_contigs = []
#for genome in all_fasta_files:
#    genomes_contigs.append(import_contigs(genome))

#
#
#
#
# Import ORF files with the start and end positions
# extract the CDSs from the contigs, create lines with CDSs info 
#genomes_cds = []
#protid = 1
#for g in range(len(all_fasta_files)):
#    orf_file_path = os.path.join(temp_path, os.path.basename(genomes_identifiers[g]) + '.fna_ORF.txt')
#    genomes_cds.append(extract_coding_sequences(orf_file_path, genomes_contigs[g], protid, genomes_identifiers[g]))
#    protid = genomes_cds[g][4]
#    
#
#
#all_seqs_dict = {}
#for g in genomes_cds:
#    all_seqs_dict = {**all_seqs_dict, **g[0]}
#    
#
#all_genomes_prots = translate_coding_sequences(all_seqs_dict)
#
#print("Contig ID: " + str(genomes_cds[2][16055]))
#
#print("Genome: " + str(genomes_cds[3][16055]))

#
#protid_genome_ids = {}
#for g in genomes_cds:
#    protid_genome_ids = {**protid_genome_ids, **g[3]}
#
## remove proteins in common
#all_genomes_repeats = determine_repeated(all_genomes_prots[0])
## from protein dictionary
#all_genomes_prots[0] = remove_repeated(all_genomes_prots[0], all_genomes_repeats)
## from DNA dictionary
#all_seqs_dict = remove_repeated(all_seqs_dict, all_genomes_repeats)
#
#all_genomes_small = determine_small(all_genomes_prots[0], 67)
## from protein dictionary
#all_genomes_prots[0] = remove_small(all_genomes_prots[0], all_genomes_small)
## from DNA dictionary
#all_seqs_dict = remove_small(all_seqs_dict, all_genomes_small)
#
## invert dictionaries
#iverted_protein_dict = invert_dictionary(all_genomes_prots[0])
##iverted_dna_dict = invert_dictionary(all_seqs_dict)
#
## determine proteins contained in others
#contained = determine_contained(iverted_protein_dict)
#
## remove contained from proteins
#finally_proteins = remove_contained(iverted_protein_dict, contained)
#all_genomes_prots[0] = finally_proteins
#
## remove contained from DNA sequences
#finally_cds = remove_contained(all_seqs_dict, list(contained.values()))
#finally_cds = invert_dictionary(finally_cds)
#
## write output files
#genomes_lines = []
#for genome in genomes_cds:
#    genomes_lines += genome[1]
#
## write protein info table
## add translated strand info to info
#for l in genomes_lines:
#    l.append(all_genomes_prots[1][int(l[-1])])
#
#protein_table_path = os.path.join(os.getcwd(), 'proteinID2_Genome.tsv')
#write_protein_table(protein_table_path, genomes_lines)
#
## write protein FASTA file
#protein_file = os.path.join(protogenome_path, protogenome_name + '_proteins.fasta')
#protein_fasta_lines = create_fasta_lines(all_genomes_prots[0], protid_genome_ids)
#write_fasta(protein_fasta_lines, protein_file)
#
## write DNA FASTA file
#dna_file = os.path.join(protogenome_path, protogenome_name + '.fasta')
#dna_fasta_lines = create_fasta_lines(finally_cds, protid_genome_ids)
#write_fasta(dna_fasta_lines, dna_file)
#
#########################
## CreateSchema.py Part #
#########################
#
#from Bio.Blast.Applications import NcbiblastpCommandline
#from createschema import init_schema_4_bbaca,CommonFastaFunctions
#
#protein_file = '/home/rfm/Lab_Software/chewBBACA/CHEWBBACA/temp/protogenome0/protogenome0_proteins.fasta'
##extended_protein_file = '/home/rfm/Lab_Software/chewBBACA/CHEWBBACA/temp/protogenome0/extended.fasta'
#Gene_Blast_DB_name = CommonFastaFunctions.Create_Blastdb(protein_file, 1, True)
#BlastpPath = '/usr/bin/blastp'
#geneF = os.path.splitext(protein_file)[0]
#blast_out_file = geneF + '.xml'
#bsr = 0.6
#
#cline = NcbiblastpCommandline(cmd=BlastpPath, query=protein_file, 
#                              db=Gene_Blast_DB_name, evalue=0.001,
#                              out=blast_out_file, outfmt=5, num_threads=8)
#
## generator
#
#start_time = time.time()
#
#blast_records = CommonFastaFunctions.runBlastParser(cline, blast_out_file)
#
#end_time = time.time()
#delta = end_time - start_time
#print(delta)
#
## file with CDSs FASTA sequences
#fasta_dict = {}
#for gene in SeqIO.parse('/home/rfm/Lab_Software/chewBBACA/CHEWBBACA/temp/protogenome0/protogenome0.fasta', 'fasta', generic_dna):
#    fasta_dict[gene.id] = str(gene.seq)
#
#
########################
## BSR Function1
#
## Always check if the result with this function is equal to the original chewBBACA
#
#to_remove = []
#to_keep = []
#all_hits = {}
#for blast_record in blast_records:
#
#    allele_name = blast_record.query
#    allele_length = len(fasta_dict[allele_name])
#    try:
#        if allele_name not in to_remove:
#            to_keep.append(allele_name)
#            
#            alignments = blast_record.alignments
#            alignments_dict = {}
#            for hit in alignments:
#                hit_id = hit.hit_def
#                hit_blast_score = hit.hsps[0].score
#                hit_length = len(fasta_dict[hit_id])
#                
#                alignments_dict[hit_id] = [hit_blast_score, hit_length]
#            #all_hits[allele_name] = alignments_dict
#            
#            self_blast_score = alignments_dict[allele_name][0]
#            
#            for hit in alignments_dict:
#                if hit != allele_name:
#                
#                    hit_blast_score = alignments_dict[hit][0]
#                    hit_length = alignments_dict[hit][1]
#                    blast_score_ratio = float(hit_blast_score) / float(self_blast_score)
#                    
#                    if blast_score_ratio > bsr and hit not in to_remove:
#                        all_hits[allele_name] = alignments_dict
#
#                        
#                        if hit_length > allele_length:
#                            to_keep.append(hit)
#                            to_keep.remove(allele_name)
#                            to_remove.append(allele_name)
#                            
#                        elif hit in to_keep:
#                            to_keep.remove(hit)
#                            to_remove.append(hit)
#    except:
#        pass
#
########################
## BSR Function2
#
## This strategy gives different results than original but might be more correct...
#
#to_remove2 = []
#all_hits2 = {}
#for blast_record in blast_records:
#
#    allele_name = blast_record.query
#    allele_length = len(fasta_dict[allele_name])
#
#    if allele_name not in to_remove2:
#        
#        alignments = blast_record.alignments
#        alignments_dict = {}
#        for hit in alignments:
#            hit_id = hit.hit_def
#            hit_blast_score = hit.hsps[0].score
#            hit_length = len(fasta_dict[hit_id])
#            
#            alignments_dict[hit_id] = [hit_blast_score, hit_length]
#        #all_hits2[allele_name] = alignments_dict
#
#        self_blast_score = alignments_dict[allele_name][0]
#        
#        for hit in alignments_dict:
#            if hit != allele_name:
#            
#                hit_blast_score = alignments_dict[hit][0]
#                hit_length = alignments_dict[hit][1]
#                blast_score_ratio = float(hit_blast_score) / float(self_blast_score)
#                
#                if blast_score_ratio > bsr and hit not in to_remove2:
#                    all_hits2[allele_name] = alignments_dict
#                    
#                    if hit_length > allele_length and allele_name not in to_remove2:
#                        to_remove2.append(allele_name)
#                    
#                    # remove hits smaller or of equal length than the current allele
#                    #elif hit_length < allele_length:
#                    elif hit_length <= allele_length:
#                        to_remove2.append(hit)

# The CreateSchema process creates FASTA files for each unique gene that only 
# have one sequence/allele? The AlleleCall process is the one that adds new Alleles???










