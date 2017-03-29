# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 19:23:25 2017

@author: Devin
"""


import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
string_nucleotides="ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
def translate_function(string_nucleotides): 
    mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
    #mito_table = CodonTable.unambiguous_dna_by_id[2]
    stop_codons= mito_table.stop_codons
    aa_seq_string=""
    isStopCodon=None
    for i in range(0, len(string_nucleotides),3):
        codon=string_nucleotides[i]+string_nucleotides[i+1]+string_nucleotides[i+2]
        for x in stop_codons:
            if codon==x:
                isStopCodon=True
                break
            elif codon!=x:
                isStopCodon=False
                
        if isStopCodon==True:
            break
        elif isStopCodon==False:
            amino_acid=mito_table.forward_table[codon]
            aa_seq_string=aa_seq_string+amino_acid
        
    return(aa_seq_string)

#print(translate_function("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"))

