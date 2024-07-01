#! /usr/bin/env python

import pandas as pd
import numpy as np
import math
import time
import glob
import json
import os
import requests
import re

# create a dictionary of each AA and corresponding codon
aa_dict = {'A':["GCT", "GCC", "GCA", "GCG"],
           'C':["TGT", "TGC"],
           'D':["GAT", "GAC"],
           'E':["GAA", "GAG"],
           'F':["TTT", "TTC"],
           'G':["GGT", "GGC", "GGA", "GGG"],
           'H':["CAT", "CAC"],
           'I':["ATT", "ATC", "ATA"],
           'K':["AAA", "AAG"],
           'L':["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"],
           'M':["ATG"],
           'N':["AAT", "AAC"],
           'P':["CCT", "CCC", "CCA", "CCG"],
           'Q':["CAA", "CAG"],
           'R':["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
           'S':["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
           'T':["ACT", "ACC", "ACA", "ACG"],
           'V':["GTT", "GTC", "GTA", "GTG"],
           'W':["TGG"],
           'Y':["TAT", "TAC"]}
# dictionary that doesn't count stops and single codon AA
aa_dict_sub = {'A':["GCT", "GCC", "GCA", "GCG"],
           'C':["TGT", "TGC"],
           'D':["GAT", "GAC"],
           'E':["GAA", "GAG"],
           'F':["TTT", "TTC"],
           'G':["GGT", "GGC", "GGA", "GGG"],
           'H':["CAT", "CAC"],
           'I':["ATT", "ATC", "ATA"],
           'K':["AAA", "AAG"],
           'L':["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"],
           'N':["AAT", "AAC"],
           'P':["CCT", "CCC", "CCA", "CCG"],
           'Q':["CAA", "CAG"],
           'R':["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
           'S':["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
           'T':["ACT", "ACC", "ACA", "ACG"],
           'V':["GTT", "GTC", "GTA", "GTG"],
           'Y':["TAT", "TAC"]}
# list of codons
codon_list = ['TTT','TTC','TTA','TTG','CTT','CTC',
                       'CTA','CTG','ATT','ATC','ATA','ATG',
                       'GTT','GTC','GTA','GTG','TAT','TAC',
                       'TAA','TAG','CAT','CAC','CAA','CAG',
                       'AAT','AAC','AAA','AAG','GAT','GAC',
                       'GAA','GAG','TCT','TCC','TCA','TCG',
                       'CCT','CCC','CCA','CCG','ACT','ACC',
                       'ACA','ACG','GCT','GCC','GCA','GCG',
                       'TGT','TGC','TGA','TGG','CGT','CGC',
                       'CGA','CGG','AGT','AGC','AGA','AGG',
                       'GGT','GGC','GGA','GGG']
# list of start stop and single
codons_to_exclude = ["TGG", "ATG", "TAG", "TAA", "TGA"]
# list of codons fro aa_dict_sub
codons_to_use = ['TTT','TTC','TTA','TTG','CTT','CTC',
               'CTA','CTG','ATT','ATC','ATA',
               'GTT','GTC','GTA','GTG','TAT','TAC',
               'CAT','CAC','CAA','CAG',
               'AAT','AAC','AAA','AAG','GAT','GAC',
               'GAA','GAG','TCT','TCC','TCA','TCG',
               'CCT','CCC','CCA','CCG','ACT','ACC',
               'ACA','ACG','GCT','GCC','GCA','GCG',
               'TGT','TGC','CGT','CGC',
               'CGA','CGG','AGT','AGC','AGA','AGG',
               'GGT','GGC','GGA','GGG']
# nucleotide list
nuc_list = ["A", "T", "C", "G"]

# First step is to read in file and get codon frequencies
#These are the codon counts
input_file = "genes_codon_frequency.csv"
codon_counts = pd.read_csv(input_file)
# Now calculate frequencies for all. Start by summing
codon_sums = codon_counts.sum()
# Get annotations
annotations = pd.read_csv("annotations.tsv", delimiter="\t")
# gather ribosomal genes from kegg
ribosomal_gene_names = pd.read_csv("https://rest.kegg.jp/link/ko/map03010/",
                                   delimiter = "\t", header=None, names = ["path","KO"])
ribosomal_gene_names["KO"] = ribosomal_gene_names["KO"].str.extract("(K[0-9]+)")
# isolate those genes
ribo_annoations = annotations[annotations["kegg_id"].isin(ribosomal_gene_names["KO"])]
# get the sum of codons from all ribosomal genes
rib_codons = codon_counts[codon_counts["ID"].isin(ribo_annoations['Unnamed: 0'])]
ribo_sums = rib_codons.sum()

# Get nucleotide frequenceies for ENC'
input_file_nuc = "genes_nucl_counts.csv"
nuc_counts = pd.read_csv(input_file_nuc)
# Now calculate frequencies for all. Start by summing
nuc_counts = nuc_counts
#Sum them all
nuc_counts_sum = nuc_counts[["A", "T", "C", "G"]].sum()
#While we're at it, turn these counts into frequencies for each gene and save that file
# I should do this when I originally parse them, but this will be our little secret ;)
nuc_counts["total_nuc"] = nuc_counts["A"]+nuc_counts["T"]+nuc_counts["C"]+nuc_counts["G"]
nuc_counts["A"] = nuc_counts["A"]/nuc_counts["total_nuc"]
nuc_counts["T"] = nuc_counts["T"]/nuc_counts["total_nuc"]
nuc_counts["C"] = nuc_counts["C"]/nuc_counts["total_nuc"]
nuc_counts["G"] = nuc_counts["G"]/nuc_counts["total_nuc"]
nuc_counts.to_csv("genes_nucl_frequency.csv")
# Nucleotide frequencies for each
nuc_counts_sum = nuc_counts_sum/nuc_counts_sum.sum()
nuc_counts_sum

# I think the easiest way to do this is to calculate all the expected frequencies up front
output_nuc_ei = dict()
# This loops through each possible combination and uses the background nucleotide frequencies to calculate expected codon values
# These are then added to the dictionary
for first in nuc_list:
    for second in nuc_list:
        for third in nuc_list:
            codon_ei_name = first+second+third
            codon_ei = nuc_counts_sum[first]*nuc_counts_sum[second]*nuc_counts_sum[third]
            output_nuc_ei[codon_ei_name] = codon_ei

# Now to get these as relative synonymous codon usages
## (normalize so the sum of frequencies for each codon in an AA is 1)
rscu_ei = dict()
for AA in aa_dict:
    codons_2_get = aa_dict[AA]
    aa_total_codons = float(0)
    for codon_to_sum in codons_2_get:
        aa_total_codons = aa_total_codons+output_nuc_ei[codon_to_sum]
    for codon_to_sum in codons_2_get:
        rscu_ei[codon_to_sum] = output_nuc_ei[codon_to_sum]/aa_total_codons

#Then go through each amino acid and calculate the relative frequency
# Writing a function for this.
# The input will be the codon counts (codon_sums_in), and an amino acid dictionary
# Keeping the aa_dict local will allow us to add in alternative coding schemes if desired
def get_codon_frequencies(aa_dict_in, codon_sums_in):
    #create dictionaries for both the codon frequencies and cai frequencies
    codon_frequencies = dict()
    cai_frequencies = dict()
    # loop through each amino acid
    for aa in aa_dict_in:
        # get the corresponding codons
        codons_to_calc = aa_dict_in[aa]
        codons_to_calc_cts = codon_sums_in[codons_to_calc]
        # calculate total codons for each amino acid
        sum_codon_cts = codons_to_calc_cts.sum()
        # get the relative frequency of each
        freq_codon_cts = codons_to_calc_cts/sum_codon_cts
        # get the cai frequency
        cai_freq = freq_codon_cts/freq_codon_cts.max()
        # append to dictionary
        for co in codons_to_calc:
            codon_frequencies[co] = freq_codon_cts[co]
            cai_frequencies[co] = cai_freq[co]
    return codon_frequencies, cai_frequencies



def calculate_MILC(input_freq, base_freq):
    Msum = 0
    aa_included = []
    L = 0
    for aa in aa_dict:
        Ma = 0
        codons_to_get = aa_dict[aa]
        sum_codons = float(sum(input_freq[codons_to_get]))
        if sum_codons == 0:
            next
        else:
            aa_included.append(aa)
            for codon in codons_to_get:
                Mcodon = 0
                Oc = float(input_freq[codon])
                L = L + Oc
                if Oc == 0:
                    next
                else:
                    fc = Oc/sum_codons
                    #print("fc "+str(fc))
                    gc = float(base_freq[codon])
                    #print("gc "+str(gc))
                    #print("Oc " +str(Oc))
                    Mcodon = Oc * np.log(fc/gc)
                    Ma = Ma + Mcodon
            Msum = Msum + Ma
    rsum = 0
    for aaM in aa_included:
        ra = len(aa_dict[aaM])
        rsum = rsum + (ra-1)
    Cor = (rsum/L)-0.5
    MILC = (Msum/L)-Cor
    return MILC


def get_CAI_and_FOP(input_freq, expected_cai_freqs):
    cai = float("NaN")
    fop = float("NaN")
    # total # of optimized codons
    fop_optimized_codon = float()
    # total # of codons (excluding start, stop, and single-degen)
    total_codon_ct = int()
    # list of all cai values
    total_cai_adjusted_ct = list()
    for codon in codons_to_use:
        # Blank a list for cai values for codon, blank codon abundance
        cai_adjusted_count = list()
        codon_abundance = float()
        # get abundance of codon
        codon_abundance = input_freq[codon]
        # add to total codon ct
        total_codon_ct = total_codon_ct + codon_abundance
        # get the cai frequency (expected) of that codon
        codon_cai_freq = expected_cai_freqs[codon]
        # make a list thats the length of the codon count for the cai expected frequency
        if (codon_abundance != 0) and (codon_cai_freq != 0):
            cai_adjusted_count = [codon_cai_freq] * codon_abundance
            # add that list to the running total list
            total_cai_adjusted_ct = total_cai_adjusted_ct + cai_adjusted_count
        else:
            pass
        # if the cai frequency is == 1, that means it's the optimal codon, so count it for the FOP value
        if codon_cai_freq == 1:
            fop_optimized_codon = fop_optimized_codon + codon_abundance
    fop = fop_optimized_codon/total_codon_ct
    cai = math.exp(np.log(total_cai_adjusted_ct).mean())
    if cai == 0:
        #print("""
        #Error: CAI equals zero, indicated underflow problem
        # - Setting to NA
        #""")
        CAI = float("NaN")
    else:
        next
    return fop, cai

def get_ENC_values(input_counts):
    F2 = float(0)
    F3 = float(0)
    F4 = float(0)
    F6 = float(0)
    for AA in aa_dict:
        x2 = float(0)
        codon_to_get = aa_dict[AA]
        N_a = input_counts[codon_to_get].sum()
        if N_a > 1:
            k = len(codon_to_get)
            codon_total_aa = float()
            for codon in codon_to_get:
                pi = input_counts[codon]/N_a
                ei = rscu_ei[codon]
                x2 = x2 + (N_a*((pi-ei)**2))/ei
            Fa = (x2 + N_a - k)/(k*(N_a-1))
            if k == 1:
                next
            if k == 2:
                F2 = F2+Fa
            if k == 3:
                F3 = F3+Fa
            if k == 4:
                F4 = F4+Fa
            if k ==6:
                F6 = F6+Fa
    F2 = F2/9
    F4 = F4/5
    F6 = F6/3
    # If there is no 3n degeneracy, average 2n and 4n
    if F3 == 0:
        F3 = (F2 + F4)/2
    else:
        F3 = F3/1
    # If we have to divide by 0, it means there are very few codons. Add an NA
    try:
        ENC = 2 + (9/F2) + (1/F3) + (5/F4) + (3/F6)
    except:
        ENC = float('NaN')
    return ENC

# Calculate relative and CAI requencies for all genes and ribosomal genes
DNA_freqs, DNA_cai_freqs = get_codon_frequencies(aa_dict, codon_sums)
RIBO_freqs, RIBO_cai_freqs = get_codon_frequencies(aa_dict, ribo_sums)
ENCall = get_ENC_values(codon_sums)
ENCribo = get_ENC_values(ribo_sums)
# Calculate delta (as decribed https://doi.org/10.1371/journal.pgen.1000808)
delta_ENC = (ENCall-ENCribo)/ENCall

output_df = pd.DataFrame(columns=["ID", "ENC", "fop", "cai", "milc", "fop_ribo", "cai_ribo", "milc_ribo", "length_bp"])
run_ribo_cai = "yes"
for aa in aa_dict:
    testAA = 0
    codons_to_test = aa_dict[aa]
    for codon in codons_to_test:
        testAA = testAA + RIBO_cai_freqs[codon]
    #print(testAA)
    if testAA == 0:
        run_ribo_cai = "no"
        testAA = 0
    else:
        testAA = 0

for row in range(len(codon_counts["ID"])):
    gene = codon_counts.loc[row]
    enc_to_add = float("NaN")
    fop_to_add = float("NaN")
    cai_to_add = float("NaN")
    Rfop_to_add = float("NaN")
    Rcai_to_add = float("NaN")
    milc = float("NaN")
    Rmilc = float("NaN")
    length = float("NaN")
    enc_to_add = get_ENC_values(gene)
    fop_to_add, cai_to_add = get_CAI_and_FOP(gene, DNA_cai_freqs)
    if run_ribo_cai == "yes":
        Rfop_to_add, Rcai_to_add = get_CAI_and_FOP(gene, RIBO_cai_freqs)
    else:
        pass
    milc = calculate_MILC(gene, DNA_freqs)
    try:
        Rmilc = calculate_MILC(gene, RIBO_freqs)
    except:
        pass
    length = gene[codon_list].sum()*3
    output_df.loc[row] = [gene["ID"], get_ENC_values(gene), fop_to_add, cai_to_add, milc, Rfop_to_add, Rcai_to_add, Rmilc, length]

# get MAG name from working directory. Goes from the end of dir string and gets last dir
wd_str = os.getcwd()
mag = str(re.search("([^/]+$)", wd_str)[1])
output_df.to_csv(mag+"_optimization_output.txt", sep ='\t')
number_of_ribo_genes = len(rib_codons["ID"])
output_list = [mag, ENCall, ENCribo, delta_ENC, number_of_ribo_genes]
enc_output_name = mag + "_ENC_output.txt"
with open(enc_output_name, 'w+') as f:
    for item in output_list:
        f.write("%s\t" % item)
f.close()
