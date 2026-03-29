# -*- coding: utf-8 -*-
"""
Created on March 18 2016 



"""



from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq

import numpy as np
from matplotlib import pyplot as plt

big_list = []

for n in SeqIO.parse("first.ffn", "fasta"):
    big_list.append(str(n.seq))


genome=''.join(big_list)
print(len(genome))

codons_list = [genome[i:i+3] for i in range(0, len(genome), 3)]    
print(len(codons_list))
# print(codons_list)

TTT=codons_list.count('ttt')  # (Phe/F) Phenylalanine
TTC=codons_list.count('ttc')

TTA=codons_list.count('tta')  # (Leu/L) Leucine
TTG=codons_list.count('ttg')
CTT=codons_list.count('ctt')
CTC=codons_list.count('ctc')
CTA=codons_list.count('cta')
CTG=codons_list.count('ctg')

ATT=codons_list.count('att') # (Ile/I) Isoleucine
ATC=codons_list.count('atc')
ATA=codons_list.count('ata')

ATG=codons_list.count('atg') # (Met/M) Methionine

GTT=codons_list.count('gtt') # (Val/V) Valine
GTC=codons_list.count('gtc')
GTA=codons_list.count('gta')
GTG=codons_list.count('gtg')

TCT=codons_list.count('tct') # (Ser/S) Serine
TCC=codons_list.count('tcc')
TCA=codons_list.count('tca')
TCG=codons_list.count('tcg')
AGT=codons_list.count('agt') # (Ser/S) Serine
AGC=codons_list.count('agc')

CCT=codons_list.count('cct') # (Pro/P) Proline
CCC=codons_list.count('ccc')
CCA=codons_list.count('cca')
CCG=codons_list.count('ccg')

ACT=codons_list.count('act') # (Thr/T) Threonine
ACC=codons_list.count('acc')
ACA=codons_list.count('aca')
ACG=codons_list.count('acg')

GCT=codons_list.count('gct') # (Ala/A) Alanine
GCC=codons_list.count('gcc')
GCA=codons_list.count('gca')
GCG=codons_list.count('gcg')

TAT=codons_list.count('tat') # (Tyr/Y) Tyrosine
TAC=codons_list.count('tac')

TAA=codons_list.count('taa') # Stop
TAG=codons_list.count('tag')
TGA=codons_list.count('tga')

CAT=codons_list.count('cat') # (His/H) Histidine
CAC=codons_list.count('cac')

CAA=codons_list.count('caa') # (Gln/Q) Glutamine
CAG=codons_list.count('cag')

AAT=codons_list.count('aat') # (Asn/N) Asparagine
AAC=codons_list.count('aac')

AAA=codons_list.count('aaa') # (Lys/K) Lysine
AAG=codons_list.count('aag')

GAT=codons_list.count('gat') # (Asp/D) Aspartic acid
GAC=codons_list.count('gac')

GAA=codons_list.count('gaa') # (Glu/E) Glutamic acid
GAG=codons_list.count('gag')

TGT=codons_list.count('tgt') # (Cys/C) Cysteine
TGC=codons_list.count('tgc')

TGG=codons_list.count('tgg') # (Trp/W) Tryptophan  

CGT=codons_list.count('cgt') # (Arg/R) Arginine
CGC=codons_list.count('cgc')
CGA=codons_list.count('cga')
CGG=codons_list.count('cgg')
AGA=codons_list.count('aga') # (Arg/R) Arginine
AGG=codons_list.count('agg')

# AGT=codons_list.count('agt') # (Ser/S) Serine
# AGC=codons_list.count('agc')

# AGA=codons_list.count('aga') # (Arg/R) Arginine
# AGG=codons_list.count('agg')

GGT=codons_list.count('ggt') # (Gly/G) Glycine
GGC=codons_list.count('ggc')
GGA=codons_list.count('gga')
GGG=codons_list.count('ggg')

print('ATG', ATG, 'TGG', TGG)

print('ATG%', ATG/len(codons_list))
print('TGG%', TGG/len(codons_list))



