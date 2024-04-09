# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 16:13:57 2024

@author: Pauline
"""
from heapq import heapify, heappop, heappush
import numpy as np
from Crypto.Util.number import *

def DNA_to_AA(DNA_sequence):
    genetic_code = { 
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}  # * is the STOP codon
    
    len_sequence = len(DNA_sequence)
    AA_sequence = ""

    for i in range(0, len_sequence, 3):
        codon = DNA_sequence[i:i+3]
        if len(codon) == 3:
            AA_sequence += genetic_code.get(codon)

    return AA_sequence


#Read the fna file and get the DNA sequence
DNA_file = "DNA_CFTR.txt"
with open(DNA_file, 'r') as file : 
    lines = file.readlines() 
lines = lines[1:]#1st line contains additionnal information, DNA sequence begins at the 2nd line
DNA = ''.join(lines).replace("\n", "") #concatenate all the lines and remove the \n

DNA_seq = list(DNA) 
print("The DNA sequence:", DNA_seq, "\n")

codon_seq = [DNA[i:i+3] for i in range(0, len(DNA_seq), 3)]
print("The codon sequence:", codon_seq, "\n")

#Translate the DNA sequence into the AA sequence
AA_seq = list(DNA_to_AA(DNA))
print("The amino-acid sequence:", AA_seq, "\n")



def count_occurrence(sequence):
    occurrences = {}
    for i in range(0, len(sequence)):
        character = sequence[i]

        if character in occurrences:
            occurrences[character] += 1
        else:
            occurrences[character] = 1

    return occurrences


occ_DNA = count_occurrence(DNA_seq)
print("Occurences in the DNA sequence : ",occ_DNA , "\n")

occ_codon = count_occurrence(codon_seq)
print("Occurences in the DNA sequence (by codons) : ",occ_codon, "\n")

occ_AA = count_occurrence(AA_seq)
print("Occurences in the amino-acid sequence : ",occ_AA, "\n")

def print_tree(node, depth=0, prefix=''):
    if node is not None:
        print(' ' * depth * 4 + prefix + str(node.ch) + ':' + str(node.freq))
        print_tree(node.left, depth + 1, 'L-- ')
        print_tree(node.right, depth + 1, 'R-- ')

class Node:
    def __init__(self, ch, freq, left=None, right=None):
        self.ch = ch
        self.left = left
        self.right = right
        self.freq = freq
    def __lt__(self, other):
        return self.freq < other.freq

def visit_node(node, path="", codes=None):
    if codes is None:
        codes = {}

    # If the node is not a leaf node, recursively visit its children
    if node.left is not None:
        visit_node(node.left, path=path + "1", codes=codes)
    if node.right is not None:
        visit_node(node.right, path=path + "0", codes=codes)

    # If the node is a leaf node, store its code in the dictionary
    if node.left is None and node.right is None:
        codes[node.ch] = path

    return codes


def Huffman(occ_seq):
    s = sum(occ_seq.values(),0.0)
    pq = [Node(k, v/s) for k,v in occ_seq.items()]
    heapify(pq)
    while len(pq) > 1:
        left, right = heappop(pq), heappop(pq)
        newFreq = left.freq + right.freq
        heappush(pq, Node(None, newFreq, left, right))
    
    root = pq[0]

    return root

def get_avg(code):
    total_symbols = 0
    total_len = 0
    for k,v in code.items():
        total_symbols += 1
        total_len += len(v)
    print("Average bits per symbol: ", total_len / total_symbols)

def get_entropy(seq):
    entropy = 0.0
    total = sum(seq.values(),0)
    for k,v in seq.items():
        v = v / total
        entropy += - np.log2(v)*v
    print("Entropy per symbol: ", entropy)

def to_binary(seq, Huff_seq):
    encoded_seq = ""
    for k in seq:
        code = Huff_seq[k]
        encoded_seq+=code
    return encoded_seq



# DNA SEQUENCE
root = Huffman(occ_DNA)
#print_tree(root)

codes = visit_node(root)
print(codes)
get_entropy(occ_DNA)
get_avg(codes)

seq = to_binary(DNA_seq, codes)
with open("DNA.bin", "wb") as f:
    f.write(long_to_bytes(int(seq,2)))



# CODON SEQUENCE
root = Huffman(occ_codon)
#print_tree(root)

codes = visit_node(root)
print(codes)
get_entropy(occ_codon)
get_avg(codes)

seq = to_binary(codon_seq, codes)
with open("codon.bin", "wb") as f:
    f.write(long_to_bytes(int(seq,2)))


# AA SEQUENCE
root = Huffman(occ_AA)
#print_tree(root)

codes = visit_node(root)
print(codes)
get_entropy(occ_AA)
get_avg(codes)

seq = to_binary(AA_seq, codes)
with open("AA.bin", "wb") as f:
    f.write(long_to_bytes(int(seq,2)))



  




