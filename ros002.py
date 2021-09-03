#!/usr/bin/env python3
import sys


def translate_dna_to_rna(dna: str):
    rna = ''
    for cha in dna:
        if cha.upper() == 'T':
            rna += 'U'
        else:
            rna += cha.upper()

    return rna

def agtc_content(dna):
    a = 0
    c = 0
    t = 0
    g = 0
    for cha in dna.upper():
        if cha == "A":
            a += 1
        elif cha == "C":
            c += 1
        elif cha == "T":
            t += 1
        elif cha == "G":
            g += 1
    print(a, c, t, g)

def count_chars(string):
    char_count = {}
    for cha in string.upper():
        if cha in char_count.keys():
            char_count[cha] += 1
        else:
            char_count[cha] = 1
    return char_count

def reverse_complement(dna_string):
    result = ""
    for cha in dna_string.upper()[::-1]:
        if cha == 'A':
            result += 'T'
        elif cha == 'T':
            result += 'A'
        elif cha == 'C':
            result += 'G'
        elif cha == 'G':
            result += 'C'

    return result


def fibonacci_rabbits (n, k):
    big_rabbits = 0
    small_rabbits = 1

    for m in range(n-1):
        prev_big_rabbits = big_rabbits
        big_rabbits = small_rabbits + big_rabbits
        small_rabbits = prev_big_rabbits * k

    return big_rabbits + small_rabbits


def mortal_fibonacci_rabbits(n, m):
    rabbits = {age: 1 if age == 0 else 0 for age in range(m)}

    for k in range(n-1):
        babies = sum([rabbits[age] if age > 0 else 0 for age in range(m)])
        for age in range(1, m)[::-1]:
            rabbits[age] = rabbits[age-1]
        rabbits[0] = babies

    '''
    for age, number_of_rabbits in rabbits.items():
        print(f'there are {number_of_rabbits} aged {age} months')
    '''

    return sum(rabbits.values())

if __name__ == '__main__':

    a= 'GCGACTCGATATCGACTGTATCAAAGCTCACGAGCCTTGGTGCTGTAGTCTACGGCGAGGCTCCTTCCTCCTTTAGGCGACGCTCCTAATTTTGAAATGGTACAAAATCACGCAACCGGCTCACGGCGGTACCAAATTAAGAGGACGAAGACCTGATGAGCAACCACATGTACAGTGCTGTCACGATGCTTAAGTCATCCGAGAAAAGCAACAGGCCCTCCATTCTAAAGGCATTTCGGCAAGGGCCGCGTGTACTCGTTTCGGTCACCGCGTGTCCAAGGCTAAAAGCGTACGGGCACGAAGGTCAATCGCGTGGGGCGTACTCCGCTACGTTGTTCGGACAGAGAAGGCACACATCAGGGTCATGCACCCACGCGGTTCATCTTGTGCATTCGGTTTTCATTTCGTGTCGTGACTGCCGCGATGCGATAATATATTTTCTTATACATTAGAAGTACTGTAATTTGAGTGTCAACAGTAATAAGGATGAGCGACGAGCCGAAGACGATAGACTTAGGAAGGTCTTTTCAATCCTACAACCTTAGGAACTGTATCTTAGGTGCAGGGGCGCCCTTTGGTTTTACCGTCAGAGACATCGCCTAAGTCATTCGTCATGGCCGCGTACGCACGATACTATTTGACCAACCCGTTGACCTTTCCGGTCCTTTGCCTTCACGGCGGTTATACGGGATACAAGTCCCGAATTGGTATGCAACGATGCTAACCACGTGAACTGTGCAAGTACCCAGATGGGTGGTACTTACTGTATACTCCACCTACTCGGTCCGGGTGAGAAGTCAGCTCCCTGCCCATATACAGCCCAGGTTAGCGCGCCCTAGAATGGGCGCGACTCTACATTTCGGGGTTGTGGGTCGGAAACAAAAATGACCTCCGTTACCGGAATAACTCTAAGACCGGAC'
    print(mortal_fibonacci_rabbits(96, 16))



