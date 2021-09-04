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


def fibonacci_rabbits(n, k):
    big_rabbits = 0
    small_rabbits = 1

    for m in range(n - 1):
        prev_big_rabbits = big_rabbits
        big_rabbits = small_rabbits + big_rabbits
        small_rabbits = prev_big_rabbits * k

    return big_rabbits + small_rabbits


def mortal_fibonacci_rabbits(n, m):
    rabbits = {age: 1 if age == 0 else 0 for age in range(m)}

    for k in range(n - 1):
        babies = sum([rabbits[age] if age > 0 else 0 for age in range(m)])
        for age in range(1, m)[::-1]:
            rabbits[age] = rabbits[age - 1]
        rabbits[0] = babies

    '''
    for age, number_of_rabbits in rabbits.items():
        print(f'there are {number_of_rabbits} aged {age} months')
    '''

    return sum(rabbits.values())


def parse_fasta(path):
    with open(path, 'r') as file:
        whole = [line.strip() for line in file]

    grouped_dna = ''
    label = ''
    result = {}
    for element in whole:
        if element.startswith('>'):
            if label != '':
                result[label] = grouped_dna
            label = element
            grouped_dna = ''
        else:
            grouped_dna += element
    return result


def gc_content(path):
    result = {}
    fasta = parse_fasta(path)
    for label, dna in fasta.items():
        char_counts = count_chars(dna)
        gc_cont = (char_counts['G'] + char_counts['C']) / len(dna) * 100
        result[label] = gc_cont
    return result

def find_key_of_max_value_in_dict(dictionary):
    for key, value in dictionary.items():
        if value == max(dictionary.values()):
            return key, value


def hamming_distance(path):
    with open(path, 'r') as file:
       [str1, str2] = [line.strip() for line in file]

    hamming_count = 0
    for (a, b) in zip(str1, str2):
        if a != b:
            hamming_count += 1
    return hamming_count

def remove_list(list, a):
    return list.remove(a)

def enumerating_gene_orders(n):
    sequence = list(range(n))
    subset = []
    rest = sequence[:]
    for a in sequence:
        subset.append(a)
        rest.remove(a)


def generate_permutations(sequence: list):
    permutations = []
    for element in sequence:
        permutation = []
        permutation.append(element)
        rest = sequence[:]
        rest.remove(element)
        permutations.append((permutation, rest))
    print(permutations)

    pass
codon_table = '''UUU F      CUU L      AUU I      GUU V
UUC F      CUC L      AUC I      GUC V
UUA L      CUA L      AUA I      GUA V
UUG L      CUG L      AUG M      GUG V
UCU S      CCU P      ACU T      GCU A
UCC S      CCC P      ACC T      GCC A
UCA S      CCA P      ACA T      GCA A
UCG S      CCG P      ACG T      GCG A
UAU Y      CAU H      AAU N      GAU D
UAC Y      CAC H      AAC N      GAC D
UAA Stop   CAA Q      AAA K      GAA E
UAG Stop   CAG Q      AAG K      GAG E
UGU C      CGU R      AGU S      GGU G
UGC C      CGC R      AGC S      GGC G
UGA Stop   CGA R      AGA R      GGA G
UGG W      CGG R      AGG R      GGG G'''

def parse_codon_string(codon_table):
    codon_dict = {}
    codon_list = codon_table.replace('\n', ' ').split(' ')
    for i in range(len(codon_list)):
        if len(codon_list[i]) == 3:
            key = codon_list[i]
            value = codon_list[i+1]
            codon_dict[key] = value
    return codon_dict


CODON_DICT = {'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V', 'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V',
              'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V', 'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V',
              'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A', 'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
              'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A', 'UCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
              'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D', 'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
              'UAA': '_', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E', 'UAG': '_', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
              'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G', 'UGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
              'UGA': '_', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G', 'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}

def translating_rna_to_protein(path, cut_after_stop_codon=True):
    with open(path, 'r') as file:
        rna = ''
        for line in file:
            rna += line.strip()

    codons = [rna[n:n+3] for n in range(0, len(rna), 3)]
    result = ""
    for codon in codons:
        result += CODON_DICT[codon]
    if cut_after_stop_codon:
        result, _ = result.split('_')
    return result


if __name__ == '__main__':
    # gc_contents = gc_content('/Users/testtest/Downloads/rosalind_gc.txt')
    #label, percentage = find_key_of_max_value_in_dict(gc_contents)

    print(translating_rna_to_protein('/Users/testtest/Downloads/rosalind_prot.txt'))
    #print(translating_rna_to_protein())