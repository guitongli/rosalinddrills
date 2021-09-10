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
        fasta_label_list = [line.strip() for line in file]

    grouped_dna = ''
    label = ''
    result = {}
    for element in fasta_label_list:
        if element.startswith('>'):
            if label != '':
                result[label] = grouped_dna
            label = element[1:]
            grouped_dna = ''
        else:
            grouped_dna += element
    result[label] = grouped_dna
    return result


def gc_content(fasta):
    result = {}
    for label, dna in fasta.items():
        char_counts = count_chars(dna)
        gc_cont = (char_counts['G'] + char_counts['C']) / len(dna) * 100
        result[label] = gc_cont
    return result


def find_key_of_max_value_in_dict(dictionary):
    for key, value in dictionary.items():
        if value == max(dictionary.values()):
            return key, value


def hamming_distance(str1, str2):
    distance = 0
    for (a, b) in zip(str1, str2):
        if a != b:
            distance += 1
    return distance

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

def translating_rna_to_protein(rna, cut_after_stop_codon=True):
    codons = [rna[n:n+3] for n in range(0, len(rna), 3)]
    protein = ""
    for codon in codons:
        protein += CODON_DICT[codon]
    if cut_after_stop_codon:
        protein, _ = protein.split('_')
    return protein

def mendels_first_law(k, m, n):
    all_possibilities = (k+m+n)*(k+m+n-1)
    displaying_possibilities = k*(k-1) + 2*k*m + 2*k*n + m*(m-1)*0.75 + 2*m*n*0.5
    return displaying_possibilities/all_possibilities

def locating_motif(s, t):
    positions = ''
    for i in range(len(s)):
        if s[(i-1):(i+len(t)-1)] == t:
           positions += str(i) + ' '
    return positions

def finding_consensus(fasta):
    matrix = [value for key, value in fasta.items()]

    consensus = ''

    profile = {
        'A': [0 for a in range(len(matrix[0]))],
        'C': [0 for c in range(len(matrix[0]))],
        'G': [0 for g in range(len(matrix[0]))],
        'T': [0 for t in range(len(matrix[0]))],
    }
    for dna in matrix:
        for i in range(len(dna)):
            profile[dna[i]][i] += 1
    line_a = 'A: '
    line_c = 'C: '
    line_g = 'G: '
    line_t = 'T: '
    cha_sum = profile['A']
    for i in range(len(cha_sum)):
        consensus_cha = max([profile['A'][i], profile['C'][i], profile['G'][i], profile['T'][i]])
        for char in ['A', 'C', 'G', 'T']:
            if profile[char][i] == consensus_cha:
                consensus += char
                break
    for cha, cha_sum in profile.items():
        if cha == 'A':
            for sub_sum in cha_sum:
                line_a = line_a + str(sub_sum) + ' '
        elif cha == 'C':
            for sub_sum in cha_sum:
                line_c = line_c + str(sub_sum) + ' '
        elif cha == 'G':
            for sub_sum in cha_sum:
                line_g = line_g + str(sub_sum) + ' '
        elif cha == 'T':
            for sub_sum in cha_sum:
                line_t = line_t + str(sub_sum) + ' '

    return consensus + '\n' + line_a + '\n' + line_c + '\n' + line_g + '\n' + line_t + '\n'


def has_overlap(s, t, k):
    # TODO this breaks when k > len(s) or len(t)
    suffix = s[-k:]
    prefix = t[:k]
    return suffix == prefix


def get_overlap_graph(fasta_dict, k):
    nodes = list(fasta_dict.keys())
    edges = []
    for label, dna in fasta_dict.items():
        for other_label, other_dna in [
            (other_label, other_dna) for other_label, other_dna in fasta_dict.items() if other_label != label
        ]:
            if has_overlap(dna, other_dna, k):
                edges.append((label, other_label))
    return nodes, edges



example_fasta_dict = {
    'Rosalind_0498': 'AAATAAA',
    'Rosalind_2391': 'AAATTTT',
    'Rosalind_2323': 'TTTTCCC',
    'Rosalind_0442': 'AAATCCC'
}

def calculating_expected_offspring(a, b, c, d, e, f):
    couple_list = [a, b, c, d, e, f]
    return couple_list[0]*2 + couple_list[1]*2 + couple_list[2]*2 + couple_list[3]*1.5 + couple_list[4]


if __name__ == '__main__':
    # gc_contents = gc_content('/Users/testtest/Downloads/rosalind_gc.txt')
    #label, percentage = find_key_of_max_value_in_dict(gc_contents)


    print(calculating_expected_offspring(17904, 17242, 17346, 17234, 19258, 17371))
