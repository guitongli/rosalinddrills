#!/usr/bin/env python3
import sys


def translate_dna_to_rna(dna: str):
    rna = ''
    for cha in dna:
        if cha.lower() == 't':
            rna += 'u'
        else:
            rna += cha.lower()

    return rna


if __name__ == '__main__':
    print(translate_dna_to_rna(sys.argv[1]).upper())



