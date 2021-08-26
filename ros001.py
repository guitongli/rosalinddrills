#!/usr/bin/env python3
import sys


def no1(a, b):
    return int(a) ** 2 + int(b) ** 2


if __name__ == '__main__':
    print(no1(sys.argv[1], sys.argv[2]))



