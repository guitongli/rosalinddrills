#!/usr/bin/env python3

import sys

#5
def no5():
    string=open(sys.argv[1],'r').readlines()[0].strip()
    wordlist=string.split (' ')
    mydict= {i:wordlist.count(i) for i in wordlist}
    
    for key, value in mydict.items():
        print(key + ' ' + str(value))
        
if __name__ == "__main__":
    no5()



"""
#2

s='rNqQtWs3qtj0ScPTEF0wZfGPkAiz4s2NjVJFdBK68gGUOUE3Candoia99p60sJxV43I61TJlCPVYp9rF5yfUo90lQfpIdi9rTejclv53gpVL818yzyLk5fnbyJleTEn6HcPLxnci7jgDxLvwIfKbWcqekPkP3LdnBk8sxgoZussuriensisiOlAFB6r6PHHpBy.'
a=22
b=27
c=97
d=102
print([s[a:b+1],s[c:d+1]])


#3
sum=0
for num in range(4481, 9145):
    if num%2==1:
        sum=sum+num

print(sum)

#4
def no4():
    f=open('rosalind_ini5.txt','r').readlines()
    i=1
    while i< len(f):
        print (f[i])
        i=i+2
"""
