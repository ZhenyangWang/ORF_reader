#! /usr/bin/env python

"""
Last updated on 2016-12-19
@author: ZhenyangWang
"""

from Bio import SeqIO

# The First Part
# define a function called look_for_orf
# translate by standard table in python
# set the minimum length for queried results
# load the fasta file
def look_for_orf(fasta_file):
    table = 1
    min_length = 400
    for myseq in SeqIO.parse(fasta_file, "fasta"):

# define the direction of reading sequence
# define the order of orf
# multiple of three
# translate the sequence to amino acid
# select the open reading frame larger than minimum length
# select the longest open reading frame
# print the result
        for strand, nuc in [(+1, myseq.seq), (-1, myseq.seq.reverse_complement())]:
            for orf in range(1, 4):
                length = 3 * ((len(myseq) - orf) // 3)
                for aminoacid in nuc[orf:orf + length].translate(table).split("*"):
                    if len(aminoacid) >= min_length:
                        min_length = len(aminoacid)
                        print("strand %i, orf %i, length %i, %s" \
                              % (strand, orf, len(aminoacid), aminoacid))
# print the result of individual file
look_for_orf("synseq.fasta")

# The Second Part
# define a function called decode
# load the fasta file
# assign string to s
def decode(fasta_file):
    for myseq in SeqIO.parse(fasta_file, "fasta"):
        s = myseq.seq

# convert A or C to 0 and G or T to 1 by using a while loop
        n = 0
        mi = ""
        while n < len(s):
            if (s[n] == 'A' or s[n] == 'C'):
                mi += '0'
            elif (s[n] == 'T' or s[n] == 'G'):
                mi += '1'
            n += 1
# intercept the characters from fourth to the end
        mi = mi[3:]
# using while loop to check every 8 bits binary number to the end
# using chr function convert the binary number to ASCII characters
        bincode = ""
        while mi != "":
            bincode += chr(int(mi[:8], 2))
            mi = mi[8:]
# converte the characters from 4000 to the 14000
# the characters from 0 to 4000 in this file could be converted to
# an unkonw terminator or clear command
        return bincode[4000:14000]

print decode("synseq.fasta")
