#!/usr/bin/python

# remove SNV with smaller frequencies in the same gene

import sys

with open(sys.argv[1]) as f:
    nrAnatomicalSites = int(f.readline().split()[0])
    nrSamples = int(f.readline().split()[0])
    nrCharacters = int(f.readline().split()[0])
    header = f.readline().rstrip("\n")

    d = {}
    for line in f:
        line = line.rstrip("\n")
        s = line.rstrip("\n").split("\t")
        symbol = s[5].split(":")[-1]
        vaf = float(s[6])
        if symbol not in d:
            d[symbol] = []
        d[symbol].append((vaf, s))

print(nrAnatomicalSites, "#anatomical sites")
print(nrSamples, "#samples")
print(len(d), "#characters")
print(header,)
for idx, symbol in enumerate(d):
    s = sorted(d[symbol])[-1][1]
    s[4] = idx
    print("\t".join(map(str, s)))
