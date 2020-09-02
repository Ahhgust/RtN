#!/usr/bin/env python3

import random
import fileinput
import sys


def readFasta(filename):
  seq = []
  with open(filename) as fh:
    for line in fh:
      if line[0] == ">":
        seq = []
      else:
        seq.append(line.rstrip())
      

  return "".join(seq)

def decode(tok):
  '''
  Takes in a token from empop (73R, 151.1C)
  and returns a triple of position
  allele (nucleic acid)
  and type 
  X=mismatch, I=insertion, D=deletion
  THis is likely not exhaustive
  '''
  idx = tok.find(".")
  # insertion
  if idx>-1:
    allele=tok[idx+1:].upper()
    pos= int(tok[ :idx])
    t = "I"
    q=0
    i=0
    # convert .1C into C
    # and .2C into CC
    while allele[i] >= '0' and allele[i] <= '9':
      q = 10*q + ord(allele[i]) - ord('0')
      i += 1

    if i > 0:
      allele = allele[i:] * q

  elif tok.upper().endswith("DEL"): # deletion case
    t = "D"
    allele = ""
    pos = int(tok[:-3])
  else:
    allele = tok[-1].upper()
    t="X" #mismatch
    pos = int(tok[:-1])

  return(pos, allele, t)

ref = readFasta(sys.argv[1])

iupac = {
  "R" : "AG",
  "Y" : "TC",
  "K" : "GT",
  "M" : "AC",
  "S" : "GC",
  "W" : "AT"}


for line in fileinput.input(sys.argv[2:]):
    if line[0] == '#':
      continue

    scaffold = list(ref)
    haps = line.rstrip().split("\t")
    who = haps[0]
    haps = haps[3:]
    
    for h in haps:
      (p,a,t) = decode(h)
      p -= 1 # 0-based coordinate
      if t == 'I':
        scaffold[p] += a
      elif t == 'D':
        scaffold[p] = ""
      elif a in iupac:
        if random.random() < 0.5:
          scaffold[p] = iupac[a][0]
        else:
          scaffold[p] = iupac[a][1]

      elif a in "GACT":
        scaffold[p] = a
      else:
        print("Oops! " , h, p, a, t, file=sys.stderr)
        sys.exit(1)

    print(">" , who,sep="")
    if scaffold[3106]=='N':
      print("".join(scaffold[0:3106]) , "".join(scaffold[3107:]) , sep="")
    else:
      print("".join(scaffold))
    

