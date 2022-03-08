#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Thomas J. Creedy
"""
# Imports
import sys
import urllib.request
# Functions
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))
# Main
acclist = sys.stdin.read().splitlines()
chunks = 100
n = 0
for chunk in chunker(acclist, chunks):
    start = chunks * n + 1
    sys.stderr.write("getting genbanks chunk %s-%s\n" % 
                     (start, start + len(chunk) - 1))
    accstring = ','.join(chunk)
    url = ("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
          + "db=nuccore&id=" + accstring + "&rettype=gb")
    sys.stdout.write(urllib.request.urlopen(url).read().decode('utf-8'))
    n += 1
exit()