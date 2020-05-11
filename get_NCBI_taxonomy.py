#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Thomas J. Creedy
"""

# Imports

import sys
import os
import re
import argparse
import time
import pickle
import json
from Bio import Entrez


import textwrap as _textwrap

# Class definitions

class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width,
                                                 initial_indent=indent,
                                                 subsequent_indent=indent
                                                 ) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text

# Function definitions


def required_multiple(multiple):
    class RequiredMultiple(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not len(values) % multiple == 0:
                msg = 'argument "{f}" requires a multiple of {multiple} values'
                msg = msg.format(f=self.dest, multiple=multiple)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredMultiple

def parse_inputfile(path):
    data = dict()
    with open(path, 'r') as fh:
        ln = 0
        for line in fh:
            ln += 1
            items = line.strip().split('\t')
            if len(items) < 2:
                sys.stderr.write(f"Line number {ln} in {path} does not have",
                                 "at least 2 tab-separated elements")
                exit()
            data[items[0]] = items[1]
    return(data)

def retrieve_ncbi_local(taxids, database):
    #database = args.localdb
    out = dict()
    absent = set()
    db = dict()
    
    if os.path.exists(database):
        with open(database, 'r') as fp:
            db = json.load(fp)
        
        for i in taxids:
            if i in db:
                out[i] = db[i]
            else:
                absent.add(i)
    else:
        # Check that we can write to the file before spending time doing 
        # ncbi stuff
        with open(database, 'w') as fp:
            None
        os.remove(database)
        
        absent = taxids
    
    return(out, absent, db)

def absent_authentication():
    
    out = ("# To automatically pull NCBI data from the internet, you need an "
          + "NCBI account and API key.\n"
          + "# Go to https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/\n"
          + "# to find out more, register for an account and create an API "
          + "key.\n#\n"
          + "# Then create a text document with two lines, one for your "
          + "registered email address,\n# one for your api key, as follows:\n"
          + "my.email@address.com\n"
          + "my9989api807key87qrd67ui7d9eu6di98eu\n#\n"
          + "# Supply this file to the --ncbiauth or -n argument\n")
    with open("ncbi_authentication.txt", 'w') as o:
        o.write(out)
    os.chmod("ncbi_authentication.txt", 0o600)
    sys.stderr.write("\n\nTo pull NCBI data from the internet, you need to "
                     + "provide authentication.\nTo help, I have just made "
                     + "the file \"ncbi_authentication.txt\"\nPlease read "
                     + "this for more information.\n\n")
    exit()

def get_authentication(path):
    auth = dict()
    emailregex = "^(\D)+(\w)*((\.(\w)+)?)+@(\D)+(\w)*((\.(\D)+(\w)*)+)?(\.)[a-z]{2,}$"
    with open(path, 'r') as fh:
        for line in fh:
            line.strip()
            if line[0] == '#' or line == '': continue
            if re.match(emailregex, line):
                auth['email'] = line
            elif re.match("^[a-z\d]{36}$", line):
                auth['key'] = line
            else:
                sys.stderr.write(f"Error: don't recognise {line} in {path}")
                exit()
    if len(auth) != 2:
        sys.stderr.write(f"Error: could not find email and api key in {path}")
        exit()
    return(auth)

def chunker(seq, size):
    if type(seq) is set:
        seq = list(seq)
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

def retrieve_ncbi_remote(taxids, auth):
    """Search NCBI for lineage information given a tax id.
    """
    #taxids, auth = [absent, get_authentication(args.ncbiauth)]
    Entrez.email = auth['email']
    Entrez.api_key = auth['key']
    
    out = dict()
    
    for taxidset in chunker(taxids, 200):
        query = '&'.join(taxidset)
        handle = Entrez.efetch(db="taxonomy", id=query)
        records = Entrez.read(handle)
        handle.close()
        query.count('&')
        
        for r in records:
            out[r['TaxId']] = r
        
        time.sleep(0.1)
        
    
    absent = set(i for i in taxids if i not in out)
    
    return(out, absent)

def write_new_local(newdb, database):
    with open(database, 'w') as fp:
        json.dump(newdb, fp)

def get_standard_lineage(record, ranks):
    
    names = {r['Rank']: r['ScientificName'] for r in record['LineageEx'] 
             if r['Rank'] != 'no rank'}
    
    return([names[r] if r in names else '' for r in ranks])

def write_to_file(data, db, path, ranks):
    #data, db, path, ranks = [data, taxonomy, args.output, args.ranks]
    ranks = ranks.lower().split(',')
    with open(path, 'w') as oh:
        #oh = open(path, 'w')
        oh.write('\t'.join(['name', 'taxid'] + ranks) + '\n')
        for n, i in data.items():
            out = [n, i] 
            if i in db and 'LineageEx' in db[i]:
                out += get_standard_lineage(db[i], ranks)
            else:
                out += [''] * len(ranks)
            oh.write('\t'.join(out) + '\n')

# Argument parser

parser = argparse.ArgumentParser(
        description="""
description:
|n
This script retrieves taxonomy information for NCBI taxonomy id numbers 
assigned to a set of sequences, such as those output by MEGAN in the format 
\'readName_to_taxonId\'. The input should be a tab delimited table with the 
sequence name in the first column and NCBI taxid in the second column. All 
other columns will be ignored. The table should not have a header row.
|n
To efficiently retrieve taxonomy, the script uses a local database, supplied 
to --localdb. If you might already have one of these, ask your administrator 
for the path. If you do not, the script will create a new database at the 
location supplied.
|n
To retrieve information from NCBI, you must authenticate with an email address 
and api key. This is only necessary if your local database isn't complete. 
Supply a file containing this information to --ncbiauth. The script will output
 an example file if NCBI authentication is required but absent.
|n
The taxonomy information will be output in a table with the first two columns 
providing the input sequence ids and taxonomy ids. The subsequent columns will 
provide values for the standard seven taxonomic ranks. Cells will be blank if 
the given taxid is higher than the lowest taxonomic rank or if no taxonomy 
could be retrieved; in the latter case, a warning will also be printedt to 
terminal. A custom set of taxonomic ranks can be supplied to --ranks or -r in 
the format \'rank,rank,rank\'.""",
                                 formatter_class=MultilineFormatter)

parser._optionals.title = "arguments"

parser.add_argument('-i', '--input', help='input file path', type=str,
                    metavar = '', required = True)
parser.add_argument('-l', '--localdb', help='local database file path',
                    type=str, metavar = '', required = True)
parser.add_argument('-o', '--output', help='output file path', type=str,
                    metavar = '', required = True)
parser.add_argument('-n', '--ncbiauth', help='ncbi authentication path',
                    type=str, metavar = '')
parser.add_argument('-r', '--ranks', help='comma-separated list of ranks',
                    metavar='rank,rank,rank',
                    default='superkingdom,kingdom,phylum,class,order,family,'\
                             + 'genus,species')

# Main

if __name__ == "__main__":
    args = parser.parse_args()
    #arglist = "-i otus_nt_blast-MEGAN.txt -l pyNCBI.db -o otus_MEGAN_tax.tsv -n tjc_ncbi_authentication.txt".split(' ')
    #args = parser.parse_args(arglist)
    #os.chdir('/home/thomas/Documents/NHM_postdoc/Sequencing/NHMMar2019')
    
    # Parse the input file to dict
    data = parse_inputfile(args.input)
    
    # Get the unique taxids
    taxids = set(data.values())
    
    # Retrieve taxonomy from local if available
    taxonomy, absent, local = retrieve_ncbi_local(taxids, args.localdb)
    
    
    # If any absent, get from ncbi remote
    if len(absent) > 0:
        if not args.ncbiauth: absent_authentication()
        ncbi, ncbiabsent = retrieve_ncbi_remote(absent, 
                                             get_authentication(args.ncbiauth))
        
        if len(ncbi) > 0:
            local.update(ncbi)
            write_new_local(local, args.localdb)
            del local
            taxonomy.update(ncbi)
        
        if len(ncbiabsent) > 0:
            absent = ', '.join(list(ncbiabsent))
            text = "Warning: could not retrieve NCBI taxonomy data for the "\
                   + "following taxids."
            sys.stderr.write(f"{text}\n{absent}\n")
        
        
    
    write_to_file(data, taxonomy, args.output, args.ranks)
    exit()

