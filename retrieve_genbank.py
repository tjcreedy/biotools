#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Description of what the script does"""

# Imports

import os
import sys
import re
import time
import argparse


import numpy as np
import textwrap as _textwrap

from Bio import Entrez, SeqIO

# Class definitions

# This reclasses the argparse.HelpFormatter object to have newlines in the help text for paragrahps
class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width, initial_indent=indent,
                                                 subsequent_indent=indent
                                                 ) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text

# Global variables


# Function definitions

def str_is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def absent_authentication():
    out = ("# To automatically pull NCBI data from the internet, you need an NCBI account and API "
           "key.\n"
           "# Go to "
           "https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/\n"
           "# to find out more, register for an account and create an API key.\n"
           "#\n"
           "# Then create a text document with two lines, one for your registered email address,\n"
           "# one for your api key, as follows:\n"
           "my.email@address.com\n"
           "my9989api807key87qrd67ui7d9eu6di98eu\n"
           "#\n"
           "# Supply this file to the --ncbiauth or -n argument\n")
    with open("ncbi_authentication.txt", 'w') as o:
        o.write(out)
    os.chmod("ncbi_authentication.txt", 0o600)
    sys.stderr.write("\n\nTo pull NCBI data from the internet, you need to provide authentication."
                     "\nTo help, I have just made the file \"ncbi_authentication.txt\"\nPlease "
                     "read this for more information.\n\n")
    exit()


def get_authentication(path):
    # path = args.ncbiauth
    if not path:
        absent_authentication()
    auth = dict()
    emailregex = "^(\D)+(\w)*((\.(\w)+)?)+@(\D)+(\w)*((\.(\D)+(\w)*)+)?(\.)[a-z]{2,}$"
    with open(path, 'r') as fh:
        for line in fh:
            line = line.strip()
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
    return auth

def parse_searchstring(path):
    # Read file
    lines = path.readlines()
    path.close()
    
    # Catch issues
    if len(lines) == 0:
        sys.exit(f"Error: {path} is empty\n")
    elif len(lines) > 1:
        sys.stderr.write(f"Warning: {path} has more than one line, the search string should be on "
                         "one line; all the rest will be ignored")

    # Return string
    return(lines[0].strip())

def search(path, auth):
    
    # Parse the search
    searchstring = parse_searchstring(path)

    # Set up for Entrez
    Entrez.email = auth['email']
    Entrez.api_key = auth['key']
    Entrez.tool = "biotools/retrieve_genbank.py:search"

    # Carry out the search and extract results
    sys.stderr.write(f"Searching NCBI nt using the supplied search term\n")
    searchresult = Entrez.read(Entrez.esearch(db='nuccore', term=searchstring, retmax=999999999))
    gbaccessions = searchresult['IdList']

    # Fetch the GB IDs
    fetchresult = Entrez.efetch(db='nuccore', id=','.join(gbaccessions), rettype='acc', 
                                retmode='text')
    gbids = [i.strip() for i in fetchresult.readlines()]
    fetchresult.close()

    sys.stderr.write(f"Found {len(gbids)} records\n")
    return(gbids)

def parse_exclusion(paths):

    # Read file
    lines=[]
    for p in paths:
        lines.extend(p.readlines())
        p.close()

    # Remove anything after first character not expected to be in a GB id
    excl = [re.sub(r'[^A-Za-z0-9\._\|].*$', '', e.strip()) for e in lines]
    # Remove NAs
    excl = [e for e in excl if e != 'NA']
    # Separate any over multiple lines
    exclout = []
    for e in excl:
        exclout.extend(e.split('|'))

    return(exclout)

def exclude(ids, paths, exclnc = False):
    #ids, paths = fullids, args.exclude
    # Parse the exclusion list
    
    exclusion = parse_exclusion(paths)
    sys.stderr.write(f"Read {len(exclusion)} IDs to exclude\n")
    # Filter the ids
    outids = [g for g in ids if not any(e in g for e in exclusion)]
    sys.stderr.write(f"After excluding supplied IDs, {len(outids)} records remain\n")
    return(outids)

def retrieve_sequences(ids, auth, chunksize):
    #ids, auth, chunksize = fullids, args.ncbiauth, args.chunksize
    #idset = ids[:10]
    
    # Set up retrieval
    def efetch_sequences(idset):
        sh = Entrez.efetch(db='nucleotide', id=','.join(idset), rettype='gb', retmode='text')
        gb = SeqIO.parse(sh, 'genbank')
        gbgen = ( {'id': i, 'sr': g} for i, g in zip(idset, gb) )
        return gbgen
    sys.stderr.write(f"Retrieving {len(ids)} full records from NCBI nt\n")
    # Do retrieval
    ncbigen = retrieve_ncbi_remote(ids, efetch_sequences, 'id', chunksize, auth)

    # Parse into output object
    out = []
    for outsub in ncbigen:
        # outsub = next(ncbigen)
        for vals in outsub.values():
            out.append(vals['sr'])

    return(out)

def rename(seqrecords, rename):
    #seqrecords, rename = retseqs, args.rename
    
    # Set up renaming
    basename = rename
    startval = 1
    valwidth = len(str(len(retseqs)))
    extractint = re.search(r'[0-9]*$', rename)
    if extractint:
        basename = re.sub(r'[0-9]*$', '', rename)
        startval = extractint.group(0)
        extwidth = len(startval)
        if( valwidth > extwidth or len(str(int(startval) + len(retseqs))) > extwidth ):
            sys.exit("Error: supplied starting name to -r/--rename does not contain enough digits "
                     "to generate sufficient new names")
        valwidth = extwidth
        startval = int(startval)

    # Generate new names
    names = [f"{basename}{str(r).zfill(valwidth)}" for r in range(startval, startval+len(retseqs))]

    # Do renaming
    renamemap = {}
    for sr, n in zip(seqrecords, names):
        renamemap[n] = sr.name
        sr.name = n

    return seqrecords, renamemap

def parse_output_metadata(seqrecords, renamemap, path, chunksize, auth):
    #seqrecords, renamemap, path, chunksize, auth = retseqs, renamemap, args.metadata, args.chunksize, args.ncbiauth

    ranks = ["species","subgenus","genus","subtribe","tribe","subfamily","family","superfamily",
             "infraorder","suborder","order","class"]
    head1=["db_id","institution_code"]
    head2=["ncbi_species","ncbi_id_rank","ncbi_taxid", "locality","subregion","country","authors",
           "genbank_accession","contigname"]
    
    sys.stderr.write(f"Parsing metadata from {len(seqrecords)} retrieved records\n")
    metadata = {}
    for sr in seqrecords:

        fmtaut = lambda string : string.replace(', ', ' and ').replace(',', ' ')
        srdata = {
            'db_id'             : sr.name,
            'genbank_accession' : renamemap[sr.name],
            'authors'           : fmtaut(sr.annotations['references'][0].authors),
            'ncbi_id_rank'      : 'species'
        }
        
        srcfeat = [f for f in sr.features if f.type == 'source'][0]
        for var in ['db_xref', 'specimen_voucher', 'country']:
            #var='country'
            if var in srcfeat.qualifiers.keys():
                val = srcfeat.qualifiers[var]
                if var == 'db_xref':
                    val = [v for v in val if 'taxon' in v][0]
                    val = val.lstrip('taxon:')
                    srdata['ncbi_taxid'] = val
                else:
                    val = val[0]
                    if var == 'country':
                        val = re.split('[,:]+', val, maxsplit = 3)
                        for k, v in zip(['country', 'subregion', 'locality'], val):
                            srdata[k] = v.lstrip(' ').rstrip(' ')
                    else:
                        var = 'institution_code' if var == 'specimen_voucher' else var
                        srdata[var] =  val
        metadata[sr.id] = srdata
    
    # Get taxonomies
    taxids = set(srd['ncbi_taxid'] for srd in metadata.values())
    taxonomy = retrieve_taxonomy(taxids, chunksize, auth)
    lineages = {tid: get_standard_lineage(tax, ranks) for tid, tax in taxonomy.items()}

    # Finalise output
    for k, v in metadata.items():
        for r, t in zip(ranks, lineages[metadata[k]['ncbi_taxid']]):
            metadata[k][r] = t
        for ex, nw in zip(['species', 'genbank_accession'], ['ncbi_species', 'contigname']):
            metadata[k][nw] = metadata[k][ex]
    
    # Write out
    sys.stderr.write(f"Writing {len(metadata)} lines of metadata to {args.metadata.name}\n")
    header = head1 + ranks + head2
    path.write(','.join(header)+'\n')
    for mv in metadata.values():
        path.write(','.join([mv[h] if h in mv.keys() else '' for h in header]) + '\n')
    path.close()


def retrieve_ncbi_remote(ids, searchfunc, responsekey, chunksize, auth, maxerrors=10):
    #ids, searchfunc, responsekey = taxids, efetch_read_taxonomy, 'TaxId'
    #ids, searchfunc, responsekey, auth, maxerrors = absent, efetch_read_taxonomy, 'TaxId', authdict, 10
    Entrez.email = auth['email']
    Entrez.api_key = auth['key']
    Entrez.tool = "biotools/blast2taxonomy.py:retrieve_ncbi_remote"

    done = 0
    total = len(ids)
    idset = {ids.pop() for _ in range(min([chunksize, len(ids)]))}
    errors = 0
    maxiterations = int(total / chunksize * 5)
    it = 1
    missprop = []
    while 1:
        # Attempt to retrieve summaries for this set
        try:
            records = searchfunc(idset)
        # If any errors, retry this set, unless the max errors has been reached
        except BaseException as exception:
            errors += 1
            if errors < maxerrors:
                sys.stderr.write(f"\nHit a {type(exception)} in during retrieve_ncbi_remote - will "
                                 f"retry, debug dump follows:\nargs\n{exception.args}\nfull\n"
                                 f"{exception}\n")
                time.sleep(0.4)
                continue
            else:
                sys.stderr.write(f"Reached max errors of {maxerrors} during retrieve_ncbi_remote, "
                                 f"last error:\n {exception=}, {type(exception)=}")
                raise
        else:
            # If successful, work through the results received to store records
            out = {}
            for r in records:
                out[r[responsekey]] = r
            done += len(out)
            # If any missing, add them to the to-do set
            missingids = [t for t in idset if str(t) not in out]
            if len(missingids) > 0:
                missprop.append(len(missingids)/chunksize)
                ids.update(missingids)
            # Report, take a new set for the next iteration out of what is left
            sys.stderr.write(f"\rSuccessfully retrieved {done}/{total} records (iteration {it})")
            idset = {ids.pop() for _ in range(min([chunksize, len(ids)]))}
            if len(idset) > 0 and it < maxiterations:
                yield out
                it += 1
                time.sleep(0.4)
            else:
                if done < total:
                    sys.stderr.write(", failed to retrieve the remainder")
                meanmissprop = np.mean(missprop) if len(missprop) > 0 else 0
                sys.stderr.write(f", {errors} failed NCBI calls, "
                                 f"mean {round(meanmissprop, 3)*100}% complete NCBI returns.\n")
                yield out
                break


def retrieve_taxonomy(taxids, chunksize, authdict):
    #tidtaxdbpath, chunksize, authpath, authdict = args.tidtaxdb, args.chunksize, args.ncbiauth, auth
    """Search NCBI for lineage information given a tax id.
    """
    # taxids, chunksize, authdict = taxids, args.chunksize, args.ncbiauth
    def efetch_read_taxonomy(idset):
        sh = Entrez.efetch(db='taxonomy', id=[str(t) for t in idset])
        records = Entrez.read(sh)
        sh.close()
        return records

    sys.stderr.write(f"Searching NCBI taxonomy for {len(taxids)} taxonomies\n")
    ncbigen = retrieve_ncbi_remote(taxids, efetch_read_taxonomy, 'TaxId', chunksize, authdict)
    out = {}
    for outsub in ncbigen:
        # outsub = next(ncbigen)
        out.update(outsub)

    return out

def get_standard_lineage(record, ranks):
    # record = db[i]
    names = {r['Rank']: r['ScientificName'] for r in record['LineageEx']
             if r['Rank'] != 'no rank'}
    names[record['Rank']] = record['ScientificName']

    return [names[r] if r in names else '' for r in ranks]

def getcliargs(arglist=None):
    # Initialise the parser object and give a description
    parser = argparse.ArgumentParser(description="""
        This script searches GenBank using a given search term, then after optionally excluding a 
        given set or sets of GenBank IDs using -e/--exclude, it will download the full records for 
        the found sequences and output them in one or more genbank-format flat files. 
        |n
        One or more text files may be passed to the -e/--exclude option, if using. GenBank IDs 
        should be supplied at the beginning of each line, either one per line or separated by | if 
        on the same line. Anything after the first whitespace or comma will be ignored. If a line 
        contains only NA, it will be ignored.
        |n
        Any RefSeq sequences (with IDs starting NC_ may also be excluded automatically using the 
        -q/--norefseq argument
        |n
        You may optionally parse and output selected metadata from the genbank files using the 
        -m/--metadata option.
        |n
        To helpfully split the output sequences into those that have gene annotations and those 
        that do not, use the -u/--unannotated and -a/--annotated options.
        |n
        To rename the output sequences/metadata, use the -r/--rename option. If the argument 
        supplied ends in a character, e.g. "SEARCH1SEQ", this will be used as the prefix for the 
        names. If it ends in a number, this will be used as the ID for the first sequence, and the 
        number incremented for each further sequence
        |n
        This script retrieves information from NCBI, so you must authenticate with an email address 
        and API  key. Supply a file containing this information to -n/--ncbiauth. The script will 
        output an example file if NCBI authentication is required but absent.
        |n
        By default, the script will send 1000 ids to NCBI in each request. If this seems to cause 
        errors, set -z/--chunksize to a lower value. This cannot be increased.
        """, formatter_class=MultilineFormatter)

    # Add individual argument specifications
    parser.add_argument('-s', '--search', type=argparse.FileType('r'), metavar='PATH', 
                        required=True,
                        help='path to a one-line text file containing the NCBI search term to use')
    parser.add_argument('-r', '--rename', type=str, metavar='NAME',
                        help='if desired, rename the sequences based on the supplied string')
    parser.add_argument('-e', '--exclude', type=argparse.FileType('r'), metavar='PATH', nargs='*',
                        help='path to a text file containing GenBank IDs to exclude, one per line. '
                             'Any text from the first non-alphanumeric character will be excluded, '
                             'so the file may contain comments or other data')
    parser.add_argument('-q', '--norefseq', action='store_true',
                        help='exclude any RefSeq sequences (starting NC_) from the retrieved '
                             'sequences')
    parser.add_argument('-n', '--ncbiauth', type=str, metavar='PATH',
                        help='ncbi authentication path')
    parser.add_argument('-z', '--chunksize', type=int, metavar='N', default=1000,
                        help='number of ids per request, default 1000')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='PATH',
                        help='path to write all retrieved sequences to in a single file')
    parser.add_argument('-a', '--annotated', type=argparse.FileType('w'), metavar='PATH',
                        help='path to write only sequences with gene annotation to')
    parser.add_argument('-u', '--unannotated', type=argparse.FileType('w'), metavar='PATH',
                        help='path to write only sequences without gene annotation to')
    parser.add_argument('-m', '--metadata', type=argparse.FileType('w'), metavar='PATH',
                        help='parse selected metadata from the genbank record and write to a csv '
                             'at this path')

    # Parse the arguments from the function call if specified, otherwise from the command line
    args = parser.parse_args(arglist) if arglist else parser.parse_args()

    # Do some checking of the inputs
    if args.annotated or args.unannotated:
        if args.output:
            parser.error("Error: use -o/--output OR one or both of -a/--annotated and "
                         "-u/--unannotated")
        if not (args.annotated and args.unannotated):
            sys.stderr.write(f"Warning: only {'-u/--un' if args.unannotated else '-a/--'}annotated "
                             f"used, will discard any {'un' if args.annotated else ''}annotated "
                              "sequences retrieved")
    elif not args.output:
        sys.exit("Error: supply either -o/--output or one or both of -a/--annotated and "
                 "-u/--unannotated")
    if args.chunksize > 1000:
        parser.error("-c/--chunksize should not be greater than 1000")

    # Process the arguments
    args.ncbiauth = get_authentication(args.ncbiauth)

    # If the arguments are all OK, output them
    return args

# Start the actual script
if __name__ == "__main__":

    # Get the arguments
    args = getcliargs()
    sys.stderr.write("\n")
    # Post the search string and get back a list of IDs
    fullids = search(args.search, args.ncbiauth)

    # Exclude IDs if requested
    filtids = exclude(fullids, args.exclude) if args.exclude else fullids
    if args.norefseq:
        filtids = [i for i in filtids if "NC_" not in i]
        sys.stderr.write(f"After excluding RefSeq records, {len(filtids)} remain\n")

    # Retrive the sequence data
    retseqs = retrieve_sequences(filtids, args.ncbiauth, args.chunksize)

    # Rename if requested
    if args.rename:
        retseqs, renamemap = rename(retseqs, args.rename)
    else:
        renamemap = {g.id: g.id for g in retseqs}

    # Parse metadata if requested and output CSV
    if args.metadata:
        parse_output_metadata(retseqs, renamemap, args.metadata, args.chunksize, args.ncbiauth)

    # Separate into annotated and unannotated if requested and output GB(s)
    if args.output:
        sys.stderr.write(f"Writing {len(retseqs)} sequence records to {args.output.name}\n")
        SeqIO.write(retseqs, args.output, 'genbank')
        args.output.close()
    else:
        ann, unann = [], []
        for sr in retseqs:
            (ann, unann)[len([f for f in sr.features if f.type != 'source']) > 0].append(sr)
        if args.annotated:
            sys.stderr.write(f"Writing {len(ann)} annotated sequence records to "
                             f"{args.annotated.name}\n")
            SeqIO.write(ann, args.annotated, 'genbank')
            args.annotated.close()
        if args.unannotated:
            sys.stderr.write(f"Writing {len(unann)} unannotated sequence records to "
                             f"{args.unannotated.name}\n")
            SeqIO.write(unann, args.unannotated, 'genbank')
            args.unannotated.close()
