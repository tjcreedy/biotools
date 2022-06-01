#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Description of what the script does"""

# Imports

import os
import sys
import re
import json
import time
import argparse

import textwrap as _textwrap

from collections import Counter

from Bio import Entrez
from Bio.Blast import NCBIXML
# Test blast commands
# BLASTDB=/mbl/share/workspaces/groups/database/nt-2021-09-07/nt
# blastn -query ASVsub.fasta -db $BLASTDB -num_threads 10 -outfmt 5 -out test5.xml
# blastn -query ASVsub.fasta -db $BLASTDB -num_threads 10 -outfmt 6 -out test6.tsv
# blastn -query ASVsub.fasta -db $BLASTDB -num_threads 10 -outfmt 10 -out test10.csv
# blastn -query ASVsub.fasta -db $BLASTDB -num_threads 10 -outfmt "6 std staxids sscinames sskingdoms" -out test6t.tsv
# blastn -query ASVsub.fasta -db $BLASTDB -num_threads 10 -outfmt "10 std staxids sscinames sskingdoms" -out test10t.csv





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

def parse_title(title):
    #title = hit.title
    gbregex = r"((?:[A-Z]{2}_?)\d+(?:\.\d)?)"
    # See if this is a longform title separated by |
    tsplit = title.split('|')
    if len(tsplit) > 0:
        # Find the gb marker and take the part following
        gbloc = [i for i, part in enumerate(tsplit) if part in ['gb', 'ref', 'dbj', 'emb']]
        if len(gbloc) > 0:
            i = gbloc[0]
            gbsrch = re.search(gbregex, tsplit[i + 1])
            if gbsrch:
                return gbsrch.group(0)
        # If no recognisable gb marker, search all parts for an accession number
        else:
            gbcand = [part for part in tsplit if re.search(gbregex, part)]
            if len(gbcand) == 1:
                return re.search(gbregex, gbcand[0]).group(0)
    # Check to see if the title is a complete accession number by itself
    elif re.match(f"^{gbregex}$", title):
        return title
    # Otherwise, search the title for an accession number
    else:
        gbsrch = re.search(gbregex, title)
        if gbsrch:
            return gbsrch.group(0)
    sys.stderr.write(f"Error: cannot recognise a single genbank accession number in {title}")
    sys.exit()


def parse_tabular(path, sep, idthresh, taxc):
    out = {}
    with open(path) as fh:
        for line in fh:
            row = line.strip().split(sep)
            qseqid = row[0]
            # Get id and skip row if id below threshold or does not exceed the current top hit (if
            # using
            pident = float(row[2])
            if pident < idthresh:
                continue
            # Get hit details and record
            sseqid = parse_title(row[1])
            resd = {'id': sseqid, 'data': row}
            if taxc:
                resd['tx'] = int(row[taxc - 1].split(';')[0])
            # Add to output dict
            if qseqid in out:
                out[qseqid].append(resd)
            else:
                out[qseqid] = [resd]
    return out


def parse_xml(path, idthresh):
    out = {}
    results = NCBIXML.parse(open(path))
    for res in results:
        # res = next(results)
        qseqid = res.query
        hits = []
        for hit in res.alignments:
            # hit = res.alignments[8]
            sseqid = parse_title(hit.title)
            pident = (hit.hsps[0].identities/hit.hsps[0].align_length) * 100
            if pident >= idthresh:
                hits.append({'id': sseqid, 'data': [qseqid, sseqid, pident]})
        out[qseqid] = hits
    return out


def parse_input(path, tophit, idthresh, taxidc):
    #
    #path, idthresh, taxidc = "test6.tsv", 0, None
    # Find out the input format
    with open(path) as fh:
        fline = fh.readline().strip()
    if re.match(r"^<\?xml.*\?>$", fline):
        indata = parse_xml(path, idthresh)
    else:
        sep = ',' if fline.count(',') > fline.count('\t') else '\t'
        indata = parse_tabular(path, sep, idthresh, taxidc)

    # Filter out non-top hits if using
    if tophit:
        for qseqid in indata.keys():
            # qseqid = list(indata.keys())[0]
            pidents = [float(h['data'][2]) for h in indata[qseqid]]
            indata[qseqid] = [h for h in indata[qseqid] if float(h['data'][2]) == max(pidents)]

    # Get accessions to search and taxids if present
    accs = set()
    taxids = set()
    for qseqid, hits in indata.items():
        for hit in hits:
            if taxidc and 'tx' in hit:
                taxids.add(hit['tx'])
            else:
                accs.add(hit['id'])

    return indata, accs, taxids


def retrieve_ncbi_local(taxids, database):
    # database = args.localdb
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
            pass
        os.remove(database)

        absent = taxids

    return out, absent, db


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


def chunker(seq, size):
    if type(seq) is set:
        seq = list(seq)
    for pos in range(0, len(seq), size):
        yield seq[pos:pos + size]


def retrieve_ncbi_remote(taxids, chunksize, auth):
    """Search NCBI for lineage information given a tax id.
    """
    # taxids, chunksize, auth = absent, args.chunksize, get_authentication(args.ncbiauth)
    Entrez.email = auth['email']
    Entrez.api_key = auth['key']

    out = dict()

    for taxidset in chunker(taxids, chunksize):
        # taxidset = list(chunker(taxids, chunksize))[2]
        handle = Entrez.efetch(db="taxonomy", id=taxidset)
        records = Entrez.read(handle)
        handle.close()
        for r in records:
            out[r['TaxId']] = r
        time.sleep(0.1)

    absent = set(i for i in taxids if str(i) not in out)

    return out, absent


def retrieve_taxids(ids, chunksize, auth):
    # ids, db, auth, chunksize = gbaccs, 'nt', auth, 10
    Entrez.email = auth['email']
    Entrez.api_key = auth['key']

    out = dict()

    for idset in chunker(ids, chunksize):
        sh = Entrez.esummary(db='nucleotide', id=','.join(idset))
        summaries = Entrez.read(sh)
        sh.close()
        for gbid, smry in zip(idset, summaries):
            out[gbid] = int(smry['TaxId'])

    return out


def write_new_local(newdb, database):
    with open(database, 'w') as fp:
        json.dump(newdb, fp)


def get_standard_lineage(record, ranks):
    # record = db[i]
    names = {r['Rank']: r['ScientificName'] for r in record['LineageEx']
             if r['Rank'] != 'no rank'}
    names[record['Rank']] = record['ScientificName']

    return [names[r] if r in names else '' for r in ranks]


def assign_taxonomy(data, gbtax, taxonomies, ranks):
    # data, gbtax, taxonomies, ranks = inputdata, gbtaxids, taxonomy, args.ranks
    for qseqid in data.keys():
        # qseqid = list(data.keys())[0]
        for i in range(len(data[qseqid])):
            # i = 0
            # Retrieve the taxid
            if 'tx' in data[qseqid][i]:
                tx = data[qseqid][i]['tx']
            else:
                tx = data[qseqid][i]['tx'] = gbtax[data[qseqid][i]['id']]
            # Retrieve the taxonomy for this taxid and format it
            if str(tx) in taxonomies:
                data[qseqid][i]['taxonomy'] = get_standard_lineage(taxonomies[str(tx)], ranks)
            else:
                data[qseqid][i]['taxonomy'] = ['']  * len(ranks)
    return data


def lca(data, ranks):
    # data = taxonomised
    out = {}
    for qseqid, hits in data.items():
        # qseqid = list(data.keys())[0]
        lcasets = {r: set() for r in ranks}
        for hit in hits:
            # i = 0
            for j, r in enumerate(ranks):
                # j = 0
                lcasets[r].add(hit['taxonomy'][j])
        lcataxonomy = []
        for r in ranks:
            if len(lcasets[r]) == 1:
                lcataxonomy.append(lcasets[r].pop())
        out[qseqid] = [{'taxonomy': lcataxonomy}]
    return out


def filter_tophit(data):
    #data = taxonomised
    out = {}
    for qseqid, hits in data.items():
        # qseqid, hits = list(data.items())[0]
        if len(hits) > 1:
            # Count the unique taxonomies
            taxcount = Counter()
            for hit in hits:
                taxcount[','.join(hit['taxonomy'])] += 1
            # Find the most frequent taxonomy
            toptaxonomy = [t for t, c in taxcount.items() if c == max(taxcount.values())][0]
            # Output the first hit with the most frequent taxonomy
            out[qseqid] = [[h for h in hits if ','.join(h['taxonomy']) == toptaxonomy][0]]
        else:
            out[qseqid] = hits
    return out


def writeout(data):
    for qseqid, hits in data.items():
        # qseqid, hits = list(data.items())[0]
        for hit in hits:
            output = []
            if 'data' in hit:
                output.extend(hit['data'])
            else:
                output.append(qseqid)
            if 'id' in hit:
                output.append(hit['id'])
            if 'tx' in hit:
                output.append(str(hit['tx']))
            output.extend(hit['taxonomy'])
            sys.stdout.write('\t'.join(output) + '\n')


def getcliargs(arglist=None):
    # Initialise the parser object and give a description
    parser = argparse.ArgumentParser(description="""
        Parse the output of a BLAST search against an NCBI database (e.g. nt), supplied to 
        -b/--blastresults, to assign taxonomy to queries. BLAST results can be in XML 
        (--outfmt 5), tsv (--outfmt 6) or csv (--outfmt 10). If tsv or csv format, the first three 
        columns must be query seq id, subject seq id and percent identity in that order, as in the 
        default tabular output format. If a non-default format has been run that includes the 
        column staxids, the column number for this data should be passed to -t/--taxidcolumn. Note 
        that if multiple taxids are found for a hit, the first will be used. 
        |n
        The taxonomy will be retrieved from NCBI Taxonomy based on the taxid of the hit(s). By 
        default, the taxonomy of all hits will be retrieved. Optionally, the hits can be filtered, 
        either keeping only the/a top hit using the -p/--tophit option, or finding the lowest 
        common ancestor of all hits using the -a/--lca option. 
        |n
        The script will output a tsv to STDOUT. With the default option, the script will output 
        the information each hit for each query id, in the same format as supplied, followed by 
        the parsed NCBI accession ID, the NCBI taxid and taxonomy columns for each hit. When using 
        the -p/--tophit option, only a top hit for each query will be output. If multiple top 
        hits are found, the first one with the most frequent taxonomy will be output. When using 
        the -a/--lca option, the script will output a single row for each query id and the LCA  
        taxonomy columns. Note that if the input is XML, the only fields parsed from the XML and 
        output are qseqid, sseqid and pident. Control the taxonomic ranks output in all cases 
        using -r/--rank. 
        |n
        BLAST results can optionally be filtered to remove hits below a percent identity threshold 
        prior to retrieving taxonomy. Supply an id threshold to -i/--idthreshold, default 0. It is 
        advisable to use this option when filtering with the -a/--lca option.
        |n
        To efficiently retrieve taxonomy, the script uses a local database, supplied to 
        -l/--localdb, that may have been created by previous runs of this script. If this is your 
        first run, the script will create a new database at the location supplied.
        |n
        To retrieve information from NCBI, you must authenticate with an email address and API 
        key. This is only necessary if your local database isn't complete. Supply a file 
        containing this information to -n/--ncbiauth. The script will output an example file if 
        NCBI authentication is required but absent.
        |n
        By default, the script will send 200 id numbers to NCBI in each request. If this seems to 
        cause errors, set --chunksize to a lower value. 
        """, formatter_class=MultilineFormatter)

    # Add individual argument specifications
    parser.add_argument('-b', '--blastresults', required=True, metavar='PATH',
                        help='path to blast results file')
    parser.add_argument('-p', '--tophit', action='store_true',
                        help='return the taxonomy of the top hit(s) by percent id')
    parser.add_argument('-a', '--lca', action='store_true',
                        help='return the lowest common ancestor of all hits for each query')
    parser.add_argument('-i', '--idthreshold', type=float, metavar='N', default=0,
                        help='minimum percent id to retrieve taxonomy for a hit, default 0')
    parser.add_argument('-l', '--localdb', type=str, metavar='PATH', required=True,
                        help='local database file path')
    parser.add_argument('-n', '--ncbiauth', type=str, metavar='PATH',
                        help='ncbi authentication path')
    parser.add_argument('-t', '--taxidcolumn', type=int, metavar='N',
                        help='the column number of the staxid field if included in an input tsv '
                             'or csv')
    parser.add_argument('-r', '--ranks', metavar='rank,rank,rank',
                        help='comma-separated list of ranks',
                        default='superkingdom,kingdom,phylum,class,order,family,genus,species')
    parser.add_argument('-c', '--chunksize', type=int, metavar='N', default=500,
                        help='number of ids per request, default 500')

    # Parse the arguments from the function call if specified, otherwise from the command line
    args = parser.parse_args(arglist) if arglist else parser.parse_args()

    # Do some checking of the inputs
    if args.tophit and args.lca:
        parser.error("select only one of -p/--tophit or -a/--lca")

    # Process arguments
    args.ranks = args.ranks.lower().split(',')

    # If the arguments are all OK, output them
    return args

# Start the actual script
if __name__ == "__main__":

    # Get the arguments
    args = getcliargs()
    auth = None

    # Parse the inputs, filtering out results below idthreshold
    inputdata, gbaccs, taxids = parse_input(args.blastresults, args.tophit, args.idthreshold,
                                            args.taxidcolumn)

    # If taxids not present, search the unique accession numbers to retreive taxids
    gbtaxids = dict()
    if len(gbaccs) > 0:
        auth = get_authentication(args.ncbiauth)
        gbtaxids = retrieve_taxids(gbaccs, args.chunksize, auth)
        # Add to the master list of taxids
        taxids.update(set(gbtaxids.values()))

    # Retrieve taxonomy from local if available
    taxonomy, absent, local = retrieve_ncbi_local(taxids, args.localdb)

    # If any absent, get from ncbi remote
    if len(absent) > 0:
        if not auth:
            auth = get_authentication(args.ncbiauth)
        ncbi, ncbiabsent = retrieve_ncbi_remote(absent, args.chunksize, auth)

        if len(ncbi) > 0:
            local.update(ncbi)
            write_new_local(local, args.localdb)
            del local
            taxonomy.update(ncbi)

        if len(ncbiabsent) > 0:
            absent = ', '.join(list(ncbiabsent))
            text = "Warning: could not retrieve NCBI taxonomy data for the " \
                   + "following taxids."
            sys.stderr.write(f"{text}\n{absent}\n")

    # Assign taxonomy to the input data
    taxonomised = assign_taxonomy(inputdata, gbtaxids, taxonomy, args.ranks)

    # Do LCA analysis
    if args.lca:
        taxonomised = lca(taxonomised, args.ranks)

    # Filter top hits
    if args.tophit:
        taxonomised = filter_tophit(taxonomised)
    # Output
    writeout(taxonomised)

