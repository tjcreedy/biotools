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

import numpy as np
import textwrap as _textwrap

from Bio import Entrez
from Bio.Blast import NCBIXML

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

class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end

class MinimumInteger(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values < 0:
            parser.error("Minimum value for {0} is 1".format(option_string))
        setattr(namespace, self.dest, values)


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
    gbregex = r"(?:^|\s)((?:[A-Z]{1,2}_?)\d{3,}(?:\.\d+)?)(?:$|\s)"
    # If this is a longform title with multiple entries, remove all but the first
    title = title.split(' >', 1)[0]
    # See if this is a longform title separated by |
    tsplit = title.split('|')
    if len(tsplit) > 0:
        # Find the gb marker and take the part following
        marker = ['gb', 'ref', 'dbj', 'emb', 'tpg', 'tpe']
        gbloc = [i for i, part in enumerate(tsplit) if part in marker]
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
    elif re.match(gbregex, title):
        return title
    # Otherwise, search the title for an accession number
    else:
        gbsrch = re.search(gbregex, title)
        if gbsrch:
            return gbsrch.group(0)
    sys.stderr.write(f"Warning: cannot recognise a single genbank accession number in {title} - "
                     f"this hit will be ignored")
    #sys.exit()


def parse_tabular(path, sep, taxc, minid, minscore, minalen):
    out = {}
    with open(path) as fh:
        for line in fh:
            row = line.strip().split(sep)
            qseqid = row[0]
            # Get id and skip row if id below threshold or does not exceed the current top hit (if
            # using
            pident = float(row[2])
            bitscore = float(row[11])
            length = int(row[3])
            if pident <  minid or bitscore < minscore or length < minalen:
                continue
            # Get hit details and record
            sseqid = parse_title(row[1])
            resd = {'id': sseqid,
                    'data': row,
                    'pident': pident,
                    'evalue': float(row[10]),
                    'bitscore': bitscore,
                    'length': length}
            if taxc:
                resd['tx'] = int(row[taxc - 1].split(';')[0])
            # Add to output dict
            if qseqid in out:
                out[qseqid].append(resd)
            else:
                out[qseqid] = [resd]
    return out


def parse_xml(path, minid, minscore, minalen):
    out = {}
    results = NCBIXML.parse(open(path))
    for res in results:
        # res = next(results)
        qseqid = res.query
        hits = []
        for hit in res.alignments:
            # hit = res.alignments[1]
            sseqid = parse_title(hit.title)
            hsp = hit.hsps[0]
            pident = (hsp.identities/hsp.align_length) * 100
            evalue = hsp.expect
            bitscore = hsp.bits
            length = hsp.score
            if pident >= minid and bitscore >= minscore and length >= minalen :
                hits.append({'id': sseqid,
                             'data': [qseqid, sseqid, pident, evalue, bitscore],
                             'pident': pident,
                             'evalue': evalue,
                             'bitscore': bitscore,
                             'length': length})
        out[qseqid] = hits
    return out


def parse_input(path, taxidc, minscore = 0, minid = 0, minalen = 0):
    # path, taxidc, minscore, minid, minalen = args.blastresults, args.taxidcolumn, args.minscore, args.minid, args.minalen
    # Find out the input format
    with open(path) as fh:
        fline = fh.readline().strip()
    if re.match(r"^<\?xml.*\?>$", fline):
        indata = parse_xml(path, minid, minscore, minalen)
    else:
        sep = ',' if fline.count(',') > fline.count('\t') else '\t'
        indata = parse_tabular(path, sep, taxidc, minid, minscore, minalen)

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


def retrieve_ncbi_remote(ids, searchfunc, responsekey, chunksize, auth, maxerrors=10):
    #ids, searchfunc, responsekey, chunksize, auth = absent, esummary_read_taxids, 'AccessionVersion', chunksize, authdict
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
                sys.stderr.write(f", {errors} failed NCBI calls"
#                                 f", mean {round(meanmissprop, 3)*100}% complete NCBI returns"
                                 f".\n")
                yield out
                break


def retrieve_ncbi_local(ids, dbpath):
    # dbpath = gbtiddbpath
    # ids, dbpath = taxids, tidtaxdbpath
    out = dict()
    absent = set()

    if os.path.exists(dbpath):
        with open(dbpath, 'r') as fp:
            db = json.load(fp)

        for i in ids:
            if str(i) in db:
                out[str(i)] = db[str(i)]
            else:
                absent.add(i)
    else:
        # Check that we can write to the file for future use
        with open(dbpath, 'w') as fp:
            pass
        os.remove(dbpath)
        absent = {i for i in ids}

    return out, absent


def update_ncbi_local(newdb, dbpath):

    olddb = dict()

    if os.path.exists(dbpath):
        with open(dbpath, 'r') as fp:
            olddb = json.load(fp)

    for k, v in newdb.items():
        olddb[k] = v

    with open(dbpath, 'w') as fp:
        json.dump(olddb, fp)


def retrieve_taxids(gbids, gbtiddbpath, chunksize, authpath=None, authdict=None, searchreps = 10):
    # gbids, gbtiddbpath, chunksize, authpath, auth = gbaccs, args.gbtiddb, args.chunksize, args.ncbiauth, None

    out, absent = retrieve_ncbi_local(gbids, gbtiddbpath)
    if len(out) > 0:
        sys.stderr.write(f"Retrieved {len(out)} taxids from local database\n")

    if len(absent) > 0:
        if not authdict:
            if not authpath:
                sys.exit("Supply one of authpath or auth to retrieve_taxids")
            authdict = get_authentication(authpath)

        def esummary_read_taxids(idset):
            sh = Entrez.esummary(db='nucleotide', id=','.join(idset))
            summaries = Entrez.read(sh, validate = False)
            sh.close()
            return summaries
        sys.stderr.write(f"Searching NCBI nt for {len(absent)} taxids\n")
        ncbigen = retrieve_ncbi_remote(absent, esummary_read_taxids, 'AccessionVersion', chunksize,
                                       authdict, maxerrors = searchreps)
        rem = {}
        for outsub in ncbigen:
            # outsub = next(ncbigen)
            for gb, smry in outsub.items():
                if type(smry['TaxId']) is Entrez.Parser.IntegerElement:
                    rem[gb] = int(smry['TaxId']) 

        update_ncbi_local(rem, gbtiddbpath)

        out.update(rem)

    absent = set(i for i in gbids if str(i) not in out)

    return out, absent, authdict


def retrieve_taxonomy(taxids, tidtaxdbpath, chunksize, authpath=None, authdict=None, 
                      searchreps = 10):
    #tidtaxdbpath, chunksize, authpath, authdict = args.tidtaxdb, args.chunksize, args.ncbiauth, auth
    """Search NCBI for lineage information given a tax id.
    """
    # taxids, chunksize = absent, 1000

    out, absent = retrieve_ncbi_local(taxids, tidtaxdbpath)
    if len(out) > 0:
        sys.stderr.write(f"Retrieved {len(out)} taxonomies from local database\n")

    if len(absent) > 0:
        if not authdict:
            if not authpath:
                sys.exit("Supply one of authpath or auth to retrieve_taxonomy")
            authdict = get_authentication(authpath)

        def efetch_read_taxonomy(idset):
            sh = Entrez.efetch(db='taxonomy', id=[str(t) for t in idset])
            records = Entrez.read(sh, validate = False)
            sh.close()
            return records

        sys.stderr.write(f"Searching NCBI taxonomy for {len(absent)} taxonomies\n")
        ncbigen = retrieve_ncbi_remote(absent, efetch_read_taxonomy, 'TaxId', chunksize, authdict, 
                                       maxerrors = searchreps)
        rem = {}
        for outsub in ncbigen:
            # outsub = next(ncbigen)
            rem.update(outsub)

        update_ncbi_local(rem, tidtaxdbpath)

        out.update(rem)

    absent = set(i for i in taxids if str(i) not in out)

    return out, absent


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
            # i = 462
            # Retrieve the taxid
            if 'tx' in data[qseqid][i]:
                tx = data[qseqid][i]['tx']
            else:
                gbid = data[qseqid][i]['id']
                if gbid in gbtax:
                    tx = data[qseqid][i]['tx'] = gbtax[gbid]
                else:
                    tx = data[qseqid][i]['tx'] = None
            # Retrieve the taxonomy for this taxid and format it
            if str(tx) in taxonomies:
                data[qseqid][i]['taxonomy'] = get_standard_lineage(taxonomies[str(tx)], ranks)
            else:
                data[qseqid][i]['taxonomy'] = ['']  * len(ranks)
    return data

def lca(data, ranks):
    # data, ranks = taxonomised, args.ranks
    out = {}
    for qseqid, hits in data.items():
        # qseqid, hits = list(data.items())[422]
        lcasets = {r: set() for r in ranks}
        for hit in hits:
            # i = 0
            for j, r in enumerate(ranks):
                # j = 0
                lcasets[r].add(hit['taxonomy'][j])
        lcataxonomy = []
        for r in ranks:
            lcataxonomy.append(lcasets[r].pop() if len(lcasets[r]) == 1 else '')
        out[qseqid] = {'taxonomy': lcataxonomy}
    return out

def megan_naive_lca(data, ranks, minscore, maxexp, minid, toppc, winid, minhitpc, minhitn, minlen, 
                    minsupppc):
    # data, ranks, minscore, maxexp, minid, toppc, winid, minhitpc, minhitn, minlen, minsupppc = taxonomised, args.ranks, args.minscore, args.maxexp, args.minid, args.winid, args.toppc, args.minhitpc, args.minhitn, args.minalen, args.minsupppc
    
    out = {}
    for qseqid, hits in data.items():
        # qseqid, hits = list(data.items())[2]
        # print(json.dumps(hits[0:20], indent = 4))
        # Filter out any hit rejected for the various filters
        scores = [h['bitscore'] for h in hits]
        mintopscore = (100 - toppc)/100 * max(scores)
        winidactive = winid and max([h['pident'] for h in hits]) >= winid
        for hit in hits:
            if hit['bitscore'] < minscore:
                hit['status'] = 'reject-minscore'
            elif hit['bitscore'] < mintopscore:
                hit['status'] = 'reject-outsidetoppc'
            elif hit['pident'] < minid:
                hit['status'] = 'reject-minid'
            elif hit['evalue'] > maxexp:
                hit['status'] = 'reject-maxexp'
            elif hit['length'] < minlen:
                hit['status'] = 'reject-minlen'
            elif winidactive and hit['pident'] < winid:
                hit['status'] = 'reject-belowwinid'
            else:
                hit['status'] = 'acceptforlca'
        
        # Do LCA
        lcaids, lcascores = [], []
        lcasets = {r: dict() for r in ranks}
        for hit in hits:
            # hit = hits[0]
            if hit['status'] != 'acceptforlca':
                continue
            
            lcaids.append(hit['pident'])
            lcascores.append(hit['bitscore'])
            for j, r in enumerate(ranks):
                # j = 1
                taxlist = hit['taxonomy'][:(j+1)]
                taxhash = str(hash(tuple(taxlist)))
                if taxhash in lcasets[r].keys():
                    lcasets[r][taxhash]['n'] += 1
                else:
                    lcasets[r][taxhash] = {'taxonomy': taxlist, 'n': 1}
        
        # Do LCA
        suppn = 0
        lcataxonomy = []
        if len(lcaids) > 0:
            minsuppn = minsupppc/100 * len(lcaids)
            for r in ranks[::-1]:
                # r = 'superkingdom'
                # Filter out taxonomies without sufficient support
                maxn = max(d['n'] for t, d in lcasets[r].items())
                if maxn < minsuppn:
                    continue
                sets = {t:d for t, d in lcasets[r].items() if d['n'] == maxn}
                if len(sets) > 1:
                    continue
                else:
                    finaltaxonomy = list(sets.values())[0]
                    lcataxonomy = finaltaxonomy['taxonomy']
                    suppn = finaltaxonomy['n']
                    break
        
        # Pad LCA taxonomy
        lcataxonomy += [''] * (len(ranks) - len(lcataxonomy))

        # Start output
        out[qseqid] = {'hits': len(hits),
                    'considered': len(lcaids),
                    'consideredpc': 100 * len(lcaids)/len(hits),
                    'supportingn': suppn}

        # Check sufficient support and report taxonomy and support
        if len(lcaids) < minhitn or len(lcaids)/len(hits) < minhitpc/100:
            out[qseqid].update({'taxonomy': ['' for r in ranks],
                                'support': 'reject-insufficient'})
        else:
            out[qseqid].update({'taxonomy': lcataxonomy,
                                'support': 'accept-sufficient'})

        # Report stats for considered hits
        if len(lcaids) > 0:
            out[qseqid].update(
                {'consideredminpident': min(lcaids),
                 'consideredmaxpident': max(lcaids),
                 'consideredminscore': min(lcascores),
                 'consideredmaxscore': max(lcascores)}
                 )
        else:
            out[qseqid].update(
                {k: '' for k in ['consideredminpident', 'consideredmaxpident', 
                                   'consideredminscore', 'consideredmaxscore']}
                )
    
        # Add info for top hit
        out[qseqid].update(
            {'toppident': hits[0]['pident'],
             'topscore': hits[0]['bitscore'],
             'toplen': hits[0]['length'],
             'topstatus': hits[0]['status']}
             )

    return out, data

def filter_tophit(data, method):
    # data, method = taxonomised, "bitscore"
    altmethod = "pident" if method == "bitscore" else "bitscore"
    #data = taxonomised
    out = {}
    for qseqid, hits in data.items():
        # qseqid, hits = list(data.items())[0]
        max1 = max(h[method] for h in hits)
        max2 = max(h[altmethod] for h in hits if h[method] == max1)
        selected = [h for h in hits if h[method] == max1 and h[altmethod] == max2]
        out[qseqid] = {'taxonomy': selected[0]['taxonomy'], 'ntop': len(selected)}
    return out

def writetaxonomy(data, path, ranks):
    #data, path, ranks = taxonomies, args.outtaxonomy, args.ranks
    header = list(data.values())[0].keys()
    header = [h for h in header if h != "taxonomy"]
    fh = open(path, 'w') 
    fh.write(','.join(['qseqid'] + ranks + list(header)) + '\n')
    for q, v in data.items():
        # q, v = list(data.items())[0]
        fh.write(','.join([q] + v['taxonomy'] + [str(v[h]) for h in header]) + '\n')
    fh.close()

def writehits(data, path, ranks):
    #data, path, ranks = taxonomised, args.outhits, args.ranks
    header = list(data.values())[0][0].keys()
    header = [h for h in header if not h in ["data", "taxonomy"]]
    fh = open(path, 'w')
    fh.write(','.join(['qseqid'] + list(header) + ranks) + '\n')
    for q, hits in data.items():
        # q, hits = list(data.items())[0]
        for hit in hits:
            # hit = hits[462]
            fh.write(','.join([q] + [str(hit[h]) for h in header] + hit['taxonomy']) + '\n')
    fh.close()

def getcliargs(arglist=None):
    # Initialise the parser object and give a description
    parser = argparse.ArgumentParser(description="""
        Parse the output of a BLAST search against an NCBI database (e.g. nt), supplied to 
        -b/--blastresults, to assign taxonomy to queries. BLAST results can be in XML 
        (--outfmt 5), tsv (--outfmt 6) or csv (--outfmt 10). If tsv or csv format, the query seq 
        id, subject seq id, percent identity, evalue and bitscore must be in columns 1, 2, 3 11 and 
        12 respectively, as in the default tabular output format. If a non-default format has been 
        run that includes the column "staxids", the column number for this data should be passed to 
        -c/--taxidcolumn. Note that if multiple taxids are found for a hit, the first will be used. 
        |n
        To skip searching for taxonomy for poor quality hits, these can optionally be filtered out 
        first according to percent identity (-i/--minid), bitscore (-s/--minscore) and/or alignment 
        length (-l/--minalen). By default, no filtering is done.
        |n
        The taxonomy will be retrieved from NCBI Taxonomy based on the taxid of the (remaining) 
        hit(s). To efficiently retrieve taxonomy, the script uses a pair of local datasets in json 
        files, one recording taxids for GenBank accession numbers, one recording taxonomy for 
        taxids. These might have been created by a previous run of this script, or the latter with 
        get_NCBI_taxonomy.py. Supply a path to the accession-taxid json with -g/--gbtiddb and the 
        taxid-taxonomy json with -x/--tidtaxdb. If this is your first run, the script will create 
        new databases at the locations supplied.
        |n
        Optionally, hits can then be processed to select a single taxonomy for each query using the 
        option -p/--process. Currently, three options are available. The two simple options are Top 
        Hit and LCA. Top Hit returns the taxonomy of the hit with the highest bitscore (-p 
        "top_score") or highest percent identity (-p "top_id"). Where there are ties, 
        these are broken by percent identity and bitscore respectively; if ties remain, the first 
        hit used. LCA (-p lca) simply returns the lowest common ancestor of all hits for which 
        taxonomy has been retrieved. 
        |n
        The final processing option, MEGAN LCA (-p megan), follows closely to the MEGAN LCA 
        algorithm as originally described here: https://genome.cshlp.org/content/17/3/377.full. 
        When this process is selected, no hit filtering prior to taxonomy retrieval is performed, 
        so taxonomy is retrieved for all hits. Then, more detailed filtering is performed to reject 
        hits prior to LCA analysis. Hits are reject if they fail the following conditions, in order:
        1. The bitscore is less than the minimum supplied to -s/--minscore (default 1); 
        2. The bitscore is less than a given percentage (-t/--toppc) below the top scoring hit 
        (default 15);
        3. The percent identity is less than the minimum supplied to -i/--minid (default 80);
        4. The expected value (e-value) is greater than the maximum supplied to -e/--maxexp 
        (default 0.01);
        5. The alignment length is less than the minimum supplied to -l/--minalen (default 100).
        Then, if -w/--winid is set, only hits meeting or exceeding this percent identity are 
        retained, and the rest discarded; by default this is not set but this can be set to 99 or 
        100 to achieve better species-level assignments.
        The final number of remaining hits is then compared against the supplied minimum number 
        (-n/--minhitn, default 0) and percentage of all hits (-a/--minhitpc, default 1). If 
        sufficient hits remain, the lowest common ancestor of these hits is found. The lowest 
        common ancestor procedure is parameterised by -f/--minsupppc; the LCA procedure finds the 
        taxonomy shared by at least this percentage of hits. Set this to 100 to find the exact 
        lowest common ancestor of all hits; it's strongly suggested that this isn't set lower than 
        85.
        |n
        The script has two outputs. Details of all hits not filtered out in the first stage, along 
        with their taxonomies, will be output to a csv if -h/--outhits is specified. If the MEGAN 
        LCA process is used, all hits are output, and the fate of the hit in MEGAN LCA filtering is 
        reported. If any processing is performed, the final taxonomy assigned to each query will be 
        output to -y/--outtaxonomy; if MEGAN LCA processing is performed, summary details of the 
        filtering are also reported. Control the taxonomic ranks output in all cases 
        using -r/--rank. 
        |n
        To retrieve information from NCBI, you must authenticate with an email address and API 
        key. This is only necessary if your local database isn't complete. Supply a file 
        containing this information to -n/--ncbiauth. The script will output an example file if 
        NCBI authentication is required but absent.
        |n
        NCBI requests and parsing by biopython can be buggy. By default, the script will send 
        chunks of 1000 ids to NCBI in each request, and repeat a specific search up to 10 times if 
        any errors are reported. If this doesn't seem to be sufficiently error-tolerant, set 
        -z/--chunksize to a lower value and or -u/--searchreps to a higher value. 
        """, formatter_class=MultilineFormatter)

    # Add individual argument specifications
    parser.add_argument('-b', '--blastresults', required=True, metavar='PATH',
                        help='path to blast results file')
    parser.add_argument('-i', '--minid', type=float, metavar='N', choices=[Range(0,100)],
                        help='filter out hits below this percent identity')
    parser.add_argument('-s', '--minscore', type=float, metavar='N',
                        help='filter out hits below this bitscore')
    parser.add_argument('-l', '--minalen', type=int, metavar='N',
                        help='filter out hits below this BLAST alignment length')
    parser.add_argument('-p', '--process', type=str, metavar='P',
                        choices=('top_score', 'top_id', 'lca', 'megan'),
                        help='if desired, process hits using the given method and return one '
                             'taxonomy per query')
    parser.add_argument('-g', '--gbtiddb', type=str, metavar='PATH', required=True,
                        help='path to accession-taxid json, will be created if absent')
    parser.add_argument('-x', '--tidtaxdb', type=str, metavar='PATH', required=True,
                        help='path to taxid-taxonomy json, will be created if absent')
    parser.add_argument('-n', '--ncbiauth', type=str, metavar='PATH',
                        help='ncbi authentication path')
    parser.add_argument('-c', '--taxidcolumn', type=int, metavar='N',
                        help='the column number of the staxid field if included in an input tsv '
                             'or csv')
    parser.add_argument('-r', '--ranks', metavar='rank,rank,rank',
                        help='comma-separated list of ranks',
                        default='superkingdom,kingdom,phylum,class,order,family,genus,species')
    parser.add_argument('-w', '--winid', type=float, metavar = 'N', choices=[Range(0,100)],
                        help = 'if any hits meet or exceed this percentage identity, only consider '
                               'these hits for MEGAN LCA (no default)')
    parser.add_argument('-e', '--maxexp', type=float, metavar='N', default = 0.01,
                        help="maximum e-value for a hit to be considered in MEGAN LCA")
    parser.add_argument('-t', '--toppc', type = float, metavar='N', default=15, 
                        choices=[Range(0,100)],
                        help="only consider hits with a bitscore within this percentage of the "
                             "highest bitscore for this query in MEGAN LCA")
    parser.add_argument('-a', '--minhitpc', type=float, metavar='N', default = 0, 
                        choices=[Range(0,100)],
                        help="minimum percentage of all hits remaining after filtering for LCA in "
                        "MEGAN LCA")
    parser.add_argument('-d', '--minhitn', type=int, metavar='N', default=1,
                        help="minimum number of hits remaning after filtering for LCA in MEGAN LCA")
    parser.add_argument('-f', '--minsupppc', type = float, metavar='N', default=90, 
                        choices=[Range(0,100)],
                        help='minimum percentage of the hits remaining after filtering that must '
                             'share the final LCA taxonomy')
    parser.add_argument('-o', '--outhits', type=str, metavar='PATH',
                        help="path to write a csv recording the taxonomy of hits")
    parser.add_argument('-y', '--outtaxonomy', type=str, metavar='PATH',
                        help='path to write a csv recording the taxonomy of queries after '
                             'processing')
    parser.add_argument('-z', '--chunksize', type=int, metavar='N', default=1000,
                    help='number of ids per request, default 1000')
    parser.add_argument('-u', '--searchreps', type=int, metavar='N', default=10,
                    help='number of times to repeat a search if errors, default 10')


    # Parse the arguments from the function call if specified, otherwise from the command line
    args = parser.parse_args(arglist) if arglist else parser.parse_args()

    # Do some checking of the inputs

    if args.chunksize > 1000:
        parser.error("-c/--chunksize should not be greater than 1000")
    if args.chunksize < 1:
        parser.error("-c/--chunksize shouldn't be less than 1")
    if args.searchreps < 1:
        parser.error("-u/--searchreps shouldn't be less than 1")

    # Process arguments
    args.ranks = args.ranks.lower().split(',')

    if args.process:
        if not args.outtaxonomy:
            parser.error("if processing taxonomies, supply a path to -y/--outtaxonomy")
        if args.process == "megan":
            args.minid = args.minid if args.minid else 80
            args.minscore = args.minscore if args.minscore else 1
            args.minalen = args.minalen if args.minalen else 100
            if not args.outhits:
                sys.stderr.write(f"Warning: no hit details will be output; supply a path to "
                                 f"-h/--outhits to review hit taxonomies and details\n")
    elif not args.outhits:
        parser.error("supply a path to -h/--outhits")

    # If the arguments are all OK, output them
    return args

# Start the actual script
if __name__ == "__main__":

    # Get the arguments
    args = getcliargs()
    auth = None
    
    # Parse the inputs, filtering out results below idthreshold
    if args.process == 'megan':
        inputdata, gbaccs, taxids = parse_input(args.blastresults, args.taxidcolumn)
    else:
        inputdata, gbaccs, taxids = parse_input(args.blastresults, args.taxidcolumn, args.minscore,
                                                args.minid, args.minalen)

    sys.stderr.write(f"Parsed {args.blastresults}, found {len(inputdata)} BLAST queries, "
                     f"identified {len(gbaccs)} unique GenBank accession values without known "
                     f"taxids, {len(taxids)} known taxids\n")
    
    # If taxids not present, search the unique accession numbers to retreive taxids
    gbtaxids = dict()
    if len(gbaccs) > 0:
        gbtaxids, absent, auth = retrieve_taxids(gbaccs, args.gbtiddb, args.chunksize,
                                            authpath=args.ncbiauth, searchreps = args.searchreps)
        # Add to the master list of taxids
        taxids.update(set(gbtaxids.values()))
    if len(absent) > 0:
        sys.stderr.write(f"Failed to get taxids for {len(absent)} GenBank accession values, these "
                        f"may have been withdrawn from GenBank:\n"
                        f"{','.join(absent)}\n")

    sys.stderr.write(f"Total {len(taxids)} unique taxids to retrieve taxonomy for\n")

    # Retrieve taxonomy from local if available
    taxonomy, absent = retrieve_taxonomy(taxids, args.tidtaxdb, args.chunksize,
                                    authpath=args.ncbiauth, authdict=auth, 
                                    searchreps = args.searchreps)
    if len(absent) > 0:
        sys.stderr.write(f"Failed to get taxids for {len(absent)} taxids, these may have been "
                         f"withdrawn from GenBank:\n"
                         f"{','.join(absent)}\n")

    # Assign taxonomy to the input data
    sys.stderr.write(f"Assigning taxonomy to BLAST hits\n")
    taxonomised = assign_taxonomy(inputdata, gbtaxids, taxonomy, args.ranks)
    
    # Process hits if requested
    if args.process == 'top_score':
        sys.stderr.write(f"Reporting taxonomy of the top hit by bitscore\n")
        taxonomies = filter_tophit(taxonomised, "bitscore")
    elif args.process == 'top_id':
        sys.stderr.write(f"Reporting taxonomy of the top hit by percent identity\n")
        taxonomies = filter_tophit(taxonomised, "pident")
    elif args.process == 'lca':
        sys.stderr.write(f"Performing basic LCA analysis\n")
        taxonomies = lca(taxonomised, args.ranks)
    elif args.process == 'megan':
        sys.stderr.write(f'Performing MEGAN naive LCA analysis\n')
        taxonomies, taxonomised = megan_naive_lca(taxonomised, args.ranks, args.minscore, 
                                                  args.maxexp, args.minid, args.toppc, args.winid,
                                                  args.minhitpc, args.minhitn, args.minalen,
                                                  args.minsupppc)

    # Output
    if args.outhits:
        writehits(taxonomised, args.outhits, args.ranks)
    if args.outtaxonomy:
        writetaxonomy(taxonomies, args.outtaxonomy, args.ranks)

    sys.stderr.write(f"Done\n")
