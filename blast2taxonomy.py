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

import scipy.stats
import random

import numpy as np
import textwrap as _textwrap
from collections import Counter, defaultdict

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
    gbregex = r"(?:^|\s)((?:[A-Z]{1,2}_?)\d{3,}(?:\.\d)?)(?:$|\s)"
    # If this is a longform title with multipl entries, remove all but the first
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
            resd = {'id': sseqid,
                    'data': row,
                    'pident': pident,
                    'evalue': float(row[10]),
                    'bitscore': float(row[11])}
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
            # hit = res.alignments[1]
            sseqid = parse_title(hit.title)
            hsp = hit.hsps[0]
            pident = (hsp.identities/hsp.align_length) * 100
            evalue = hsp.expect
            bitscore = hsp.bits
            score = hsp.score
            if pident >= idthresh:
                hits.append({'id': sseqid,
                             'data': [qseqid, sseqid, pident, evalue, bitscore],
                             'pident': pident,
                             'evalue': evalue,
                             'bitscore': bitscore,
                             'score': score})
        out[qseqid] = hits
    return out


def parse_input(path, process, idthresh, taxidc):
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
    if process == 'top':
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


def retrieve_taxids(gbids, gbtiddbpath, chunksize, authpath=None, authdict=None):
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
            summaries = Entrez.read(sh)
            sh.close()
            return summaries
        sys.stderr.write(f"Searching NCBI nt for {len(absent)} taxids\n")
        ncbigen = retrieve_ncbi_remote(absent, esummary_read_taxids, 'AccessionVersion', chunksize,
                                       authdict)
        rem = {}
        for outsub in ncbigen:
            # outsub = next(ncbigen)
            for gb, smry in outsub.items():
                rem[gb] = int(smry['TaxId'])

        update_ncbi_local(rem, gbtiddbpath)

        out.update(rem)

    absent = set(i for i in gbids if str(i) not in out)

    return out, absent, authdict


def retrieve_taxonomy(taxids, tidtaxdbpath, chunksize, authpath=None, authdict=None):
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
            records = Entrez.read(sh)
            sh.close()
            return records

        sys.stderr.write(f"Searching NCBI taxonomy for {len(absent)} taxonomies\n")
        ncbigen = retrieve_ncbi_remote(absent, efetch_read_taxonomy, 'TaxId', chunksize, authdict)
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

# TODO: implement meganlca option, which filters according to MEGAN algorithm
# See methods section at end https://genome.cshlp.org/content/17/3/377.full
# def filter_hits(data, minscore, toppercent, winscore, minsupport):
#     out = {}
#     for qseqid, hits in data.items():
#


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
        out[qseqid] = [{'taxonomy': lcataxonomy}]
    return out


def p_greater_single(test, series, nperm=10000):
    def firstrank(lis):
        return scipy.stats.rankdata(lis)[0]
    fullseries = [test] + series
    testrank = firstrank(fullseries)
    permranks = []
    for _ in range(nperm):
        random.shuffle(fullseries)
        permranks.append(firstrank(fullseries))
    return len([r for r in permranks if r < testrank])/nperm


def p_greater_pool(test, series, nperm=10000, usestatmean=False):

    def statmax(x, y, axis=0):
        return np.max(x, axis=axis) - np.max(y, axis=axis)

    def statmean(x, y, axis=0):
        return np.mean(x, axis=axis) - np.mean(y, axis=axis)

    res = scipy.stats.permutation_test((test, series),
                                       statmean if usestatmean else statmax,
                                       vectorized=True,
                                       n_resamples=nperm,
                                       alternative='greater')
    return res.pvalue


def perm_test(x, y, stat, alternative, n=5000, converge=None, converge_dp=3):
    # x, y, n, stat, alternative, converge, converge_dp = staxscores, otherscores, 300, 'max', 'greater', 100, 2
    # x, y = [1, 3, 2, 5, 4], [1,2,5,3,4]
    if alternative not in ['greater', 'less']:
        sys.exit("Error: value of alternative in perm_test should be 'greater' or "
                 "'less'")

    def get_random_sample(lis, n):
        new_lis = list(lis)
        random.shuffle(new_lis)
        return new_lis[:n]

    def diffmax(a, b):
        return max(a) - max(b)

    def diffmean(a, b):
        return np.mean(a) - np.mean(b)

    if stat == 'mean':
        fun = diffmean
    elif stat == 'max':
        fun = diffmax
    else:
        sys.exit("Error: value of stat in perm_test should be 'max' or 'mean'")

    def do_perm(x, y, sample, function, alternative):
        xsamp, ysamp = get_random_sample(x, sample), get_random_sample(y, sample)
        exp = fun(xsamp, ysamp)
        complete = list(xsamp + ysamp)
        random.shuffle(complete)
        obs = fun(complete[:sample], complete[sample:])
        return (alternative == 'greater' and exp > obs) or (alternative == 'less' and obs < exp)

    def convergence(x, n, dp):
        if len(x) < n:
            return False
        z = [round(i, dp) for i in x[-n:]]
        return len(set(z)) == 1

    sample = min(len(x), len(y))
    res = 0
    if converge:
        ps = []
        while not (len(ps) >= n and convergence(ps, converge, converge_dp)):
            if do_perm(x, y, sample, fun, alternative):
                res += 1
            ps.append(res/(len(ps)+1))
        return res/len(ps)
    else:
        for _ in range(n):
            if do_perm(x, y, sample, fun, alternative):
                res += 1
        return res/n


def best_hit(data, ranks, speciesid, rejectunknownequal=False, bestp=True, taxp=False, binspecies=False):
    # NOTE: binspecies not currently working
    # data, ranks, process, speciesid, rejectunknownequal, taxp, bestp, binspecies = taxonomised, args.ranks, args.process, args.speciesid, False, False, True, False
    out = {}
    total = len(data)
    for n, (qseqid, hits) in enumerate(data.items()):
        # qseqid, hits = list(data.items())[3]
        sys.stderr.write(f"Performing besthit analysis for query {n+1}/{total}\r")
       # [f"{h['pident']}: {','.join(h['taxonomy'])}" for h in hits]
        # Set up containers for processing and taxonomy reference. If species are present, we want
        # to collapse the hits to a single value per species, to avoid biasing by coverage
        if 'species' in ranks and binspecies:
            speciesi = [i for i, r in enumerate(ranks) if r == 'species'][0]
            taxsets = {r: defaultdict(dict) for r in ranks}
        else:
            speciesi = None
            taxsets = {r: defaultdict(list) for r in ranks}
        parents = dict()
        # Sort the hits by bitscore
        hits = sorted(hits, key=lambda x: x['bitscore'], reverse=True)
        # Parse the hits, assessing each taxonomic level down
        # Store the highest %id currently achieved by a species-level hit. Drop any unknown taxon
        # levels if their %id is not less than this value
        for hit in hits:
            # hit = hits[7]
            lastknowni = 0
            revisedtaxonomy = []
            # First, correct the taxonomy and record parents
            for j, (r, taxon) in enumerate(zip(ranks, hit['taxonomy'])):
                # Retrieve key details
                parent = None if j == 0 else revisedtaxonomy[j - 1]
                # If the taxon is blank, or in the special case that this is a species and it's an
                # undetermined species
                if taxon == '' or (r == 'species' and re.search(r" (sp|aff)\.?( |$)", taxon))\
                        or (r == 'genus' and re.search(r" gen\.?( |$)", taxon)):
                    # Check the last known taxon and see if there's an existing hit
                    lastknown = '' if lastknowni == 0 else f"{hit['taxonomy'][lastknowni]}_"
                    taxon = f"unknown_{lastknown}{r}"
                else:
                    lastknowni = j
                # Record the corrected taxonomy and output
                revisedtaxonomy.append(taxon)
                parents[taxon] = parent
            # Second, check and score
            pident = hit['pident']
            species = revisedtaxonomy[speciesi] if 'species' in ranks and binspecies else None
            for j, (r, taxon) in enumerate(zip(ranks, revisedtaxonomy)):
                # j = 6; r = ranks[j]; taxon = hit['taxonomy'][j]
                # If we're at species level, we don't want to consider the species if it's not a
                # high enough %id to be likely to be an exact species hit
                if r == 'species' and pident < speciesid:
                    break
                # If the taxon is blank, or in the special case that this is a species and it's an
                # undetermined species
                if 'unknown_' in taxon:
                    # Check the last known taxon and see if there's an existing hit
                    if j > 0:
                        lastknown = revisedtaxonomy[lastknowni]
                        if lastknown in taxsets[ranks[j-1]]:
                            if species:
                                currbest = max(taxsets[ranks[j - 1]][lastknown].values())
                            else:
                                currbest = max(taxsets[ranks[j - 1]][lastknown])
                            if currbest == 100 or currbest > pident:
                                # If the current best for this taxon is a 100 match, this is either
                                # an unidentified conspecific (if pident == 100), or it will just
                                # be rejected in favour - either way, reject this. Similarly if the
                                # current best for this taxon is better than this pident, this
                                # hit has already contributed its pident as far as possible so
                                # reject as adding unknowns will just muddy the water
                                break
                            elif currbest == pident:
                                # However in the unique situation where the current best has the
                                # same identity as this unknown, leave it to the user. Chances are
                                # the unknown is the same species as the identificied ones, but
                                # if we're conservitive, we should simply say unknown at this level
                                if rejectunknownequal:
                                    break
                else:
                    lastknowni = j
                # Calculate the score
                if species:
                    if r == 'species':
                        species = taxon
                    if species in taxsets[r][taxon]:
                        taxsets[r][taxon][species] = max(pident, taxsets[r][taxon][species])
                    else:
                        taxsets[r][taxon][species] = pident
                else:
                    taxsets[r][taxon].append(pident)

        # If species are present, collect the %id by species together into a list
        if 'species' in ranks and binspecies:
            for r in ranks:
                taxa = list(taxsets[r].keys())
                for t in taxa:
                    taxsets[r][t] = list(taxsets[r][t].values())
        # Reorganise the parent dict to a children dict - done after the above to ensure that any
        # unconsidered species don't show up
        children = defaultdict(set)
        for child, parent in parents.items():
            if parent:
                children[parent].add(child)
        # Output containers
        bestfittaxonomy = []
        bestfitpvalue = []
        pident = None
        for i, r in enumerate(ranks):
            # i = 0 ; r = ranks[i]
            # If this is not the root, remove any other taxa that are not valid children
            possibletaxa = list(taxsets[r].keys())
            if i > 0:
                for t in possibletaxa:
                    if t not in children[bestfittaxonomy[i - 1]]:
                        del taxsets[r][t]
            # If no hits remain, or the remaining hit is unknown, fill the rest of the output lists
            # with nothing and break
            if len(taxsets[r]) == 0 or all('unknown_' in t for t in taxsets[r].keys()):
                break
            # If only one hit remains, return it with the probability of its parent, or 1 if this
            # is the root
            elif len(taxsets[r]) == 1:
                pident = max(list(taxsets[r].values())[0])
                bestfittaxonomy.append(list(taxsets[r].keys())[0])
                bestfitpvalue.append(1.0 if i == 0 else bestfitpvalue[-1])
                continue
            # Otherwise, calculate the summary scores
            scores = {'taxa': [], 'max': [], 'mean': []}
            for t in taxsets[r].keys():
                scores['taxa'].append(t)
                # Sum: good for getting consensus, but too overwhelmed by coverage
                #lcasets[r][t] = sum(lcasets[r][t])
                # Mean: undersupports mid level taxa
                #lcasets[r][t] = sum(lcasets[r][t])/len(lcasets[r][t])
                # Max: isn't biased by coverage, nor by midlevel taxa
                scores['max'].append(max(taxsets[r][t]))
                scores['mean'].append(np.mean(taxsets[r][t]))

            # Find the taxon with the best score
            besttaxa = [t for t, x in zip(scores['taxa'], scores['max']) if x == max(scores['max'])]
            # Break ties using mean
            if len(besttaxa) > 1:
                besttaxa = [t for t, m in zip(scores['taxa'], scores['mean'])
                            if m == max(scores['mean']) and t in besttaxa]
                # If ties remain, or the best taxon is unknown, output nothing
                if len(besttaxa) > 1 or 'unknown_' in besttaxa[0]:
                    break
            # Calculate the probability that we would select this taxon based on its score if
            # we had the same number of hits as the next best taxa do, or vice versa
            stax = besttaxa[0]
            staxscores = taxsets[r][stax]
            otherscores = []
            for t, s in taxsets[r].items():
                if t != stax and 'unknown_' not in t:
                    otherscores.extend(s)
            pvalue = 1.0 if len(otherscores) == 1 else perm_test(staxscores, otherscores, 'max',
                                                                 'greater', n=1000, converge=300)
            pident = max(staxscores)
            # If this is the root, output the taxon and p value
            bestfittaxonomy.append(stax)
            if i == 0:
                bestfitpvalue.append(pvalue)
            else:
                # Otherwise, calculate the probability that this taxon is maximum and that its parental
                # taxa are also maximum
                bestfitpvalue.append(bestfitpvalue[i-1] * pvalue)

        if taxp:
            outtaxonomy = [f"{t}_{str(round(s, 2))}" for t, s in zip(bestfittaxonomy, bestfitpvalue)]
        else:
            outtaxonomy = bestfittaxonomy
        outtaxonomy += ([''] * (len(ranks) - len(outtaxonomy)))
        out[qseqid] = [{'taxonomy': outtaxonomy}]
        if bestp:
            out[qseqid][0]['pident'] = pident
    sys.stderr.write("\nCompleted besthit analysis\n")
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
                output.append(hit['tx'])
            if 'pident' in hit and 'data' not in hit:
                output.append(hit['pident'])
            output.extend(hit['taxonomy'])
            output = [str(o) for o in output]
            sys.stdout.write('\t'.join(output) + '\n')


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
        The taxonomy will be retrieved from NCBI Taxonomy based on the taxid of the hit(s). By 
        default, the taxonomy of all hits will be retrieved and reported. Optionally, the hits can 
        be processed, using the option -p/--process, either keeping only the/a top hit using 
        "-p top", finding the lowest common ancestor of all hits using "-p lca", or attempting to 
        find the best fit taxonomy using -p "bestfit", which will consider the hit properties to 
        select a single taxon for each rank and, if -w/--showscores is set, report the score (0-1) 
        for that taxon. Note that scores are currently not calculated well where multiple higher 
        taxa have high percent ID hits.
        |n
        The script can filter hits based on percent identity in two ways. To speed up NCBI 
        searches, any hits below a minimum percent identity supplied to -m/--minidthreshold will be 
        discarded prior to taxonomy assignment, default 0. When using the best fit option only, 
        species level taxonomic assignment will only consider hits with a percent identity greater 
        than or equal to the value of -s/--speciesidthreshold, default 95. By default, during best 
        fit processing, hits with unknown species identity are excluded only if they score worse
        than hits with known identity, the conservative option. To reject unknown hits that score 
        worse or equal to hits with known identity, use the flag '-j/--rejectunknowns. This is 
        only recommended when it is expected that NCBI contains the species in the query.
        |n
        The script will output a tsv to STDOUT. With the default option, the script will output 
        the information each hit for each query id, in the same format as supplied, followed by 
        the parsed NCBI accession ID, the NCBI taxid and taxonomy columns for each hit. When using 
        the -t/--tophit option, only a top hit for each query will be output. If multiple top 
        hits are found, the first one with the most frequent taxonomy will be output. When using 
        the -l/--lca option, the script will output a single row for each query id and the LCA  
        taxonomy columns. Note that if the input is XML, the only fields parsed from the XML and 
        output are qseqid, sseqid and pident. Control the taxonomic ranks output in all cases 
        using -r/--rank. 
        |n
        To efficiently retrieve taxonomy, the script uses a pair of local datasets in json files, 
        one recording taxids for GenBank accession numbers, one recording taxonomy for taxids. 
        These might have been created by a previous run of this script, or the latter with 
        get_NCBI_taxonomy.py. Supply a path to the accession-taxid json with -g/--gbtiddb and the 
        taxid-taxonomy json with -x/--tidtaxdb. If this is your first run, the script will create 
        new databases at the locations supplied.
        |n
        To retrieve information from NCBI, you must authenticate with an email address and API 
        key. This is only necessary if your local database isn't complete. Supply a file 
        containing this information to -n/--ncbiauth. The script will output an example file if 
        NCBI authentication is required but absent.
        |n
        By default, the script will send 1000 ids to NCBI in each request. If this seems to cause 
        errors, set -z/--chunksize to a lower value. This cannot be increased.
        """, formatter_class=MultilineFormatter)

    # Add individual argument specifications
    parser.add_argument('-b', '--blastresults', required=True, metavar='PATH',
                        help='path to blast results file')
    parser.add_argument('-p', '--process', type=str, metavar='P',
                        choices=('top', 'lca', 'bestfit'),
                        help='if desired, the name of the method to process the hits to return '
                             'one taxonomy per query')
    parser.add_argument('-g', '--gbtiddb', type=str, metavar='PATH', required=True,
                        help='path to accession-taxid json, will be created if absent')
    parser.add_argument('-x', '--tidtaxdb', type=str, metavar='PATH', required=True,
                        help='path to taxid-taxonomy json, will be created if absent')
    parser.add_argument('-s', '--speciesid', type=float, metavar='N',
                        help="only consider species from hits above this percent identity during "
                             "best fit analysis")
    parser.add_argument('-j', '--rejectunknowns', action='store_true',
                        help="do not consider unknown taxonomies with the same score as known "
                             "taxonomies during best fit analysis")
    parser.add_argument('-w', '--showscores', action='store_true',
                        help="show individual taxon scores in the best hit analysis output")
    parser.add_argument('-n', '--ncbiauth', type=str, metavar='PATH',
                        help='ncbi authentication path')
    parser.add_argument('-m', '--minid', type=float, metavar='N', default=0,
                        help='minimum percent id to retrieve taxonomy for a hit, default 0')
    parser.add_argument('-c', '--taxidcolumn', type=int, metavar='N',
                        help='the column number of the staxid field if included in an input tsv '
                             'or csv')
    parser.add_argument('-r', '--ranks', metavar='rank,rank,rank',
                        help='comma-separated list of ranks',
                        default='superkingdom,kingdom,phylum,class,order,family,genus,species')
    parser.add_argument('-z', '--chunksize', type=int, metavar='N', default=1000,
                        help='number of ids per request, default 1000')

    # Parse the arguments from the function call if specified, otherwise from the command line
    args = parser.parse_args(arglist) if arglist else parser.parse_args()

    # Do some checking of the inputs

    if args.process == 'bestfit':
        args.speciesid = args.speciesid if args.speciesid else 95
    else:
        msg = "is only relevant when using '-p bestfit'/'--process bestfit'"
        if args.speciesid:
            parser.error(f"Error: -s/--speciesid {msg}")
        if args.rejectunknowns:
            parser.error(f"Error: -j/--rejectunknowns {msg}")
        if args.showscores:
            parser.error(f"Error: -w/--showscores {msg}")
    if args.chunksize > 1000:
        parser.error("-c/--chunksize should not be greater than 1000")

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
    inputdata, gbaccs, taxids = parse_input(args.blastresults, args.process, args.minid,
                                            args.taxidcolumn)

    sys.stderr.write(f"Parsed {args.blastresults}, found {len(inputdata)} BLAST queries, "
                     f"identified {len(gbaccs)} unique GenBank accession values without known "
                     f"taxids, {len(taxids)} known taxids\n")
    # If taxids not present, search the unique accession numbers to retreive taxids
    gbtaxids = dict()
    if len(gbaccs) > 0:
        gbtaxids, _, auth = retrieve_taxids(gbaccs, args.gbtiddb, args.chunksize,
                                            authpath=args.ncbiauth)
        # Add to the master list of taxids
        taxids.update(set(gbtaxids.values()))

    sys.stderr.write(f"Total {len(taxids)} unique taxids to retrieve taxonomy for\n")

    # Retrieve taxonomy from local if available
    taxonomy, _ = retrieve_taxonomy(taxids, args.tidtaxdb, args.chunksize,
                                    authpath=args.ncbiauth, authdict=auth)

    # Assign taxonomy to the input data
    sys.stderr.write(f"Assigning taxonomy to BLAST hits\n")
    taxonomised = assign_taxonomy(inputdata, gbtaxids, taxonomy, args.ranks)

    # Process hits if requested
    if args.process == 'top':
        sys.stderr.write(f"Filtering hits to report top hit taxonomy\n")
        taxonomised = filter_tophit(taxonomised)
    elif args.process == 'lca':
        sys.stderr.write(f"Performing LCA analysis\n")
        taxonomised = lca(taxonomised, args.ranks)
    elif args.process == 'bestfit':
        #data, ranks, speciesid, rejectunknownequal = False, bestp = True, taxp = False, binspecies = False
        taxonomised = best_hit(taxonomised, args.ranks, args.speciesid,
                               args.rejectunknowns, True, args.showscores)

    # Output
    writeout(taxonomised)
    sys.stderr.write(f"Done\n")
