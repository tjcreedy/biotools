#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Description of what the script does"""

# Imports
import os
import sys
import shutil
import argparse
import functools
import multiprocessing
import subprocess
import time
import datetime
import re
import urllib
import tarfile

import textwrap as _textwrap

from io import StringIO
from collections import defaultdict, Counter
import Bio.Seq
from Bio import SeqIO, SeqFeature
from Bio.Data import IUPACData
from Bio.Align.Applications import MafftCommandline

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


class NullQueue:
    def __init__(self, ):
        pass

    def put(self, value):
        print(value)

    def get(self):
        pass


# Monkey Patches


Bio.SeqFeature.FeatureLocation.__hash__ = lambda self : hash(f"{self.start.__str__()}"
                                                             f"{self.end.__str__()}"
                                                             f"{self.strand}")


# Global variables

def loadnamevariants(source=None, simplifytrnas=False):

    conversion = {}
    url = "https://raw.githubusercontent.com/tjcreedy/constants/master/gene_name_variants.txt"
    fullparse = {}
    alltypes = set()
    #  Identify source

    if source is None:
        source = urllib.request.urlopen(url)
    else:
        source = open(source, 'r')

    # Read source

    for line in source:
        #line = list(source)[36]
        line = line.decode('utf-8').strip()
        description, variants = line.split(":")
        name, annotype, fullname = description.split(";")
        variants = variants.split(',')

        outname, fulloutname = name, fullname
        if simplifytrnas and annotype == 'tRNA':
            outname = name[:4]
            fulloutname = fullname[:8]

        fullvariants = set()
        for v in set([name, outname] + variants + [fullname.upper(), fulloutname.upper()]):
            for g in ['', ' ']:
                v = v.replace(g, '')
                for s in ['', ' GENE', ' ' + annotype.upper()]:
                    fullvariants.add(v + s)
                    conversion[v + s] = outname

        alltypes.add(annotype)
        if outname not in fullparse:
            fullparse[outname] = {'type': {annotype},
                                  'variants': fullvariants,
                                  'product': {fulloutname}}
        else:
            fullparse[outname]['type'].add(annotype)
            fullparse[outname]['variants'].update(fullvariants)
            fullparse[outname]['product'].add(fulloutname)


    for name in fullparse.keys():
        for i in ['type', 'product']:
            fullparse[name][i] = list(fullparse[name][i])[0]

    # Close handle
    source.close()
    return conversion, alltypes, fullparse

nameconvert, alltypes, variants = loadnamevariants()

formatregex = {'fasta': "^>",
               'genbank': "^LOCUS"}

# Function definitions


def get_feat_name(feat):
    featname = "unknown"
    nametags = ['gene', 'product', 'label', 'standard_name']
    if any(t in feat.qualifiers.keys() for t in nametags):
        for t in nametags:
            if t in feat.qualifiers.keys():
                featname = feat.qualifiers[t][0].upper()
                break
    return featname


def set_feat_name(feat, name):
    nametags = ['gene', 'product', 'label', 'standard_name']
    if any(t in feat.qualifiers.keys() for t in nametags):
        for t in nametags:
            if t in feat.qualifiers.keys() :
                feat.qualifiers[t][0] = name
    return feat


def write_log(logq):
    # Set up the logging output handle
    dt = datetime.datetime.now()
    logh = open(f"fixtrnas_log_{dt.strftime('%Y-%m-%d-%H%M%S')}.txt", 'w')
    # Start a constant process that waits to receive data in the form of a file
    # specifying the source of the seqrecord, and the new seqrecord
    while 1:
        logline = logq.get()
        if logline is None:
            break
        now = datetime.datetime.now()
        logh.write(f"{str(now-dt)}\t{logline}")
        logh.flush()
    logh.close()


def write_genbank(seqq):
    """Sets up the output file writing then waits to receive data
       requesting printing of log info and seqrecord to the specified file"""

    # Start a constant process that waits to receive data in the form of a file
    # specifying the source of the seqrecord, and the new seqrecord
    while 1:
        # queueitem = (outname, seqrecord, filetotal)
        queueitem = seqq.get()
        if queueitem is None:
            break
        seqrecord = queueitem
        # Check if the file handle is already open for this file
        sys.stdout.write(seqrecord.format("genbank"))
        sys.stdout.flush()


def write_term(prinq):
    # filenames = args.genbank
    # Set up to print
    start = time.perf_counter()
    done = 0
    # Print
    while 1:
        sys.stderr.write(f"\r{' ' * 70}\r")
        line = f"\rProcessed {done} records"
        sys.stderr.write(line)
        sys.stderr.flush()
        queueitem = prinq.get()
        if queueitem is None:
            break
        done += 1
    now = time.perf_counter()
    elapsed = now - start
    elapsedper = datetime.timedelta(seconds=elapsed / done)
    elapsed = datetime.timedelta(seconds=round(elapsed))
    line = f"\nFinished in {elapsed}, {elapsedper} per record\n\n"
    sys.stderr.write(line)
    sys.stderr.flush()


def start_accessories(pool, manager):

    if pool is None or manager is None:
        sys.stderr.write("start_accessories has been passed null pool and/or manager: assuming "
                         "interactive\n")
        return (NullQueue(), NullQueue(), NullQueue()), None
    else:
        seqq = manager.Queue()
        logq = manager.Queue()
        termq = manager.Queue()

        seqwrite = pool.apply_async(functools.partial(write_genbank), (seqq,))
        logwrite = pool.apply_async(functools.partial(write_log), (logq,))
        termwrite = pool.apply_async(functools.partial(write_term), (termq,))

        return (seqq, logq, termq), (seqwrite, logwrite, termwrite)


def sequence_generator(current, handle, seqformat):
    for nextline in handle:
        # If next line is the start of a new record, output the current record
        if re.match(formatregex[seqformat], nextline):
            yield next(SeqIO.parse(StringIO(current), seqformat))
            current = nextline
        # Otherwise keep building the current record
        else:
            current += nextline
    # Once out of lines, output the current record
    yield next(SeqIO.parse(StringIO(current), seqformat))


def parse_input(path=None):
    # Get handle
    source = open(path, 'r') if path else sys.stdin
    # Find out format from first line
    firstline = source.readline()
    seqformat = [k for k, r in formatregex.items() if re.match(r, firstline)][0]
    # Return generator and boolean for whether checking should be done or not
    return sequence_generator(firstline, source, seqformat), seqformat == "genbank"


def parse_reference(path, map):
    # Parse map if present
    mapdict = defaultdict(list)
    if map:
        with open(map, 'r') as mh:
            for line in mh:
                ref, inseq = re.split("[\t,]", line.rstrip())
                mapdict[ref].append(inseq)

    # Parse reference and create dict:
    refdict = {}
    refrecords = SeqIO.parse(path, 'genbank')
    for refseq in refrecords:
        dictnames = [refseq.name]
        if map:
            if refseq.name not in mapdict:
                sys.exit(f"Error: reference sequence {refseq.name} is not in {map}")
            dictnames = mapdict[ref]
        for n in dictnames:
            refdict[n] = refseq

    return(refdict)


def gapped_position(seq, n):
    """Return the position in the gapped sequence corresponding to the given position in the
    ungapped sequence. If the given position exceeds the number of ungapped bases, returns
    None"""
    i = -1
    for c, b in enumerate(seq):
        if b != '-':
            i += 1
        if i == n:
            return c
    return None


def ungapped_position(seq, n):
    """Return the position in the ungapped sequence corresponding to the given alignment
    position. If the alignment position corresponds to a gap in the gapped sequence, returns
    None"""
    if seq[n] == '-':
        return None
    else:
        return n - seq[:n].count('-')


def sortextractfeats(record, extracttypes):
    if type(extracttypes) is str:
        extracttypes = {extracttypes,}
    focalfeats = {}
    otherfeats = []
    for feat in record.features:
        convtype = ''
        if feat.type == 'gene':
            tryname = get_feat_name(feat)
            if tryname in nameconvert:
                convtype = variants[nameconvert[tryname]]['type']
            elif tryname[:3].upper() == 'TRN':
                convtype = 'tRNA'
        if feat.type in extracttypes or convtype in extracttypes:
            if feat.location in focalfeats:
                focalfeats[feat.location].append(feat)
            else:
                focalfeats[feat.location] = [feat]
        else:
            otherfeats.append(feat)
    return focalfeats, otherfeats


def closestbases(seq, n, nbases=1):
        segments = {'left': seq[:n][::-1],
                    'right': seq[n:]}
        cbases = defaultdict(list)
        for i in range(0, nbases):
            # i = list(range(0, nbases))[0]
            for dir in segments.keys():
                cbases[dir].append(gapped_position(segments[dir], i))
        return([None if i is None else n-1-i for i in cbases['left']],
               [None if i is None else n+i for i in cbases['right']])


def get_simple_trnas(seqrecord, nameconvertsimple):
    seqtrnas = defaultdict(list)
    seqtrnasbyloc, seqother = sortextractfeats(seqrecord, 'tRNA')
    # Sort tRNAs by name to remove as new ones found
    for loc, feats in seqtrnasbyloc.items():
        # loc, feats = list(seqtrnasbyloc.items())[0]
        # feat = feats[0]
        names = {nameconvertsimple[get_feat_name(feat)] for feat in feats}
        if len(names) > 1:
            seqtrnas['unknown'].extend(feats)
        else:
            seqtrnas[list(names)[0]].extend(feats)
    return seqtrnas, seqother


def unreplaced_trnas(seqtrnas, trnasadded):
    outfeats = []
    for simplename, origfeats in seqtrnas.items():
        if ((simplename in ['TRNL', 'TRNS'] and trnasadded[simplename] == 4) or
                (simplename not in ['TRNL', 'TRNS'] and trnasadded[simplename] == 2)):
            pass
        else:
            outfeats.extend(origfeats)
    return outfeats


def ref_annotate(seqrecord, refdict, tempdir, queues):
    logq, termq = queues

    # Get the reference sequence from refdict
    if len(refdict) > 1:
        refname = {i for i in [seqrecord.name, seqrecord.id] if i in refdict}
        if len(refname) == 0:
            sys.exit(f"Error: input sequence {seqrecord.name} is not in supplied reference file ("
                     f"or reference map file, if supplied)")
        elif len(refname) > 1:
            sys.exit(f"Error: input sequence name {seqrecord.name} and ID {seqrecord.id} match "
                     f"different to different reference (or map) names")
        refrecord = refdict[list(refname)[0]]
    else:
        refrecord = list(refdict.values())[0]

    # If seqrecord has features, extract out tRNA and other feats from seqrecord
    seqtrnas = {}
    nameconvertsimple, *_ = loadnamevariants(simplifytrnas=True)
    if len(seqrecord.features) > 0:
        seqtrnas, outfeats = get_simple_trnas(seqrecord, nameconvertsimple)
    else:
        outfeats = []

    # If sequences match, no need to do an alignment
    reftrnas, refother = sortextractfeats(refrecord, 'tRNA')
    align = True
    refaln, seqaln = None, None
    if seqrecord.seq == refrecord.seq:
        align = False
        logq.put(f"{seqrecord.name}: reference annotation - sequences are identical")
    else:
        # Otherwise, run alignment
        fastaname = os.path.join(tempdir, f"{seqrecord.name}.fasta")
        with open(fastaname, 'w') as fh:
            for rec in [seqrecord, refrecord]:
                fh.write(f">{rec.name}\n{rec.seq}\n")
        mafftcli = MafftCommandline(input=fastaname, adjustdirection=True)
        stdout, stderr = mafftcli()
        seqaln, refaln = list(SeqIO.parse(StringIO(stdout), 'fasta'))
        os.unlink(fastaname)
        # If the target is RC, reverse it back and reverse the references to match
        if re.match('^_R_', seqaln.name):
            seqaln = seqaln.reverse_complement()
            refaln = refaln.reverse_complement()
            refrecord = refrecord.reverse_complement()
        # If the reference is RC, reverse the reference to match
        elif re.match('^_R_', refaln.name):
            refrecord = refrecord.reverse_complement()
        reftrnas, refother = sortextractfeats(refrecord, 'tRNA')
        logq.put(f"{seqrecord.name}: reference annotation - completed aligning sequences")

    replaced, unreplaced = 0, 0
    for location, feats in reftrnas.items():
        trnasadded = Counter()
        # If aligned, modify the location of the feats based on the alignment
        if align:
            # location, feats = list(reftrnas.items())[0]
            # Get the feature length
            flen = len(location)

            # Convert positions on one sequence to positions on the other using alignment
            alnpos = [gapped_position(refaln.seq, int(location.start)),
                      gapped_position(refaln.seq, int(location.end))]
            seqpos = [ungapped_position(seqaln.seq, ap) for ap in alnpos]

            # Check that both positions have matched up and generate approximate locations if not
            outpos = [None, None]
            trunc = [type(location.start) != SeqFeature.ExactPosition,
                     type(location.end) != SeqFeature.ExactPosition]
            for pos, sp, ap in zip([0, 1], seqpos, alnpos):
                # pos, sp, ap = 0, seqpos[0], alnpos[0]
                if sp is not None:
                    outpos[pos] = sp
                    continue
                # Get the closest bases in each direction
                clspos = closestbases(seqaln.seq, ap, nbases=1)
                # Get an  approxmiate position
                apxalnpos = clspos[pos][0]
                # Check that there are positions available, if not set as truncated
                if apxalnpos is None:
                    apxalnpos = clspos[1 if pos == 0 else 0][0]
                    trunc[pos] = True
                outpos[pos] = ungapped_position(seqaln.seq, apxalnpos)
            # If either position has not matched up, shift the position as needed
            if any(p is None for p in seqpos):
                # Get the feature length
                flen = len(location.extract(refrecord))
                # Calculate an approximate length
                apxlen = outpos[1] - outpos[0]
                # Check the approximate difference in length
                # +ve = new too short, -ve = new too long
                lendiff = flen - apxlen
                # If not zero, nudge the positions as needed
                if seqpos == [None, None]:
                    # If both ends were none, expand or contract equally, unless truncated
                    outpos[0] = outpos[0] - (round(lendiff/2) if not trunc[0] else 0)
                    outpos[1] = outpos[1] + (round(lendiff/2) if not trunc[1] else 0)
                    # Otherwise expand or contract the None end only
                elif seqpos[0] is None and not trunc[0]:
                    outpos[0] = outpos[0] - lendiff
                elif seqpos[1] is None and not trunc[1]:
                    outpos[1] = outpos[0] + lendiff
                # Check that the sequence length hasn't been exceeded
                if outpos[0] < 0:
                    outpos[0] = 0
                    trunc[0] = True
                if outpos[1] >= len(seqrecord):
                    outpos[1] = len(seqrecord)-1
                    trunc[1] = True

            # Create the new location with the correct type of position
            outpost = []
            for op, t, ep, tp in zip(outpos,
                                     trunc,
                                     [SeqFeature.ExactPosition] * 2,
                                     [SeqFeature.BeforePosition, SeqFeature.AfterPosition]):
                outpost.append(tp(op) if t else ep(op))

            outl = SeqFeature.FeatureLocation(start=outpost[0], end=outpost[1],
                                              strand=location.strand)
            # Change the location on the features and append them to seqrecord
            outfeats = []
            for feat in feats:
                feat.location = outl
                outfeats.append(feat)
            feats = outfeats

        # Add the features using the standard name, record simplified names
        for feat in feats:
            name = get_feat_name(feat)
            if name in nameconvert:
                feat = set_feat_name(feat, nameconvert[name])
            if name in nameconvertsimple:
                trnasadded[nameconvertsimple[name]] += 1
            outfeats.append(feat)
        # Add back in any tRNAs from original set not found
        unrepfeats = unreplaced_trnas(seqtrnas, trnasadded)
        outfeats.extend(unrepfeats)
        if feats:
            replaced += 1
        if unrepfeats:
            unreplaced += 1
    logq.put(f"{seqrecord.name}: reference annotation - {replaced} tRNA annotations "
             f"replaced/added, {unreplaced} not replaced")
    issues = None
    seqrecord.features = outfeats
    return seqrecord, issues


def get_mitfi(tempdir):
    # Download mitfi into the temporary
    source = "https://raw.githubusercontent.com/tjcreedy/biotools/main/bin/mitfi.tar.gz"
    destination = os.path.join(tempdir, 'mitfi.tar.gz')
    urllib.request.urlretrieve(source, destination)
    tar = tarfile.open(destination, "r:gz")
    tar.extractall(tempdir)
    tar.close()

def mitfi_annotate(seqrecord, table, tempdir, queues):
    # table = args.table
    logq, termq = queues

    # If seqrecord has features, extract out tRNA and other feats from seqrecord
    seqtrnas = {}
    nameconvertsimple, *_ = loadnamevariants(simplifytrnas=True)
    if len(seqrecord.features) > 0:
        seqtrnas, outfeats = get_simple_trnas(seqrecord, nameconvertsimple)
    else:
        outfeats = []

    # Output the sequence as a fasta
    fastapath = os.path.join(tempdir, f"{seqrecord.name}.fasta")
    with open(fastapath, 'w') as fh:
        fh.write(f">{seqrecord.name}\n{seqrecord.seq}\n")

    # Generate a Mitfi config and output
    logq.put(f"{seqrecord.name}: MiTFi annotation starting")
    configpath = os.path.join(tempdir, 'mitfi_config.txt')
    with open(configpath, 'w') as fh:
        fh.write(f"infernalCall = {os.path.join(tempdir, 'cmsearch')}\n"
                 f"mpiCall = orterun\n"
                 f"cores = 1\n"
                 f"CodeFile = gc.prt\n"
                 f"code = {table}\n"
                 f"evalue = 0.001\n"
                 f"overlap = 10\n")

    # Run Mitfi
    cmd = f"java -Xmx2048m " \
          f"-jar {os.path.join(tempdir, 'mitfi.jar')} " \
          f"-onlycutoff " \
          f"{fastapath}".split(' ')
    mitfirun = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Parse output
    trnasadded = Counter()
    for line in StringIO(mitfirun.stdout.decode('UTF-8')):
        # line = list(StringIO(mitfirun.stdout.decode('UTF-8')))[4]
        if line[0] == "#":
            continue
        # Parse line
        n, start, end, score, evalue, AC, AA, model, strand = line.rstrip().split('\t')
        anticodon = AC.replace('U', 'T')
        shortname = f"TRN{AA[0]}"
        name = f"{shortname}-{anticodon}"
        trnasadded[shortname] += 2
        # Create features
        location = SeqFeature.FeatureLocation(start=SeqFeature.ExactPosition(int(start)),
                                              end=SeqFeature.ExactPosition(int(end)),
                                              strand=1 if strand == "+" else -1)
        outfeats.extend([
            SeqFeature.SeqFeature(location, type='gene', qualifiers={'gene': name}),
            SeqFeature.SeqFeature(location, type='tRNA', qualifiers={'gene': name,
                                                              'product': variants[name]['product'],
                                                              'standard_name': name})])

    # Add back in any tRNAs from original set not found
    unrepfeats = unreplaced_trnas(seqtrnas, trnasadded)
    outfeats.extend(unrepfeats)
    unreplaced = set(f.location for f in unrepfeats)

    logq.put(f"{seqrecord.name}: MiTFi annotation - {len(trnasadded)} tRNA annotations "
             f"replaced/added, {len(unreplaced)} not replaced")

    issues = None
    seqrecord.features = outfeats
    return seqrecord, issues


def add_anticodon(feats):
    # feats = trnas
    outfeats = []
    noanticodons = 0
    # Iterate
    for i, featlis in enumerate(feats.values()):
        # featlis = list(feats.values())[0]
        # Check if at least one is a tRNA, dicard all if not
        types = set(f.type for f in featlis)
        if 'tRNA' not in types:
            outfeats.extend(featlis)
            continue
        # Get the names of the features
        names = [get_feat_name(feat) for feat in featlis]
        trnanames = [get_feat_name(feat) for feat in featlis if feat.type == 'tRNA']
        #sys.stderr.write(f"{i} {names}\n")
        # Check if the names are already recognisable, discard any that are, discard any that don't
        # have the same name as any trnas
        featcont = []
        for feat, name in zip(featlis, names):
            if name in nameconvert or name not in trnanames:
                outfeats.append(feat)
            else:
                featcont.append(feat)
        if len(featcont) == 0:
            continue
        # Search for anticodon tags
        anticodons = {}
        for feat in featcont:
            # feat = featcont[0]
            if 'anticodon' in feat.qualifiers:
                name = get_feat_name(feat)
                # Check for tags
                acstring = feat.qualifiers['anticodon']
                if len(acstring) > 1:
                    sys.stderr.write(f"Warning: {feat.type} annotation {name} has multiple "
                                     f"anticodon tags, this will be skipped\n")
                    continue
                # Search for information within tag
                acstring = acstring[0]
                parser = {#'position': r"pos:([0-9]+)\.\.([0-9]+)",
                          'aminoacid': r"aa:([A-Za-z]+)",
                          'anticodon': r"seq:([A-Za-z]{3})"}
                acdata = {}
                for part, regex in parser.items():
                    m = re.search(regex, acstring)
                    if not m:
                        sys.stderr.write(f"Warning: {part} cannot be found in {feat.type} "
                                         f"annotation {name} anticodon tag, this will be skipped\n")
                        continue
                    acdata[part] = m.groups(1) if len(m.groups(1)) > 1 else m.groups(1)[0]
                if len(acdata) < 2:
                    continue
                # Construct putative new names
                protletter = IUPACData.protein_letters_3to1[acdata['aminoacid']]
                acname = f"TRN{protletter}-{acdata['anticodon'].upper()}"
                if name in anticodons:
                    anticodons[name].add(acname)
                else:
                    anticodons[name] = {acname}
        if len(anticodons) == 0:
            noanticodons += 1
            outfeats.extend(featcont)
            continue
        # Set up renamer
        rename = {}
        for fname, acname in anticodons.items():
            if len(acname) > 1:
                sys.stderr.write(f"Warning: anticodon tags do not match for the two or more "
                                 f"annotations named {fname} at the same locus")
                outfeats.extend(featcont)
                continue
            rename[fname] = list(acname)[0]
        # Rename all
        outfeats.extend([set_feat_name(feat, rename[get_feat_name(feat)]) for feat in featcont])
    # Output list of features
    return outfeats, noanticodons


def do_checking(seqrecord, queue):
    logq, termq = queue
    # Extract tRNAs
    trnas, otherfeats = sortextractfeats(seqrecord, 'tRNA')
    logq.put(f"{seqrecord.name}: checking - found {len(trnas)} tRNA locations")
    # Add anticodons
    correctedtrnas, nac = add_anticodon(trnas)
    logq.put(f"{seqrecord.name}: checking - {nac} tRNAs had no anticodon data")
    # Return all features to seqrecord
    seqrecord.features = correctedtrnas + otherfeats

    issues = (len(trnas) < len(otherfeats), nac)
    logq.put(f"{seqrecord.name}: checking - {len(otherfeats)} other features")

    return seqrecord, issues


def process_seqrecord(dochecking, args, refdict, tempdir, queues, seqrecord):
    # seqrecordgen, dochecking = parse_input("/home/thomas/work/iBioGen_postdoc/MMGdatabase/genbank_download/reprocessed_2022-06-11/AB267275.gb")
    # refdict = reference
    # seqrecord = next(seqrecordgen)
    seqq, logq, termq = queues

    logq.put(f"{seqrecord.name}: started processing")

    checkissues, refissues, mitfiissues = None, None, None
    if dochecking:
        seqrecord, checkissues = do_checking(seqrecord, (logq, termq))

    if any(checkissues):
        logq.put(f"{seqrecord.name}: processing - issues are present")

    if args.forcereannotate or any(checkissues):
        if refdict:
            seqrecord, refissues = ref_annotate(seqrecord, refdict, tempdir, (logq, termq))
        else:
            seqrecord, mitfiissues = mitfi_annotate(seqrecord, args.table, tempdir, (logq, termq))

    termq.put(len(seqrecord.features))
    seqq.put(seqrecord)
    #Bio.SeqIO.write(seqrecord, "test.gb", "genbank")
    return checkissues, refissues, mitfiissues


def getcliargs(arglist=None):
    # Initialise the parser object and give a description
    parser = argparse.ArgumentParser(description="""
        Add or correct existing tRNA annotations on a mitochondrial genome. Supply a fasta or 
        genbank file on STDIN comprising one or more sequences.
        |n
        Annotations are added either by running MitFi on the input sequence or by comparison 
        against a reference sequence, supplied as a genbank to the -r/--reference argument. If 
        this file contains a single sequence, it will be used as a reference for all inputs; if
        the reference contains the same number of sequences as the input, with the same names, the 
        matching files will be used; otherwise, supply a csv or tsv to -m/--referencemap that 
        details which reference (column one) should be used for each input sequence (column two)
        |n
        If the input is a genbank file, by default the tRNAs will be checked first, and if names 
        do not contain anticodon sequences, the qualifier tags will be checked for anticodon 
        sequences. If all tRNAs are already clearly identifiable, the sequence will not be run 
        through MitFi or compared against the reference. To force re-annotation, use 
        -f/--forcereannotate.
        """, formatter_class=MultilineFormatter)

    # Add individual argument specifications
    parser.add_argument("-r", "--reference", type=str, metavar="PATH",
                        help="path to a genbank file containing sequences with reference tRNAs")
    parser.add_argument("-m", "--referencemap", type=str, metavar="PATH",
                        help="path to a tsv or csv mapping input sequences to reference "
                             "sequences, if >1 reference sequence and sequence names don't match")
    parser.add_argument("-f", "--forcereannotate", action='store_true',
                        help="if input is a genbank file, skip checking existing tRNAs and always "
                             "reannotate using references or MitFi")
    parser.add_argument("-b", "--table", type=int, metavar="N", default=5,
                        help="translation table number")
    parser.add_argument("-t", "--threads", type=int, metavar="N", default=1,
                        help="number of threads to use")

    # Parse the arguments from the function call if specified, otherwise from the command line
    args = parser.parse_args(arglist) if arglist else parser.parse_args()

    # Do some checking of the inputs
    # if args.filepath is not "what you want"
    #     parser.error(f"{args.filpath} is not what I want!")

    # If the arguments are all OK, output them
    return args


# Start the actual script
if __name__ == "__main__":
    interactive = False
    # Get the arguments
    args = getcliargs()
    # interactive, args = True, getcliargs('-t 1'.split(' '))

    # Initialise queue manager and pool
    manager = multiprocessing.Manager() if not interactive else None
    pool = multiprocessing.Pool(args.threads + 3) if not interactive else None

    # Start the accessory threads
    queues, writers = start_accessories(pool, manager)

    # Read input
    seqrecordgen, check = parse_input()
    # seqrecordgen, check = parse_input("/home/thomas/work/iBioGen_postdoc/MMGdatabase/genbank_download/reprocessed_2022-06-11/AB267275.gb")

    # Overwrite checking if genbank but forceannotate has been used
    check = False if args.forcereannotate else check

    # Create temporary directory
    tempdir = "fixtrnas_temp"
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    # Read the references, if supplied
    reference, refalndir = None, None
    if args.reference:
        reference = parse_reference(args.reference, args.referencemap)
    # Otherwise get mitfi
    else:
        get_mitfi(tempdir)

    # Do the work
    issues = pool.map(functools.partial(process_seqrecord, check, args,
                                        reference, tempdir, queues),
                      seqrecordgen)

    # Delete the temporary reference alignment directory
    if tempdir:
        shutil.rmtree(tempdir)

    # Process issues
    #process_issues(issues)

    # Push any issues in any accessory queues to the writers then close the multiprocessors
    for q in queues:
        if q is not None:
            q.put(None)
    pool.close()
    pool.join()

    # Pull any error messages from the writers
    for w in writers:
        if w is not None:
            w.get()

    exit()
