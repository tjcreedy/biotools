#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Thomas J. Creedy
"""

# Imports

import os
import re
import argparse
import textwrap as _textwrap
import MITOcorrect_modules as mcm
import shutil
import multiprocessing
import functools

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

# Argument parser

parser = argparse.ArgumentParser(description="""description:
                                                |n
                                                Text
                                                |n
                                                Text
                                                """,
                                 formatter_class=MultilineFormatter)

parser._optionals.title = "arguments"

parser.add_argument('-t', '--threads', type = int)
parser.add_argument('-n', '--namevariants') # Additional name variants file, optional
parser.add_argument('-s', '--specifications') # Required, specifications file
parser.add_argument('-a', '--alignmentpaths') # Optional?, alignments file
parser.add_argument('-l', '--logfile') # Optional, log file
parser.add_argument('-g', '--genbank', nargs = '+') # Required, one or more genbank files
parser.add_argument('-b', '--translationtable')
parser.add_argument('-f', '--framefree', default = False, action = 'store_true')# Optional - flag to say that search strings do not start in frame.
parser.add_argument('-c', '--alignmenttype') # aa or nt
parser.add_argument('-o', '--outputdirectory')
parser.add_argument('-k', '--keepalignments', default = False, action = 'store_true')
parser.add_argument('-r', '--detailedresults', default = False, action = 'store_true')
parser.add_argument('-p', '--potentialfeatures', default = False, action = 'store_true')
parser.add_argument('-1', '--onefile', type = str) # Output all input gb files in one output file, given as argument
parser.add_argument('-i', '--alignmentweight', type = float, default = 0.5)
parser.add_argument('-j', '--overlapweight', type = float, default = 0.5)

# Main
if __name__ == "__main__":
    
    args = parser.parse_args()
    
    arglist = ("-s MITOcorrect_specs.tsv "
              "-g /home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-04-25_current/BIOD00949.gb "
              "-l testlog.txt "
              "-a aaalignfile.tsv "
              "-o testout/ "
              "-t 2 -b 5 -c aa -r -1 out.gb").split()
    #-g dir/test_multigenbank.gb
    #-g /home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-04-25_current/CCCP00094.gb
    #os.chdir('/home/thomas/seqtesting/MITOcorrect_testing')
    #args = parser.parse_args(arglist)
    
    # Parse the arguments into the main utility variables
    utilityvars = mcm.initialise(args)
    
    # Initialise the queue manager and pool
    manager = multiprocessing.Manager()
    pool =  multiprocessing.Pool(args.threads + 4)
    
    # Start the writers first in their own threads
    seqq, statq, logq, prinq, watch = mcm.start_writers(pool, manager, args)
    
    # Do the work
    seqrecordgen = mcm.get_seqrecords(args.genbank, args.onefile)
    out = pool.map(functools.partial(mcm.process_seqrecord, args, utilityvars, 
                                     seqq, statq, logq, prinq),
                   seqrecordgen)
    
    
    # TODO: Process issues dict
    
    seqq.put('kill')
    if args.detailedresults: statq.put('kill')
    logq.put('kill')
    prinq.put('kill')
    pool.close()
    pool.join()
    
    # Delete temporary alignment directory
    if not args.keepalignments:
        shutil.rmtree(utilityvars[3])
        # delete output
    
    exit()
