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
parser.add_argument('-1', '--onefile', default = False, action = 'store_true') # Output all input gb files in one output file
parser.add_argument('-i', '--alignmentweight', default = 0.5)
parser.add_argument('-j', '--overlapweight', default = 0.5)

# Main
if __name__ == "__main__":
    
    args = parser.parse_args()
    
    arglist = re.sub('dir', '/home/thomas/MITOcorrect_testing',
            """-s dir/testspecs.tsv
               -g /home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-04-25_current/BIOD00411.gb
               -l dir/testlog.txt
               -a dir/test_ntalignfile.tsv
               -o dir/testout/ 
               -t 2 -b 5 -c nt -k -r -p -1""").split()
    #-g dir/test_multigenbank.gb
    #os.chdir('/home/thomas/MITOcorrect_testing')
    #args = parser.parse_args(arglist)
    
    # Parse the arguments into the main utility variables
    namevars, annotypes, specs, temp, log, stath, statw = mcm.initialise(args)
    
    # Set up the first file and output set
    gblist, gbname, seqn = mcm.get_file_details(args.genbank.pop(0))
    outrecords = []
    seqn = 0
    
    # Loop through genbank files and entries within each.
    while True:
        
        # Retreive the next record to work on
        seqrecord = next(gblist, None)
        
        # If generator is empty, write out the results file unless user 
        # requests all outputs to be in the same output, then 
        # move onto the next file if available, or exit.
        if seqrecord is None:
            if len(outrecords) > 0:
                if not args.onefile:
                    outrecords = mcm.write_genbank_file(outrecords,
                                                        args.outputdirectory,
                                                        gbname)
            if len(args.genbank) > 0:
                gblist, path, gbname, ext = mcm.get_file_details(
                        args.genbank.pop(0))
                continue
            else:
                if args.onefile:
                    outrecords = mcm.write_genbank_file(outrecords,
                                                        args.outputdirectory,
                                                        gbname)
                break
        seqn += 1
        
        
        # Extract the necessary items from the seqrecord and clean
        present, cleanfeats, ofeats, issues, plog = mcm.prepare_seqrecord(
                                                    seqn, seqrecord, gbname,
                                                    namevars, annotypes, specs)
        log.write(plog)
        
        # Process the present cleanfeatures in parallel
        
        with multiprocessing.Pool(processes=args.threads) as pool:
             poolout = pool.map(functools.partial(mcm.correct_feature,
                                                  cleanfeats, specs, gbname, 
                                                  seqrecord, args, temp), 
                                present, 1)
        
        flatten = lambda l: [item for sublist in l for item in sublist]
        outfeats, logl, statsl = map(flatten, zip(*poolout))
        
        log.write(''.join(logl))
        
        if args.detailedresults:
            for l in statsl:
                statw.writerow(l)
        
        # Replace all features with the new ones and add on the others
        if len(outfeats) > 0:
            seqrecord.features = outfeats + ofeats
        
        # Append completed seqrecord to outputs records
        outrecords.append(seqrecord)
        
        # End of loop on gb
    
    # TODO: Process issues dict
    
    
    # Delete temporary alignment directory
    if not args.keepalignments:
        shutil.rmtree(temp)
        # delete output
    
    # Close logfile
    log.close()
    
    # Close stats file
    if args.detailedresults:
        stath.close()
    
    exit()
