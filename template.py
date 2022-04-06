#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Description of what the script does"""

# Imports

import argparse
import textwrap as _textwrap

from Bio import SeqIO

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

def getcliargs(arglist=None):
    # Initialise the parser object and give a description
    parser = argparse.ArgumentParser(description="""
        Here describe what the code does. I often find it useful to draft this before I even write
        the code, because it helps me think things through from the end-user perspective (but it 
        may be that a little experience is needed before you can think like this and that's OK!)
        |n
        Separate paragraphs with the symbol above. It doesn't have to be on its own line but I find
        it nicer
        """, formatter_class=MultilineFormatter)

    # Add individual argument specifications
    parser.add_argument("-f", "--flag", action='store_true',
                        help="this argument stores True if used and False otherwise")
    parser.add_argument("-p", "--filepath", type=str, metavar="PATH", required=True,
                        help="this argument records the path to a file, it's required")
    parser.add_argument("-n", "--number", type=int, metavar="N", default=2,
                        help="this argument records an integer, it must be 1, 2, 3 or 4. If not "
                             "supplied, 2 will be used",
                        choices=[1,2,3,4])

    # Parse the arguments from the function call if specified, otherwise from the command line
    args = parser.parse_args(arglist) if arglist else parser.parse_args()

    # Do some checking of the inputs
    if args.filepath is not "what you want"
        parser.error(f"{args.filpath} is not what I want!")

    # If the arguments are all OK, output them
    return args

# Start the actual script
if __name__ == "__main__":

    # Get the arguments
    args = getcliargs() # Try to read from the command line, won't work interactively!
    args = getcliargs('-f -p mypath.txt -n 4'.split(' ')) # Read from a string, good for testing
