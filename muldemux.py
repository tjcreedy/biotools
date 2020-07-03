#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""""

# Imports

import sys
import argparse
import csv

import os
import re
import subprocess
import textwrap as _textwrap

from Bio import SeqIO

# Function definitions

def required_multiple(multiple):
    class RequiredMultiple(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not len(values)%multiple == 0:
                msg='argument "{f}" requires a multiple of {multiple} values'.format(
                    f=self.dest,multiple=multiple)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredMultiple

def str_is_int(s):
	try: 
		int(s)
		return True
	except ValueError:
		return False


def std_well(well):
	well_search = re.search('^' + wellregex + '$', well)
	if well_search and str_is_int(well_search.group(2)):
		if 0<int(well_search.group(2))<13 :
			return (well_search.group(1) + '%02d' % int(well_search.group(2)))


def find_regex_in_str(string, regex, regex_borders = ['^', '_', '-', '$']):
	#string, regex = [name, directregex]
	regexlist = {b1+regex+b2 for b1 in regex_borders if b1 != '$' for b2 in regex_borders if b2 != '^' and ({b1, b2} != {''})}
	regexhits = {regex:re.findall(regex, string) for regex in regexlist if len(re.findall(regex, string)) > 0}
	matchlist = {m for ml in regexhits.values() for m in ml}
	
	if(len(matchlist) > 1):
		return('ambiguous')
	elif(len(matchlist) == 0):
		return('none')
	else:
		regexhit = max(regexhits.keys(), key = len)
		
		[match] = matchlist
		return(re.sub(regexhit, '', string), ''.join(match), int(match[1]) if len(match)>1 and int(match[1]) else None)

def parse_input_files(inputfiles, conversion):
	#inputfiles = args.input
	
	data = {}
	conversion_needed = False
	ext = None
	warned = False
	n_wells_detected = 0
	
	for path in inputfiles:
		#path = inputfiles[0]
		
		# Extract file name parts
		file = os.path.basename(path)
		name = file
		if('.' in file):
			name, ext = file.split('.', 1)
			ext = '.'+ext
		elif(not warned):
			sys.stderr.write("Warning: file " + file + " has no detectable file extension: uncompressed fastq will be assumed but this may cause errors. This warning will not be repeated\n")
			warned = "True"
		
		# Parse salient details from name
		well, direction = [find_regex_in_str(name, wellregex), find_regex_in_str(name, directregex)]
		
		# Check salient details
		ident = None
		if(direction == "ambiguous"):
			sys.exit("Error: found multiple possible R1 or R2 in file" + path + "\n")
		elif(direction == "none"):
			sys.exit("Error: could not find R1 or R2 in file " + path + "\n")
		else:
			ident, suffix, d = direction
		
		if(well in ["ambiguous", "none"] or conversion is not None):
			conversion_needed = True
		else:
			n_wells_detected += 1
			ident = std_well(well[1])
		
		# Store into dict
		if ident in data.keys() :
			if(None in data[ident]['files']):
				data[ident]['files'][d-1] = path
			else:
				sys.exit("Error: found more than two files for " + ident + "\n")
		else:
			data[ident] = {'files':[None, None], 'name':ident}
			data[ident]['files'][d-1] = path
	
	# Check all input names have two files
	
	for ident, subdict in data.items():
		if(None in subdict['files']):
			sys.exit("Error: only one file found for " + subdict['name'] + "\n")
	
	return(data, conversion_needed, ext, n_wells_detected)

def parse_conversion(path, datadict):
	#path, data = args.conversion, data
	
	missing_names = []
	
	with open(path, 'r') as tab_file:
		for n, row in enumerate(csv.reader(tab_file, delimiter = '\t')):
			n = str(n+1)
			if(len(row) != 2):
				sys.exit("Error: line " + n + " of " + path + " does not contain two elements\n")
			
			name = row[0]
			if(name not in datadict.keys()):
				namematches = {dname for dname in datadict.keys() if re.match(name, dname)}
				
				e = " matches for " + name + ", " + path + " line " + n + " in the supplied filenames"
				if(len(namematches) > 1):
					sys.exit("Error: found multiple" + e + ":\n\t" + "\n\t".join(namematches) + "\n")
				elif(len(namematches) == 0):
					missing_names.append(name)
					continue
				
				[name] = namematches
			
			well = std_well(row[1])
			if(well):
				datadict[well] = datadict.pop(name)
			else:
				sys.exit("Error: could not understand well in " + path + " line " + n + "\n")
	
	for name in datadict.keys():
		if(std_well(name) is None):
			sys.exit("Error: no conversion was found for " + name + "\n")
	
	return(datadict, missing_names)


def parse_demuxtables(demuxtables, datadict):
	#demuxtables, datadict = [args.demuxtable, data]
	
	missing_demux_names = []
	no_well_demux = []
	namelist = []
	
	for table in demuxtables:
		#table = demuxtables[0]
		with open(table, 'r') as tab_file:
		#tab_file =  open(table, 'r')
		#close(tab_file)
			for n, row in enumerate(csv.reader(tab_file, delimiter = '\t')):
				#n, row = list(enumerate(csv.reader(tab_file, delimiter = '\t')))[0]
				n = str(n+1)
				if(len(row) != 3):
					sys.exit("Error: line " + n + " does not contain 3 elements\n")
				
				outname = row[0]
				
				if(outname in namelist):
					sys.exit("Error: \" " + outname + " appears multiple times in demultiplexing table(s)\n")
				else:
					namelist.append(outname)
				
				well = find_regex_in_str(outname, wellregex)
				e =  " in " + outname + ", " + table + " line " + n 
				if(well == "ambiguous"):
					sys.exit("Error: found multiple possible wells" + e + "\n")
				elif(well == "none"):
					sys.exit("Error: could not identify any wells" + e + "\n")
				else:
					well = std_well(well[1])
				
				if(well in datadict.keys()):
					if('demux' not in datadict[well].keys()):
						datadict[well]['demux'] = {'table':[[],[]], 'indices':[set(),set()]}
					
					if(row[1:] in datadict[well]['demux']['table'][1]):
						sys.exit("Error: the index pair " + row[1] + " and " + row[2] + " appear multiple times for well " + datadict[well]['name'] + " in demultiplexing tables(s)\n")
					
					for i, r in enumerate([outname, row[1:]]):
						datadict[well]['demux']['table'][i].append(r)
					for i, x in enumerate(row[1:]):
						datadict[well]['demux']['indices'][i].add(x)
					
				else:
					missing_demux_names.append(outname)
	
	for well in list(datadict.keys()):
		if 'demux' not in datadict[well].keys():
			no_well_demux.append(well + ": " + datadict[well]['name'])
			del datadict[well]
	
	
	return(datadict, missing_demux_names, no_well_demux)

# Class definitions

class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width, initial_indent=indent, subsequent_indent=indent) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text

# Global variables

wellregex = '([A-H])([1-9]|0[1-9]|1[0-2])'
directregex = '(R)([12])'


# Arguments

parser = argparse.ArgumentParser(description= """description:
                                                 |n
                                                 Apply cutadapt paired-end trimming to a set of paired files, applying indices from a supplied table to demultiplex the files. The paired files should be supplied to --input, and pairs will be identified based on file names containing '_R1' and '_R2' denoting forward and reverse reads respectively. Files may be fasta or fastq, and may be compressed in any way that cutadapt accepts. In order to accurately generate output file names, input names should not contain '.' symbol apart from denoting file types, e.g. 'filename.fa', 'filename.fastq.gz' are permitted, 'file.name.fasta' will cause unexpected and possibly clashing output file names, as the '.name.' section will be lost.
									|n
									Indices are read from one or more demultiplexing tables supplied to --demuxtable arguments, which must be a three-column tab-separated text file each line of which contains an output name, the forward index, and the reverse index.  Indices will be matched to input read files based on well numbers. Thus, input files and output names must contain well numbers in the format 'A1' or 'A01' and must be preceded by - or _ and/or be followed by - or _ and/or be at the end of the name. Where input files do not contain well numbers, a conversion table must be supplied to --conversion in a two-column tab-separated format, where each line gives a suitably unique name for the input file followed by a well number. 
                                                 |n
                                                 The basic cutadapt command used is 'cutadapt -cores N -g O1=INDEX -G O2=INDEX [ -g O2=INDEX2 -G O2=INDEX2 ] -o IN_R1 -p IN_R2 OUT_R1 OUT_R2'. The number of cores for cutadapt to run on will be passed directly from the --threads argument. Further arguments can be supplied to cutadapt using --arguments '-v W -x Y'. These will be passed to cutadapt unchecked.
                                                 |n
                                                 Output files will be written to a directory named in --output, which will be created if necessary. A statistics file will be written to --statistics if specified.""",  
                                  formatter_class=MultilineFormatter)

parser._optionals.title = "arguments"

parser.add_argument("-i", "--input", help = "paths to two or more fastx files to demultiplex, required. These may be compressed.", type = str, metavar = "IN_R1 IN_R2", required = True, nargs ='+', action=required_multiple(2))
parser.add_argument("-t", "--threads", help = "the number of cores/threads to run cutadapt on", type = int, default = 1)
parser.add_argument("-d", "--demuxtable", help = "path to a demultiplexing table, required", type = str, metavar = "X", action='append', required=True)
parser.add_argument("-c", "--conversion", help = "path to a conversion table", type = str, metavar = "X", required=False)
parser.add_argument("-a", "--arguments", help = "further arguments to pass to cutadapt, in a single quoted list", type = str, metavar = "'-v W -x Y'")
parser.add_argument("-^", "--anchored", help = "anchor the indices to the start of the reads", action = 'store_true')
parser.add_argument("-o", "--output", help = "path to a directory to store output files", type = str, metavar = "X", required = False, default = "./")
parser.add_argument("-s", "--statistics", help = "path to a file to store demultiplexing statistics", type = str, metavar = "X", required = False)
parser.add_argument("-v", "--verbose", help = "print details of of found/missing names when parsing file list, demux tables and/or conversion table", action = 'store_true')
parser.add_argument("-k", "--keeperrors", help = "don't delete files for incorrect index combinations", action = 'store_true')
parser.add_argument("-p", "--printcutadapt", help = "print cutadapt stdout/stderr to the terminal", action = 'store_true')

# Main

if __name__ == "__main__":
	
	#args = parser.parse_args(['-i', '/home/thomas/Documents/NHM_postdoc/iBioGen/MetagenWorkshop/AMM/resources/metabarcoding/0_rawsequences/Lib1_R1.fastq', '/home/thomas/Documents/NHM_postdoc/iBioGen/MetagenWorkshop/AMM/resources/metabarcoding/0_rawsequences/Lib1_R2.fastq', '/home/thomas/Documents/NHM_postdoc/iBioGen/MetagenWorkshop/AMM/resources/metabarcoding/0_rawsequences/Lib2_R1.fastq', '/home/thomas/Documents/NHM_postdoc/iBioGen/MetagenWorkshop/AMM/resources/metabarcoding/0_rawsequences/Lib2_R2.fastq', '/home/thomas/Documents/NHM_postdoc/iBioGen/MetagenWorkshop/AMM/resources/metabarcoding/0_rawsequences/Lib3_R1.fastq', '/home/thomas/Documents/NHM_postdoc/iBioGen/MetagenWorkshop/AMM/resources/metabarcoding/0_rawsequences/Lib3_R2.fastq', '/home/thomas/Documents/NHM_postdoc/iBioGen/MetagenWorkshop/AMM/resources/metabarcoding/0_rawsequences/Lib4_R1.fastq', '/home/thomas/Documents/NHM_postdoc/iBioGen/MetagenWorkshop/AMM/resources/metabarcoding/0_rawsequences/Lib4_R2.fastq', '-d', 'TEMP/demux.txt', '-c', 'TEMP/convert.txt', '-o', 'TEMP/out/', '-s', 'TEMP/stats.txt', '-w', '-j', '20'])
	
	args = parser.parse_args()
	
	# Sort input files
	
	sys.stdout.write("Checking input files...")
	
	data, conversion_needed, ext, n_wells_detected = parse_input_files(args.input, args.conversion)
	ffmt = "fasta" if re.match('a$|a\.', ext) else "fastq"
	
	sys.stdout.write("%s file pairs parsed" % (len(data)))
	if(args.verbose):
		sys.stdout.write(":\n\t%s\n" % ("\n\t".join(data.keys())))
	else:
		sys.stdout.write("\n")
	
	
	# Parse conversion table and rename files dict
	
	if(conversion_needed):
		if(args.conversion):
			sys.stdout.write("Parsing conversion table...")
			data, convert_missing_names = parse_conversion(args.conversion, data)
			sys.stdout.write("done\n")
			
			if(args.verbose and len(convert_missing_names) > 0):
				sys.stdout.write("The following conversion table entries did not match to any supplied file pairs:\n\t" + "\n\t".join(convert_missing_names) + "\n")
		else:
			nerr = len(data)*2 - n_wells_detected
			sys.exit("Error: well numbers cannot be ascertained from file names for " + str(nerr) + ", a conversion table is required\n")
	else:
		if(args.conversion):
			sys.stderr.write("Warning: supplied conversion table not needed\n")
	
	
	# Parse demultiplexing tables
	sys.stdout.write("Parsing conversion table(s)...")
	
	data, demux_missing_names, well_missing_demux = parse_demuxtables(args.demuxtable, data)
	
	sys.stdout.write("done\n")
	
	if(args.verbose):
		if(len(demux_missing_names) > 0):
			sys.stdout.write("The following demultiplexing table(s) entries did not match to any supplied file pairs:\n\t" + "\n\t".join(demux_missing_names) + "\n")
		if(len(well_missing_demux) > 0):
			sys.stdout.write("The following file pairs had no entries in the demultiplexing table(s) and their associated files will not be processed:\n\t" + "\n\t".join(well_missing_demux) + "\n")
	
		
	# Set up output directory
	
	if not os.path.exists(args.output):
		os.makedirs(args.output)
	
	
	
	# Set up writing to statistics file
	if(args.statistics):
		stats = open(args.statistics, 'w')
		stats_write = csv.writer(stats, delimiter = '\t', quotechar = '', quoting = csv.QUOTE_NONE, escapechar = '')
	
	# Loop through input file pairs
	
	for well, specs in data.items():
		sys.stdout.write("Running cutadapt on file pair %s: %s..." % (well, specs['name']))
		
		#well, specs = list(data.items())[0]
		
		# Generate and run cutadapt command
		
		cutargs = ['cutadapt -cores ' + args.threads]
		if(args.arguments): cutargs.extend(re.split(' +', args.arguments))
		
		cutargs.extend([o for a, d in zip(['-o', '-p'],['R1', 'R2']) for o in [a, os.path.join(args.output, well + '_{name1}-{name2}_' + d + '.' + ffmt)]])
		eq = '=^' if args.anchored else '='
		cutargs.extend([o for a, d in zip(['-g', '-G'], specs['demux']['indices']) for i in d for o in [a, i + eq + i]])
		cutargs.extend(specs['files'])
		
		cutcmd = subprocess.run(cutargs, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
		
		cutstderr = cutcmd.stderr.decode("utf-8")
		cutstdout = cutcmd.stdout.decode("utf-8")
		
		# Write outputs
		if(args.printcutadapt):
			sys.stdout.write("\n")
			sys.stdout.write("Cutadapt command: %s \n" % (' '.join(cutargs)))
			sys.stdout.write(cutstdout)
			sys.stderr.write(cutstderr)
		
		sys.stdout.write("processing cutadapt outputs %s..." % ("and generating statistics " if args.statistics else ""))
		
		# Find all indices
		forindex, revindex = specs['demux']['indices']
		
		if(not args.arguments or '--discard-untrimmed' not in args.arguments):
			forindex.add('unknown')
			revindex.add('unknown')
		
		# Work through all combinations of indices for this file pair
		
		for f in forindex:
			for r in revindex:
				#f, r = [next(iter(forindex)), next(iter(revindex))]
				# Determine if this combination is a target
				demuxtab = specs['demux']['table']
				name = demuxtab[0][demuxtab[1].index([f,r])] if [f, r] in demuxtab[1] else None
				
				# Find the files for this combination
				
				files = [os.path.join(args.output, well + '_' + f + '-' + r + '_' + d + '.' + ffmt) for d in ['R1', 'R2']]
				
				nseqs = 0
				if os.path.exists(files[0]):
					nseqs = len(list(SeqIO.parse(files[0], ffmt)))
				
					# Rename the files if they are a target, otherwise delete them
					if(name):
						for file, d in zip(files, ['R1', 'R2']):
							os.rename(file, os.path.join(args.output, name + '_' + d + '.' + ffmt))
					else:
						if(not args.keeperrors):
							for file in files: os.remove(file)
				else:
					if(name):
						sys.stderr.write("No output file could be found for %s, either there was an error or no matches for the corresponding index pair were found by cutadapt\n" % (name))
				
				# Output the statistics if reporting is on
				
				if(args.statistics):
					stats_write.writerow([well, specs['name'], f, r, nseqs, name if name else ''])
		sys.stdout.write("done\n")
	
	
	if args.statistics : stats.close()
	
	
	exit()
