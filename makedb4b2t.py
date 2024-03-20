#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Description of what the script does"""

# Imports
import sys
import tarfile
import sqlite3
import gzip
import argparse
import textwrap as _textwrap

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
        Construct or update one or both of the tables in a SQLite reference database for 
        blast2taxonomy.py. The database should comprise two tables: gb_taxid that maps GenBank 
        accession numbers to taxids, and taxid_lineage that records lineages for a given taxid. 
        Supply the path to the SQLite database to -d/--database. If this doesn't exist, it will be 
        created. If either of the tables already exist, they will be deleted and re-created. The 
        database should have both tables for blast2taxonomy.py to function effectively, and this 
        script can import both in one run, or with separate runs if only one needs updating.
        |n
        To create/update the gb_taxid table, supply a copy of nucl_gb.accession2taxid.gz to 
        -a/--acc2taxid. This can be retrieved from 
        https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/.
        |n
        To create/update the taxid_lineage table, supply a copy of taxdump.tar.gz to -t/--taxdump. 
        This can be retrieved from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/.
        """, formatter_class=MultilineFormatter)

    # Add individual argument specifications
    parser.add_argument("-d", "--database", type=str, metavar="PATH", required=True,
                        help="path to a database")
    parser.add_argument("-a", "--acc2taxid", type=str, metavar="PATH", required=False,
                        help="path to a ncbi_gb.accession2taxid.gz file")
    parser.add_argument("-t", "--taxdump", type=str, metavar="PATH", required=False,
                        help="path to a taxdump.tar.gz file")
    
    # Parse the arguments from the function call if specified, otherwise from the command line
    args = parser.parse_args(arglist) if arglist else parser.parse_args()

    # Do some checking of the inputs
    if not args.acc2taxid and not args.taxdump:
        parser.error("nothing to do! Supply one of -a/--acc2taxid or -t/--taxdump")
    
    # If the arguments are all OK, output them
    return args

# Start the actual script
if __name__ == "__main__":

    # Get the arguments
    args = getcliargs()
    #args = getcliargs("--database /home/thomc/scratch/ncbidb/test.db --taxdump /home/thomc/scratch/ncbidb/taxdump.tar.gz".split(' '))

    # Connect to database
    connection = sqlite3.connect(args.database)
    cursor = connection.cursor()
    tables = cursor.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall()
    tables = [f[0] for f in tables]

    if args.acc2taxid:
        # Check if the table exists, and remove if it does
        if "gb_taxid" in tables:
            cursor.execute("DROP TABLE gb_taxid")
            
        # Create table
        tablesql = ("create table gb_taxid ("
                    "accession TEXT,"
                    "version TEXT,"
                    "taxid INTEGER)")
        cursor.execute(tablesql)

        # Load data into database
        with gzip.open(args.acc2taxid, 'rt') as fh:
            i = 0
            for line in fh:
                values = line.rstrip().split("\t")
                if(values[0] == "accession"):
                    continue
                i += 1
                values = values[:2]
                insertsql = (f"INSERT INTO gb_taxid VALUES ({', '.join('?' for _ in values)})")
                cursor.execute(insertsql, values)
                if i < 10000 or i%1000 == 0:
                    sys.stderr.write(f"\rRead {i} lines from {args.acc2taxid} into gb_taxid")
        sys.stderr.write("\nCreating indices\n")
        cursor.execute("CREATE INDEX index_accession ON gb_taxid(accession)")
        cursor.execute("CREATE INDEX index_version ON gb_taxid(version)")
        connection.commit()

    
    if args.taxdump:
        # Check if the table exists, and remove if it does
        if "taxid_lineage" in tables:
            cursor.execute("DROP TABLE taxid_lineage")
        
        # Create table
        tablesql = ("create table taxid_lineage ("
                    "taxid INTEGER NOT NULL,"
                    "parent INTEGER,"
                    "rank TEXT,"
                    "scientificname TEXT NOT NULL)")
        cursor.execute(tablesql)

        # Load tar
        tar = tarfile.open(args.taxdump, "r:gz")
        # Read names into dict
        namesdict = {}
        with tar.extractfile("names.dmp") as fh:
            i = 0
            for line in fh:
                values = line.decode().rstrip("\t|\n").split("\t|\t")
                namesdict[values[0]] = values[1]
        
        # Load data into database
        with tar.extractfile("nodes.dmp") as fh:
        #fh = tar.extractfile("nodes.dmp")
            i = 0
            for line in fh:
                # line = fh.readline()
                linedata = line.decode().rstrip("\t|\n").split("\t|\t")
                i += 1
                values = linedata[:3] + [namesdict[linedata[0]]]
                if linedata[0] == '1':
                    values[1] = None
                insertsql = (f"INSERT INTO taxid_lineage VALUES ({', '.join('?' for _ in values)})")
                cursor.execute(insertsql, values)
                sys.stderr.write(f"\rRead {i} lines from {args.taxdump} into taxid_lineage")
        tar.close()
        sys.stderr.write("\nCreating indices\n")
        cursor.execute("CREATE INDEX index_taxid ON taxid_lineage(taxid)")
        cursor.execute("CREATE INDEX index_parent ON taxid_lineage(parent)")
        connection.commit()
    
    sys.stderr.write("Optimising database\n")
    cursor.execute("VACUUM;")
    connection.close()
    sys.stderr.write("\nDone\n")