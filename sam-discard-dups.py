#!/usr/bin/env python

"""
sam-discard-dups.py: part of the Personal Identification Pipeline
https://github.com/TeamErlich/personal-identification-pipeline

Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)
All Rights Reserved.
This software is restricted to educational, research, not-for-profit purposes.
See LICENSE file for full details.

Version 0.3
"""

import sys,argparse, errno
from os.path import basename

version_info=\
"""
sam-discard-dups  - version 0.1
"""

def parse_command_line():
    # Define parameters
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="SAM - Discard duplicates/multiple-mappers",
        version=version_info,
        epilog="""
This script reads a SAM file, and discards ALL reads which map
multiple times (does not keep even one of the mapping locations).

TODO: Explain hard-clipped reads in minION runs

Example:
    # Remove all duplicated reads dirnames:

    $ %(prog)s input.sam > output.no-dups.sam

        """)


    # Positional parameter
    parser.add_argument('filename', metavar='FILE', help='file to process');
    args = parser.parse_args()

    return args


def get_dup_read_ids(filename):
    """
    Reads a SAM file, returns a set() of duplicated reads
    (reads which are listed more than once)
    """
    try:
        seen_ids = set()
        dup_ids = set()

        sam=file(filename,'r')

        for linenum,line in enumerate(sam):
            err = "input error in '%s' line %d: " % (filename, linenum+1)
            line = line.strip()

            if line[:1]=='@':
                continue

            flds = line.split('\t')
            read_id = flds[0]

            if read_id in seen_ids:
                dup_ids.add(read_id)

            seen_ids.add(read_id)

        return dup_ids

    except IOError as e:
        if e.errno == errno.EPIPE:
            sys.exit(0) # exit silently

        # TODO: this is a cop-out, hard to tell what's the exact error
        #       and give informative,useful error message to the user.
        sys.exit("I/O error: %s" % (str(e)))


def filter_sam_dups(filename,ids_to_discard):
    """
    Reads a SAM file, returns a set() of duplicated reads
    (reads which are listed more than once)
    """
    try:
        sam=file(filename,'r')

        for linenum,line in enumerate(sam):
            err = "input error in '%s' line %d: " % (filename, linenum+1)
            line = line.strip()

            if line[:1]=='@':
                print line
                continue

            flds = line.split('\t')
            read_id = flds[0]

            if not (read_id in ids_to_discard):
                print line

    except IOError as e:
        if e.errno == errno.EPIPE:
            sys.exit(0) # exit silently

        # TODO: this is a cop-out, hard to tell what's the exact error
        #       and give informative,useful error message to the user.
        sys.exit("I/O error: %s" % (str(e)))



if __name__ == "__main__":
    args = parse_command_line()

    dups = get_dup_read_ids(args.filename)

    filter_sam_dups(args.filename, dups)
