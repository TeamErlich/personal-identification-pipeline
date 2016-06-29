#!/bin/sh

version="
run-parallel-calc-prob.sh: part of the Personal Identification Pipeline
https://github.com/TeamErlich/personal-identification-pipeline

Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)
All Rights Reserved.
This software is restricted to educational, research, not-for-profit purposes.
See LICENSE file for full details.

Version 0.3
"

set -u

# minION Personal Identification Pipeline, Step 3.5
# Erlich Lab (http://teamerlich.org)

die()
{
    BASE=$(basename "$0")
    echo "$BASE: error: $*" >&2
    exit 1
}

show_help_and_exit()
{
    BASE=$(basename "$0")
    echo \
"$version

This wrapper scripts runs calc-match-prob.py.py on multiple 23-and-Me files
in parallel.

Usage: $BASE [OPTIONS] INPUT-SNP-FILE  23andMe-DATA-DIR

ALL files in [23-and-Me-DATA-DIR] will be compared against
[INPUT-SNP-FILE]

Options:
    -h      = This help screen.
    -v      = be verbose
              (add --debug in calc-match-prob.py.py)
    -p N    = use N CPUs (default: auto-detect)
    -b      = output basenames of files
              (add --output-basename in calc-match-prob.py.py)
    -e 0.xx = use an error rate of 0.xx (otherwise defaults to fixed 15% error rate)
    -q      = Use quality scores for error rate (otherwise defaults to fixed 15% error rate)

Example:

  \$ $BASE -b -p 10 INPUT.txt /data/ > results.txt

For further information about input/output file formats,
see 'calc-match-prob.py.py --help' .
"
    exit
}


##
## Script start, parse command-line parameteres
##
show_help=
verbose=
output_basename=
use_q=
error_rate=
cpus=$(nproc --ignore 1 2>/dev/null || echo 1)

# Parse parameters
while getopts hvbe:qp: param
do
    case $param in
    h)   show_help=1;;
    v)   verbose=1;;
    b)   output_basename=1;;
    e)   error_rate="$OPTARG";;
    q)   use_q=1;;
    p)   cpus="$OPTARG";;
    ?)   die "unknown/invalid command line option";;
    esac
done
shift $(($OPTIND - 1))

test -n "$show_help" && show_help_and_exit

# Ensure is one non-option parameters (the input file)
test $# -eq 0 && die "missing input SNP file and data directory. See -h for help."
test $# -eq 1 && die "missing 23-and-Me data directory. See -h for help."
test $# -gt 2 && die "too many parameters ($*). See -h for help."

# first parameter - filename
infile="$1"
test -z "$infile" && die "empty input filename. See -h for help."
test -e "$infile" \
    || die "input file '$infile' not found"

# second parameter - data directory with 23-and-Me files
datadir="$2"
test -z "$datadir" && die "empty data directory name. See -h for help."
test -d "$datadir" \
    || die "'$datadir' is not a valid directory"


# extra params
extra=
test -n "$output_basename" && extra="--output_basename"
test -n "$verbose" && extra="$extra --debug"
test -n "$use_q" && extra="$extra --use_q"
test -n "$error_rate" && extra="$extra --error_rate $error_rate"

find "$datadir" \( -type f -o -type l \) -print0 \
    | xargs -0 -I% -n1 -P"$cpus" \
          stdbuf -oL calc-match-probs.py $extra "$infile" %
