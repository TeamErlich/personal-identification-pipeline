#!/bin/sh

version="
Personal-ID-pipeline - run on DNA.Land samples.

https://github.com/TeamErlich/personal-identification-pipeline

Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)
All Rights Reserved.

This software is restricted to educational, research, not-for-profit purposes.
See LICENSE file for full details.
"

set -u

# Default values
cpus=$(nproc --ignore 2)

# ignore reads arriving after this time (since experiment-start-time)
cutoff_seconds=3600

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

Usage: $BASE [OPTIONS] INPUT-SNP-FILE 23-and-ME-DIR

Options:
    --help      = This help screen.
    --cpus N    = use N CPUs (default: auto-detect)

calc-match-prob.py options:
(passed as is, see 'calc-match-prob.py --help')

    --fixed-error-rate E
    --nucq-error-rate
    --read-mm-error-rate

    --prob-matrix-new
    --prob-matrix-old

    --initial-prior INITIAL_PRIOR
    --output-basename

    --cutoff-time SECONDS
    --cutoff-low-probability P

For further information about input/output file formats,
see 'calc-match-prob.py.py --help' .
"
    exit
}


parse_config()
{
    cpus=$(nproc --ignore 2)
    debug=
    xargs_params=
    params=
    while test $# != 0 ; do
        case "$1" in
            --help)
                show_help_and_exit
                ;;

            --cpus)
                shift
                cpus="$1"
                ;;
            --cpus=*)
                cpus=${1#*=}
                ;;

            --debug)
                debug=1
                xargs_params="--verbose"
                params="$params --debug"
                ;;

            --nucq-error-rate|--read-mm-error-rate|--prob-matrix-new|--prob-matrix-old|--output-basename|--debug-genotypes)
                params="$params $1"
                ;;

            --fixed-error-rate=*|--cutoff-time=*|--cutoff-low-probability=*)
                params="$params $1"
                ;;

            --fixed-error-rate)
                shift
                params="$params --fixed-error-rate '$1'"
                ;;

            --cutoff-time)
                shift
                params="$params --cutoff-time '$1'"
                ;;

            --cutoff-low-probability)
                shift
                params="$params --cutoff-low-probability '$1'"
                ;;

            --no-header)
                die "parallel calc-prob does not take --no-header."
                ;;

            --only-header)
                die "parallel calc-prob does not take --only-header."
                ;;

            --)
                shift
                break;;

            -*)
                die "Unknown option '$1'. See --help for help."
                ;;

            *) break;;
        esac
        shift
    done

    test $# -lt 2 \
        && die "missing input file and 23-and-Me directory. See --help for help."
    test $# -gt 2 \
        && die "extra operand ($3). See --help for help."

    # first parameter - filename
    input="$1"
    test -e "$input" \
        || die "input file '$input' not found"


    # Check directory with 23-and-Me files
    candidate_dir="$2"
    test -d "$candidate_dir" \
        || die "can't find inpu directory with 23-and-Me files ($candidate_dir)"
}

parse_config "$@"


output=${input%.snps}.matches
test -e "$output" && die "output file '$output' already exists, aborting."


## Find the location of the python script
script=
_dir=$(dirname "$0")
script="$_dir/calc-match-probs.py"
test -e "$script" \
    || die "can't find required script ($script)"

python "$script" $params --only-header > "$output.unsorted" \
	|| die "failed to print header line"

## NOTE: some files are invalid, and calc-match-probs.py will exit with code 1.
## don't stop.
nice find "$candidate_dir" \( -type f -o -type l \) -print0 \
    | sort -zV \
    | xargs $xargs_params -0 -I% -n1 -P"$cpus" \
          stdbuf -oL python "$script" $params --no-header "$input" % \
	    >> "$output.unsorted"

# Sort with header line:
# 1st column is sample-id,
# 5th column is read-arrival-time-since-experiment-start
( sed -u 1q ; sort -k1V,1 -k5n,5 -S 100G ) < "$output.unsorted" \
    > "$output.sorted.t" \
    || die "failed to sort file '$output.unsorted.t'"

# rename to final "*.matches" file
mv "$output.sorted.t" "$output" \
    || die "failed to rename '$output.sorted.t' to '$output'"

if test "$debug" ; then
    echo "input file = $input"
    echo "candir     = $candidate_dir"
    echo "cpus       = $cpus"
    echo "params     = $params"
    echo "output     = $output"
fi
