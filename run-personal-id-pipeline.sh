#!/bin/sh

version="
run-personal-id-pipeline.sh: part of the Personal Identification Pipeline
https://github.com/TeamErlich/personal-identification-pipeline

Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)
All Rights Reserved.
This software is restricted to educational, research, not-for-profit purposes.
See LICENSE file for full details.

Version 0.3.3
"

set -u

die()
{
    BASE=$(basename "$0")
    echo "$BASE: error: $*" >&2
    if test -n "$outputDir" ; then
        echo "  output directory: $outputDir" >&2
        echo "  logs: $outputDir/run-log" >&2
    fi
    exit 1
}

log()
{
    # Write message to log file with a timestamp
    _timestamp=$(date +"%F-%H%M%S.%N")
    echo "$_timestamp: $*" >> "$outputDir/run-log/messages"

    # if verbose, also print to STDOUT
    test -z "$verbose" && return
    echo "$_timestamp: $*"
}

show_help_and_exit()
{
    BASE=$(basename "$0")
    echo \
"$version

This script runs the MinION Personal Identification pipeline.

It extracts SNPs from a minION sequencing run, compares them against
a know set of 23-and-Me files, and reports possible matches.

Usage: $BASE [OPTIONS] Sample-Name MinION-Reads-DIR

Parameters:
    Sample-Name:         Output files will be created with this name.
    MinION-Reads-DIR:    Directory containing *.FAST5 files from a minION run.

Options:
    -1	    = run matching algorithm on only 1D minION reads (fwd,rev)
              (default is 2D only)
    -a	    = run matching algorithm on all 1D and 2D minION reads
              (default is 2D only)
    -b      = output basenames of 23-and-Me files in probability scores output.
              (defailt: full path).
              see --output-basename in calc-match-prob.py.
    -f FILE = FASTA of human reference genome file with corresponding BWA
              built. Default: search for ./hg19/hg19.fa in same directory as
              this script.
    -t DIR  = directory containing *.txt files in 23-and-Me format.
              The MinION/FAST5 sample will be compared against each file
              in this directory. Either -t or -T must be used.
    -T      = Skip the last part of the pipeline (matching against 23-and-Me
              files). Either -t or -T must be used.
    -c      = create Match plot file plotting p(match) for all candidates
    -e 0.xx = use an error rate of 0.xx  (default: 0.15)
              see -e in calc-match-prob.py.
    -h      = This help screen.
    -o DIR  = Create and write output files to DIR.
              (default: 'output-YYYY-MM-DD-HHMMSS' based on current time)
    -p N    = use N CPUs (default: auto-detect)
    -q      = Use quality scores for error rate.
              (default: fixed 15% error rate).
              see -q in calc-match-prob.py.
    -s SNP  = SNP database (snp138common.txt file).
              Default: use file in current directory (if exists).
    -v      = be verbose
              (add --debug in calc-match-prob.py)

"
    exit 0
}

parse_parameters()
{
    ##
    ## Parse and validate command-line parameteres
    ##
    show_help=
    verbose=
    use_q=
    error_rate=
    runAllP=0
    createPlot=0
    output_basename=
    cpus=$(nproc --ignore 1 2>/dev/null || echo 1)
    readType="2D"
    outputDir=
    original_parameters="$*"
    snp_database=""
    humanRefGen=
    candidatesDir=
    skip_candidate_matching=

    # Parse parameters
    while getopts 1abce:f:ho:p:qs:t:Tv param
    do
        case $param in
            1)   readType="fwd,rev" ;;
            a)   readType="all" ;;
            b)   output_basename=1;;
            c)	 createPlot=1 ;;
            e)   error_rate="$OPTARG"
                 echo "$error_rate" | grep -qE '^0\.[0-9]+$' \
                     || die "invalid error rate '-e $error_rate'"
                 ;;
            f)   humanRefGen="$OPTARG";;
            h)   show_help=1;;
            o)   outputDir="$OPTARG"
                 # being extra strict here, but it's to protect the users
                 # from being stupid...
                 echo "$outputDir" | grep -qE '^[-a-zA-Z0-9_+\/\.=]+$' \
                     || die "output directory name contains forbidden " \
                            "characters ($outputDir)"
                 ;;
            p)   cpus="$OPTARG"
                 echo "$cpus" | grep -qE '^[0-9]+$' \
                     || die "invalid cpus '-p $cpus'"
                 ;;
            q)   use_q="_usingQ";;
            s)   snp_database="$OPTARG";;
            t)   candidatesDir="$OPTARG";;
            T)   skip_candidate_matching=1;;
            v)   verbose=1;;
            ?)   die "unknown/invalid command line option";;
        esac
    done
    shift $(($OPTIND - 1))

    test -n "$show_help" && show_help_and_exit

    # Parameters passed to 'calc-match-prob.py'
    _x=
    test -n "$output_basename" && _x="$_x -b"
    test -n "$verbose" && _x="$_x -v"
    test -n "$use_q" && _x="$_x -q"
    test -n "$error_rate" && _x="$_x -e $error_rate"
    calc_prob_params=$_x

    # Check positional parameters
    test $# -eq 2 || die "expecting 2 parameters. See -h for help."

    sampleName="$1"
    fast5InputDir="$2"

    # being extra strict here, but it's to protect the users
    # from being stupid...
    echo "$sampleName" | grep -qE '^[-a-zA-Z0-9_+\.=]+$' \
        || die "sample-name ($sampleName) contains forbidden characters"
    test -d "$fast5InputDir" \
        || die "minION-reads-DIR '$fast5InputDir' is not a valid directory"

    # Human Genome Reference File + BWA Index:
    if test -z "$humanRefGen" ; then
        # autodetect if non is given on the command line
        _dir=$(dirname "$0")
        if test -e "$_dir/hg19/hg19.fa.bwt" ; then
            # found in the expected location, as generated by
            # "./setup/setup-hg19.sh"
            humanRefGen="$_dir/hg19/hg19.fa"
        else
            die "Missing human-genome-reference FASTA file and BWA index. " \
                "Create it with ./setup/setup-hg19.sh , or specify with -f " \
                "FILE. See -h for help."
        fi
    else
        # Human Genome Reference provided on command line,
        # verify it exists
        test -e "$humanRefGen" \
            || die "Genome FASTA File '$humanRefGen' does not exist"
        test -e "$humanRefGen.bwt" \
            || die "BWA index for FASTA File '$humanRefGen.bwt' not found"
    fi


    # 23-and-Me candidates directory
    test "$candidatesDir" && test "$skip_candidate_matching" \
        && die "-t and -T are mutually exclusive. Use only one of them. " \
               "See -h for help."
    test -z "$candidatesDir" && test -z "$skip_candidate_matching" \
        && die "Please use either '-t DIR' or '-T' for 23-and-Me matching or" \
               "no matching. See -h for help."
    if test "$candidatesDir" ; then
        test -d "$candidatesDir" \
            || die "'$candidatesDir' (-t) is not a valid directory"
    fi


    # dbSNP file
    if test -z "$snp_database" ; then
        # try to auto detect file in script's directory
        _dir=$(dirname "$0")
        if test -e "$_dir/snp138Common.fixed.txt.gz" ; then
            # found in the expected location, as generated by
            # "./setup/setup-snp138common.sh"
            snp_database="$_dir/snp138Common.fixed.txt.gz"
        else
            die "Missing dbSNP-common tabix file. " \
                "Create it with ./setup/setup-snp138common.sh , or specify" \
                "with -s FILE. See -h for help."
        fi
    else
        test -e "$snp_database" \
            || die "snp database file (-s $snp_database) not found"
    fi


    # if outputDir not specified, use a sane default
    test -z "$outputDir" \
        && outputDir="${sampleName}-$(date +%F-%H%M%S)"

    # Create Output Directory
    mkdir "$outputDir" \
        || die "failed to create output directory '$outputDir'"
}

log_run_parametars()
{
    # log the runtime parameters to a file, for later troubleshooting
    _dir="$outputDir/run-log"
    mkdir "$_dir" || die "failed to create run-log directory ($_dir)"

    echo "$sampleName" > "$_dir/sampleName" \
        || die "failed to write sampleName log"
    echo "$fast5InputDir" > "$_dir/minION-reads-DIR" \
        || die "failed to write minION-reads-DIR log"
    echo "$candidatesDir" > "$_dir/23-and-Me-DATA-DIR" \
        || die "failed to write 23-and-Me-DATA-DIR log"
    echo "$humanRefGen" > "$_dir/genome-reference" \
        || die "failed to write genome-reference log"
    echo "$readType" > "$_dir/readType" \
        || die "failed to write readType log"
    echo "$error_rate" > "$_dir/error_rate" \
        || die "failed to write error_rate log"
    echo "$use_q" > "$_dir/using_q" \
        || die "failed to write using_q log"
    echo "$snp_database" > "$_dir/snp-database" \
        || die "failed to write snp-database log"

    date --rfc-3339=seconds > "$_dir/start-timestamp" \
        || die "failed to write start-timestamp log"

    echo \
"#!/bin/sh

cd '$PWD'

$0 $original_parameters
" > "$_dir/re-run.sh" || die "failed to write re-run.sh script"
}

detect_sub_scripts_path()
{
    # Determine the path for the required sub-scripts
    # of this pipeline.
    # This enables running the pipeline from another directory.
    # 1. If this script ($0) was executed with relative path,
    #    the scripts should be in the same directory.
    # 2. If this script was not run with a relative directory, two options:
    #    2.1. The script was with 'sh pipeline.sh' - the scripts are in the
    #         current directory.
    #    2.2. The script was added to the PATH, expect the other scripts to be
    #         in the PATH as well.
    #
    # This function must be called with "$0" as the parameter.
    #
    # This function updates $PATH as needed.

    script_dir=$(dirname "$0")

    if test "x$script_dir" != "x." ; then
        # This script was run with a relative - expand it, add it to the PATH
        script_dir=$(realpath "$script_dir")
    else
        # Script was not run with a relative path, try to find the other
        # scripts in the current directory
        if test -e sam-to-bedseq.py ; then
            # Scripts are in the current directory
            script_dir=$PWD
        fi
    fi

    # Update PATH if needed
    if test -n "$script_dir" ; then
        PATH="$script_dir:$PATH"
        export PATH
    fi

    # Now let's try to find all scripts
    for s in \
        poretools-basenames.py \
        sam-to-bedseq.py \
        sam-discard-dups.py \
        generate-snp-list.py \
        calc-match-probs.py ;
    do
        which "$s" >/dev/null 2>&1 \
            || die "required script '$s' not found in PATH ($PATH)"
    done
}

##
## Program Starts Here
##
detect_sub_scripts_path "$0"

parse_parameters "$@"

log_run_parametars

log "pipeline starting"
log "Sample Name:         $sampleName"
log "Output dir:          $outputDir"
log "minION directory:    $fast5InputDir"
log "23-and-Me directory: $candidatesDir"
log "SNP database:        $snp_database"
log "run logs directory:  $outputDir/run-log/"

##
## Use poretools to convert FAST5 files to FASTQ,
## extract read timing and statistics.
##
log "running: poretools fastq"
poretools fastq --type "$readType" "$fast5InputDir" \
          > "$outputDir/$sampleName.fullpath.fastq" \
    || die "poretools fastq failed on '$fast5InputDir'"

# remove the fullpath from the fastq file,
# set the sequence-id to the filename.
poretools-basenames.py --fastq-id "$outputDir/$sampleName.fullpath.fastq" \
                              > "$outputDir/$sampleName.fastq" \
    || die "poretools-basename.py failed on '$sampleName.fullpath.fastq'"


log "running: poretools times"
poretools times "$fast5InputDir" > "$outputDir/$sampleName.fullpath.times" \
    || die "poretools times failed on '$fast5InputDir'"

# remove the fullpath from the times file
poretools-basenames.py --times "$outputDir/$sampleName.fullpath.times" \
                              > "$outputDir/$sampleName.times" \
    || die "poretools-basename.py failed on '$sampleName.fullpath.times'"


log "running: poretools stats"
poretools stats --type "$readType" "$fast5InputDir" \
          > "$outputDir/$sampleName.stats.txt" \
    || die "poretools stats failed on '$fast5InputDir'"

# If the generated files are empty, bail out with an error
test -s "$outputDir/$sampleName.fastq" \
    && test -s "$outputDir/$sampleName.times" \
    || die "generated FASTQ/TIMES files are empty " \
           "($outputDir/$sampleName.{fastq,times}) - " \
           "perhaps a readType mismatch (all/2d/fwd/rev), " \
           "or a quality error?"

##
## Run BWA.
##
log "running: bwa"
bwa mem -t "$cpus" -x ont2d "$humanRefGen" \
            "$outputDir/$sampleName.fastq" \
            > "$outputDir/$sampleName.dups.sam" \
    || die "bwa mem failed on '$outputDir/$sampleName.fastq'"

##
## Remove duplicates (minION artifacts?)
##
sam-discard-dups.py "$outputDir/$sampleName.dups.sam" \
                    > "$outputDir/$sampleName.sam" \
    || die "sam-discard-dups failed on '$outputDir/$sampleName.dups.sam'"

##
## Create BEDSeq file (BED file + normalized sequences),
## as pre-processing step to SNP-Calling
log "running: sam-to-bedseq"
sam-to-bedseq.py --times "$outputDir/$sampleName.times" \
                 "$outputDir/$sampleName.sam" \
                 > "$outputDir/$sampleName.unsorted.bedseq" \
    || die "sam-to-bedseq.py failed"

test -s "$outputDir/$sampleName.unsorted.bedseq" \
    || die "BEDSeq output is empty '$outputDir/$sampleName.unsorted.bedseq'" \
           " (perhaps no reads mapped, or all mapped to negative strand?)"

##
## Sort BEDSeq by arrival-number
## (sed is used to skip the header line)
( sed -u 1q ; sort -k8n,8 ) < "$outputDir/$sampleName.unsorted.bedseq" \
                            > "$outputDir/$sampleName.bedseq" \
    || die "failed to sort BEDSeq file"


## Perform "SNP Calling" by comparing the minION sequences (from the BEDSeq file)
## against the SNPs data (e.g. snp138common).
##
## Generate a list of SNPs, which will be compared against
## multiple 23-and-Me files (below).
log "running: pipeline step 2"
_v=""
test "$verbose" && _v="--verbose"
generate-snp-list.py $_v "$outputDir/$sampleName.bedseq" "$snp_database" \
                     > "$outputDir/$sampleName.unsorted.snps" \
    || die "generate-snp-list.py (step 2) failed"

test -s "$outputDir/$sampleName.unsorted.snps" \
    || die "data error: the final SNP list is empty " \
           "($outputDir/$sampleName.unsorted.snps)! " \
           "no reads matched to any SNPs?"

log "Sorting SNP list"
# column 13 is arrival-num
( sed -u 1q ; sort -k13n,13 ) < "$outputDir/$sampleName.unsorted.snps" \
                              > "$outputDir/$sampleName.snps" \
    || die "failed to sort SNP list"


##
## Compare the generated SNP list against the 23-and-Me collection.
##
if test "$candidatesDir" ; then
    log "running: pipeline step 3 (calc-probabilities)"

    find "$candidatesDir" \( -type f -o -type l \) -print0 \
        | xargs -0 -I% -n1 -P"$cpus" \
            stdbuf -oL \
             calc-match-probs.py $calc_prob_params "$outputDir/$sampleName.snps" % \
             > "$outputDir/$sampleName.unsorted.matches" \
        || die "calc-match-probs.py failed"

    log "Sorting Match Probabilities"
    # column 1 is 23-and-Me ref-id,
    # column 6 is read arrival-num.
    ( sed -u 1q ; sort -k1,1 -k6n,6 ) < "$outputDir/$sampleName.unsorted.matches" \
                                      > "$outputDir/$sampleName.matches" \
        || die "failed to sort match-probabilities"
else
    log "skipping 23-and-Me matching (-T)"
fi

log "completed."
log "output Directory: $outputDir"
log "SNP list: $outputDir/$sampleName.SNP-list.txt"
test "$candidatesDir" \
    && log "matches:  $outputDir/$sampleName.matches.tsv"
