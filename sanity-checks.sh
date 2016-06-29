#!/bin/sh

version="
sanity-checks.sh: part of the Personal Identification Pipeline
https://github.com/TeamErlich/personal-identification-pipeline

Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)
All Rights Reserved.
This software is restricted to educational, research, not-for-profit purposes.
See LICENSE file for full details.
"

set -u
set -e

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

Usage: $BASE [PREFIX]

Assumes the following files exist:
  [PREFIX].fastq    - (output of 'poretools fastq')
  [PREFIX].times    - (output of 'poretools times')
  [PREFIX].sam-dups - (output of 'bwa', before discarding duplicates)
  [PREFIX].sam      - (output of 'sam-discard-dups.py')
  [PREFIX].bedseq   - (output of 'sam-to-bedseq.py')
  [PREFIX].snps     - (output of 'generate-snp-list.py')
"
    exit 0
}


test $# -ge 1 || die "missing input prefix. See -h for help"
test "x$1" = "x-h" && show_help_and_exit

prefix="$1"
for ext in fastq times sam sam-dups bedseq snps ;
do
    test -e "$prefix.$ext" \
        || die "file '$prefix.$ext' not found (based on prefix '$prefix')"
done


# Output format string: item, value, comment
fmt="%-20s %10s   %-30s\n"


## number of reads in FASTQ file
l=$(cat "$prefix.fastq" | wc -l)
printf "$fmt" "fastq_num_reads" "$((l/4))"

## number of lines in TIMES file (-1 for header line)
l=$(cat "$prefix.times" | wc -l)
printf "$fmt" "times_num_lines" "$((l-1))" \
        "(number of data lines in $prefix.times)"

## number of reads in SAM file (with dups and non-mappers)
l=$(grep -vc '^@' "$prefix.sam-dups")
printf "$fmt" "sam_dups_num_reads" "$l" \
    "(including non-mappers and multi-mappers)"

## number of non-duplicate reads in SAM file (with non-mappers)
l=$(grep -vc '^@' "$prefix.sam")
printf "$fmt" "sam_num_reads" "$l" "(including non-mappers)"

## number of non-duplicate reads in SAM file (with non-mappers)
l=$(grep -v '^@' "$prefix.sam" | awk '$3 != "*"' | wc -l)
printf "$fmt" "sam_num_mapped_reads" "$l" \
        "(number of mapped reads in $prefix.sam)"

## number of lines in the BEDSeq file (-1 for header line)
l=$(cat "$prefix.bedseq" | wc -l)
printf "$fmt" "bedseq_num_lines" "$((l-1))" \
        "(number of data lines in $prefix.bedseq)"

## Experiment Start time
## NOTE: This assumes local timezone is the same for reporting and experiment
est=$(cut -f4 "$prefix.times" | head -n2 | tail -n1)
printf "$fmt" "exp_start_unix" "$est" "(experiment start, unix time)"
est2=$(date -d "@$est" "+%F")
printf "$fmt" "exp_start_date" "$est2" "(experiment start date)"
est3=$(date -d "@$est" "+%H:%M:%S")
printf "$fmt" "exp_start_time" "$est3" "(experiment start time)"

## Experiment End time (time of last read)
## NOTE: This assumes local timezone is the same for reporting and experiment
eet=$(cut -f7 "$prefix.times" | sed 1d | sort -n | tail -n1)
printf "$fmt" "exp_end_unix" "$eet" "(experiment end (last read), unix time)"
eet2=$(date -d "@$eet" "+%F")
printf "$fmt" "exp_end_date" "$eet2" "(experiment end date)"
eet3=$(date -d "@$eet" "+%H:%M:%S")
printf "$fmt" "exp_end_time" "$eet3" "(experiment end time)"

## Experiment Duration
dur=$((eet-est))
printf "$fmt" "exp_duration_sec" "$dur" "(experiment duration, seconds)"

dur_h=$((dur/3600))
dur_m=$(((dur%3600)/60))
dur_s=$(((dur%3600)%60))
printf "$fmt" "exp_duration_human" "${dur_h}h${dur_m}m${dur_s}s" \
    "(experiment duration)"

## Number of SNPs (-1 for header line)
l=$(cat "$prefix.snps" | wc -l)
printf "$fmt" "snps_num_lines" "$((l-1))" \
        "(number of data lines in $prefix.snps)"

mr=$(awk 'NR>1 && $4==$9 { sum=sum+1} END { print sum}' "$prefix.snps")
printf "$fmt" "snps_match_ref" "$mr"

ma=$(awk 'NR>1 && $5==$9 { sum=sum+1} END { print sum}' "$prefix.snps")
printf "$fmt" "snps_match_alt" "$ma"

af0=$(awk 'NR>1 && $11==0 {sum=sum+1} END { print sum}' "$prefix.snps")
printf "$fmt" "snps_af_zero" "$af0" \
        "(snps not matching REF/ALT marked by af=zero)"

bad1=$(awk 'NR>1 && $11==0 && ($4==$9 || $5==$9) {sum=sum+1}
           END { print sum}' "$prefix.snps")
printf "$fmt" "snps_af0_dbsnp" "$bad1" \
        "(snps that have AF=0 in dbSNP, and will mess things up)"
