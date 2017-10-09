#!/bin/sh

# Personal Identification Pipeline - setup scripts
# https://github.com/TeamErlich/personal-identification-pipeline
#
# Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)
# All Rights Reserved.
# This program is licensed under GPL version 3.
# See LICENSE file for full details.

# Downloand and pre-process the hg19 genome
# (from the UCSC Genome Browser).

die()
{
    BASE=$(basename "$0")
    echo "$BASE: error: $*" >&2
    exit 1
}

URL=http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

test -e hg19 \
    && die "a file/directory named hg19 already exist. aborting."

test -d hg19 \
    && test -e "hg19/hg19.fa" \
    && die "file hg19/hg19.fa already exists. aborting."


dir=$(mktemp -d hg19.XXXXXX) \
    || die "failed to create temp dir"

cd "$dir" \
    || die "failed to cd into '$dir'"

wget -O "chromFa.tar.gz" "$URL" \
    || die "failed to download '$URL'"

tar -xzvf chromFa.tar.gz \
    || die "failed to extract chromFa.tar.gz"

find . -type f \( -name "chr[0-9].fa" \
                  -o -name "chr[0-9][0-9].fa" \
                  -o -name "chr[XYM].fa" \) \
    | sort -V \
    | xargs cat > hg19.fa \
    || die "failed to creaete hg19.fa"

# Sanity checks: should be 25 sequences (1-22,X,Y,M)
count=$(grep -c "^>chr" < hg19.fa)
test "x$count" = "x25" \
    || die "internal error, file '$dir/hg19.fa' does not contain 25 sequences"


# Build samtools fasta index
samtools faidx hg19.fa \
    || die "samtools faidx failed on hg19.fa"

# Build BWA genome reference index
bwa index hg19.fa \
    || die "bwa index failed on hg19.fa"

# Remove source files
rm chr*.fa \
    || die "failed to clean-up chr*.fa files"

# rename the directory into the final name
cd .. || die "failed to cd .."
mv "$dir" hg19 \
    || die "failed to rename directory '$dir' to 'hg19'"
