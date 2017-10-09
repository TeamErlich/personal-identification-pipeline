#!/bin/sh

# Personal Identification Pipeline - setup scripts
# https://github.com/TeamErlich/personal-identification-pipeline
#
# Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)
# All Rights Reserved.
# This program is licensed under GPL version 3.
# See LICENSE file for full details.

#
# Downloand and pre-process the snp138Common track
# (from the UCSC Genome Browser).
#
# NOTE:
# Tabix squashes multiple tabs into one,
# so "A\t\tB" is considered 2 fields instead of 3
# (with middle field being empty string).
# Similarly, lines ending with "\t$" will have one less field.
#
# This script injects a dummy character to avoid such cases
# (e.g. 'rs7450001')


die()
{
    BASE=$(basename "$0")
    echo "$BASE: error: $*" >&2
    exit 1
}

URL=http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp138Common.txt.gz

raw=$(basename "$URL" .gz)

# Step 1: Download the original file from UCSC
if ! test -e "$raw.gz" ; then
  wget -O "$raw.gz" "$URL" \
    || die "failed to download snp file '$URL'"
fi

# Step 2: decompress it
if ! test -e "$raw" ; then
  gunzip -k "$raw.gz" \
    || die "failed to decomperess '$raw.gz'"
fi

# check whether it's the known version
# this is the md5sum of:
#    snp138Common.txt.gz   27-Apr-2014   14:42    621M
# from:
#    http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/
c="$raw"
EXP="892152efbfa1c712934dec9a4255e77b  -"
SUM=$(md5sum < "$c") \
    || die "failed to calculate md5sum on '$c'"
test "x$EXP" = "x$SUM" \
    || {
          echo \
"WARNING: downloaded snp138Common.txt file is not the same
          version as we expect. The known/tested version is:
            http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/

            file: snp138Common.txt.gz
            timestamp:  27-Apr-2014   14:42
            size: 621M

            after decompression:
              size:   3937024218
              md5sum: 892152efbfa1c712934dec9a4255e77b

          The pipeline will still work, but might produce slightly
          different results.
" >&2
       }

# Prevent sequential tabs (i.e. empty fields) from confusing tabix
sed -e 's/		/	**	/g' -e 's/	$/	0/' \
    "$c" > "$c.t" \
    || die "sed failed on '$c'"

# Filter out some funky snps:
#   anything except SNPs ( $12 == "single" )
#   anything except bi-allelic SNPs ( $22 == 2)
#   where refNCBI != refUCSC ( $8 == $9 )
#   and any exceptions found by UCSC ( $19 == "**" )
awk '($12=="single") &&
     ($8==$9) &&
     ($22==2) &&
     ($19=="**") { print > "snps138.good" ; next }
                 { print > "snps138.rejected" }' "$c.t" \
    || die "failed to filter SNPs from '$c.t'"

c="${c%.txt}.fixed.txt"
mv "snps138.good" "$c" \
    || die "failed to rename snp138.good to '$c.fixed.txt'"

# And re-compress it using 'bgzip'
# (which is required for tabix)
bgzip -c "$c" > "$c.gz" \
    || die "failed to re-compress '$c' using bgzip"

# Build tabix index
tabix -0 -s 2 -b 3 -e 4 "$c.gz" \
    || die "failed to build tabix index on '$c.gz'"

test -e "$c.gz.tbi" \
    || die "index file '$c.gz.tbi' not found (after running tabix)"

echo \
"snp138Common tabix pre-processing complete.

To query file manually, try:

   tabix '$c.gz' chr14:30000000-31000000

"
