#!/bin/sh

# Personal Identification Pipeline - setup scripts
# https://github.com/TeamErlich/personal-identification-pipeline
#
# Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)
# All Rights Reserved.
# This software is restricted to educational, research, not-for-profit purposes.
# See LICENSE file for full details.

set -e

# This script downloads Yaniv Erlich's public Genotype file.
# This file can be used as the reference database of samples to compare
# against.


URL=http://files.teamerlich.org/genomes/Yaniv_Erlich_Genome.txt.gz
FILE=$(basename "$URL")

wget -O "$FILE" "$URL"
gunzip "$FILE"
