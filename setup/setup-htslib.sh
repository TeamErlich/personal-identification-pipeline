#!/bin/sh

# Personal Identification Pipeline - setup scripts
# https://github.com/TeamErlich/personal-identification-pipeline
#
# Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)
# All Rights Reserved.
# This program is licensed under GPL version 3.
# See LICENSE file for full details.

set -e

# This script is setup downloads,builds and installs samtools,bwa
# for the personal-identification-pipeline.
# It assumes the current (non-root) user can run sudo
# (will ask for password if needed).
#

# Samtools instructions in  http://www.htslib.org/download/

cd /tmp

VER=1.3.1
URL=https://github.com/samtools/htslib/releases/download/$VER/htslib-$VER.tar.bz2
FILE=$(basename "$URL")

wget -O "$FILE" "$URL"
tar -xf "$FILE"
cd htslib-$VER
./configure
make

sudo -p "Enter password for 'sudo make install' (htslib $VER): "\
    make install
