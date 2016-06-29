#!/bin/sh

# Personal Identification Pipeline - setup scripts
# https://github.com/TeamErlich/personal-identification-pipeline
#
# Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)
# All Rights Reserved.
# This software is restricted to educational, research, not-for-profit purposes.
# See LICENSE file for full details.


set -e

# This script is setup downloads,builds and installs bwa
# for the personal-identification-pipeline.
# It assumes the current (non-root) user can run sudo
# (will ask for password if needed).
#

# BWA website:   https://github.com/lh3/bwa
# BWA releases:  https://github.com/lh3/bwa/releases

cd /tmp

VER=0.7.15
URL=https://github.com/lh3/bwa/releases/download/v$VER/bwa-$VER.tar.bz2
FILE=$(basename "$URL")

wget -O "$FILE" "$URL"
tar -xf "$FILE"
cd bwa-$VER
make
sudo -p "Enter (sudo) password copy 'bwa' ($VER) to /usr/local/bin': "\
   cp bwa /usr/local/bin/
