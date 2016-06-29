#!/bin/sh

# Personal Identification Pipeline - setup scripts
# https://github.com/TeamErlich/personal-identification-pipeline
#
# Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)
# All Rights Reserved.
# This software is restricted to educational, research, not-for-profit purposes.
# See LICENSE file for full details.

set -e

# This script is setup the prerequisites for the
# personal-identification-pipeline on an Ubuntu-14.04 system.
#
# NOTES:
#  libncurses-dev is needed for samtools tview.
#  libhdf5-dev is needed for python/pip modules.
#  libfreetype6-dev, libpng12-dev - for rebuilding matplotlib (if needed later)

apt-get update
apt-get install -y build-essential git wget zlib1g-dev pkg-config \
                   libncurses-dev libhdf5-dev \
                   libfreetype6-dev libpng12-dev \
                   python python-pip python-dev python-numpy python-scipy \
                   python-matplotlib \
