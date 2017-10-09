#!/bin/sh

# Personal Identification Pipeline - setup scripts
# https://github.com/TeamErlich/personal-identification-pipeline
#
# Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)
# All Rights Reserved.
# This program is licensed under GPL version 3.
# See LICENSE file for full details.

set -e

# This script is setup the prerequisites for the
# personal-identification-pipeline on a CentOS 7.X machine.
#
# NOTES:
#  libncurses-dev is needed for samtools tview.
#  libhdf5-dev is needed for python/pip modules.
#  freetype-devel, libpng-devel - for rebuilding matplotlib

yum groupinstall -y "Development Tools"
yum install -y ncurses-devel git wget zlib-devel python-devel curl \
               libpng-devel freetype-devel

# If python pip is not install, install it manually
# (seems like CentOS does not provide a prebuilt package)
if ! type pip 2>/dev/null 1>/dev/null ; then
    curl "https://bootstrap.pypa.io/get-pip.py" -o "get-pip.py"
    python get-pip.py
fi
