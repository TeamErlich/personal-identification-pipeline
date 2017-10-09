personal-identificaion-pipeline - example
=========================================

This directory contains an example dataset to be used with the
personal-identification-pipeline .

The file `demo.fq` contains variations (added SNPs, mismatches, insertions,
deletions) of the following sequence:

    >hg19_dna range=chr2:77870548-77870664 5'pad=0 3'pad=0 strand=+
    TTTATGGCTAGGGTTGCAGGAGTGCATGGTGGGAATGTAAACTGCCAGGT
    ATCGCTCAATTTCCCTTTCTTTGCATGGGGAACCTCTCTAGGCTCCCAGC
    CAATTACAGCGAAGTTG

The file 'demo.times' contains simulated timing information of a minION
run (simulating as if the sequences in `demo.fq` arrived from a minION
sequencer, and were processed using `poretools fastq` and `poretools times`).

The file `demo.sh` runs the individual pipeline steps.
It is equivalent to the full-blown `run-personal-id-pipeline.sh`, but lacks
any informational messages or error checking (easier to read, but not
robust enough for automated production-quality pipeline usage).


License
-------

Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)

All Rights Reserved.
This program is licensed under GPL version 3.
See LICENSE file for full details.
