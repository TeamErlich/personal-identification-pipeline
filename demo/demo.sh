#!/bin/sh

# Map to the genome reference
threads=$(nproc)
bwa mem -V -x ont2d -t "$threads" ../hg19/hg19.fa demo.fq > demo.dups.sam

../sam-discard-dups.py demo.dups.sam > demo.sam

# Create BEDSeq file (.t = temp, unsorted)
../sam-to-bedseq.py --times demo.times demo.sam > demo.bedseq.t

# sort (with header line)
( sed -u 1q ; sort -k8n,8 ) < demo.bedseq.t > demo.bedseq
rm demo.bedseq.t

# Generate candidate SNP list
../generate-snp-list.py demo.bedseq ../snp138Common.fixed.txt.gz > demo.snps.t

# sort (with header line)
( sed -u 1q ; sort -k13n,13 ) < demo.snps.t > demo.snps
rm demo.snps.t

# Match against Yaniv's 23-and-Me file
# (will return empty results, as the simulated SNP will not match)
file=../Yaniv_Erlich_Genome.txt

../calc-match-probs.py demo.snps "$file" > demo.matches
