#!/usr/bin/env Rscript

# Personal Identification Pipeline - Plotting Scripts
# https://github.com/TeamErlich/personal-identification-pipeline
#
# Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)
# All Rights Reserved.
# This software is restricted to educational, research, not-for-profit purposes.
# See LICENSE file for full details.

options(warn=-1)

suppressPackageStartupMessages(library(gplots,quietly=TRUE))
suppressPackageStartupMessages(library(RColorBrewer,quietly=TRUE))
suppressPackageStartupMessages(library(hexbin,quietly=TRUE))
suppressPackageStartupMessages(library(naturalsort,quietly=TRUE))
library(optparse,quietly=TRUE)

parser <- OptionParser(
  usage = paste("%prog FILE.bedseq",
                "This program plots information about a minION run based on",
                "BEDSeq file from 'sam-to-bedseq.py' output.",
                sep="\n"),
  epilog="Output PDF file will be named as input file + '.pdf' extension"
)

## Hack note: far from ideal, but the tryCatch+is.na()
## will print a friendlier error on parsing error
## instead of R's default cryptic message.
arguments=NA
tryCatch(
  { arguments = parse_args(parser, positional_arguments=1);},
  error = function(e) { })
if (all(is.na(arguments))) {
  stop (paste("Failed to parse command-line parameters",
              "(missing filename?).",
              "Use --help for help"))
}
file = arguments$args
if (is.na(file)) {
  stop ("missing file name. Use --help for help")
}


##
## Plotting starts here
##

bedseqfile = file
outfile = paste(bedseqfile,".pdf",sep="")

# Load the first 21 columns, skip the last 4 ones (the very large CIGAR,MD,SEQ,QUAL strings)
col_to_load = c(rep(NA,21),rep("NULL",4))
data = read.table(bedseqfile,comment.char="",sep="\t",header=T,quote="", colClasses = col_to_load)

pdf(outfile)

##
## Front Page
##


# Silly hack:
# use an empty plot with text strings to draw some information in the output PDF
plot(0:5, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(1,4,  sprintf("minION run: %s", bedseqfile), pos=4)
text(1,3,  sprintf("total mapped/usable reads: %s", 
                    format(nrow(data), big.mark=",", scientific=FALSE)),pos=4)
text(1,2,  sprintf("total aligned bases (after mapping): %s",
                   format(sum(data$nt_align), big.mark=",", scientific=FALSE)),
              pos=4)


##
## minION read length - over time
## (not accumulated reads - hence hexbin to hint at individual reads)
hexbinplot(data$orig_seq_len ~ data$arrival_time_rel, aspect=0.5,
            main="read-lengths over time\n(only successfully mapped read)",
            xlab="arrival time (seconds since expermt-start)",
            ylab="minION read length (nt)" )


##
## Accumulated counts - over time
## (2 plots in page)

par(mfrow=c(2,1)) 

##
## accumulated number of reads over time
##
c =  cumsum(rep(1,length(data$arrival_time_rel)))

plot(data$arrival_time_rel, c, pch=20,
     xlab="arrival time (seconds since expermt-start)",
     xaxt='n',
     ylab="accum. num. of minION reads")

# Hack: draw x-axis as hours:minutes:seconds
f = axTicks(1)
k = data.frame(epoch=f, hours=trunc(f/60/60), minutes=trunc((f%%3600)/60), seconds=(f%%3600)%%60)
l = apply(k, 1, function(x){ sprintf("%02d:%02d:%02d\n(%ds)", x[2],x[3],x[4],x[1]) })
axis(1,at=f, labels=l)

##
## Accumulated Genome Mapping counts - over time
## (multiple counts in one chart)
plot(data$arrival_time_rel, cumsum(data$orig_seq_len),type='l', col="black",
     xlab="arrival time (seconds since expermt-start)", xaxt='n',
     ylab="accum. num. nucleotides")
lines(data$arrival_time_rel, cumsum(data$nt_align), col="green")
lines(data$arrival_time_rel, cumsum(data$nt_match), col="red")
lines(data$arrival_time_rel, cumsum(data$nt_mismatch), col="blue")
lines(data$arrival_time_rel, cumsum(data$nt_ins), col="purple", lty=3)
lines(data$arrival_time_rel, cumsum(data$nt_del), col="brown", lty=4)
lines(data$arrival_time_rel, cumsum(data$nt_sclip), col="orange", lty=4)

legend ( 'topleft', 
         c("minION read length",
           "total mapped to ref, of which:",
           "  mapped + matches to ref",
           "  mapped + mismatch from ref (=SNP/error)",
           "insertions",
           "deletions",
           "soft-clips"
         ),
         bty="n",
         lty = c(1,1,1,1,3,4,4), 
         col = c("black","green","red","blue","purple","brown","orange"))

# Hack: draw x-axis as hours:minutes:seconds
f = axTicks(1)
k = data.frame(epoch=f, hours=trunc(f/60/60), minutes=trunc((f%%3600)/60), seconds=(f%%3600)%%60)
l = apply(k, 1, function(x){ sprintf("%02d:%02d:%02d\n(%ds)", x[2],x[3],x[4],x[1]) })
axis(1,at=f, labels=l)


##
## Four distributions (in one page)
##

par(mfrow=c(2,2))

##
## summary of reads/chrom
##
b = table(data$chrom)
# order by natsort(chrom)
b = b[naturalsort(rownames(b))]
barplot(b, main="minION reads per Chrom", ylab="num. reads", xlab="chrom", las=2)

##
## nuc. error rate
##
hist(data$nt_mm_err_rate, main="nuc. mapping error rate\nin minION reads", 
     sub="0.00 =no mismatches, 1.00=all mismatches",
     xlab="1 - ( match-count / aligned-count )",
     ylab="num. reads")

##
## alignment counts
##
hist(data$nt_align,main="alignment count distribution\nin minION reads", 
     sub = "(insertions/deletions/clippping not counted)",
     xlab="num. nucleotides aligned to ref genome",
     ylab="num. reads")

##
## mismatches counts
##
hist(data$nt_mismatch,main="mismatches count distribution\nin minION reads", 
     sub = "(insertions/deletions/clippping not counted)",
     xlab="num. nucleotides aligned but differ from ref genome (=SNPs or errors)",
     ylab="num. reads")


dummy = dev.off()
