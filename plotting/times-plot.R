#!/usr/bin/env Rscript

# Personal Identification Pipeline - Plotting Scripts
# https://github.com/TeamErlich/personal-identification-pipeline
#
# Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)
# All Rights Reserved.
# This program is licensed under GPL version 3.
# See LICENSE file for full details.

options(warn=-1)

suppressPackageStartupMessages(library(gplots,quietly=TRUE))
suppressPackageStartupMessages(library(RColorBrewer,quietly=TRUE))
suppressPackageStartupMessages(library(hexbin,quietly=TRUE))
library(optparse,quietly=TRUE)

parser <- OptionParser(
  usage = paste("%prog FILE.times",
                "This program plots information about a minION run based on",
                "timing file from 'poretools times' output.",
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

timesfile = file
t = read.table(timesfile,comment.char="",sep="\t",header=T,quote="" )

outfile = paste(timesfile,".pdf",sep="")

##
## Create output PDF file
##
pdf(outfile)

##
## Front Page
##
par(mfrow=c(2,1))

# Silly hack:
# use an empty plot with text strings to draw some information in the output PDF
plot(0:5, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(1,4,  sprintf("minION run: %s", timesfile), pos=4)
text(1,3,  sprintf("total reads: %s",
                   format(nrow(t), big.mark=",", scientific=FALSE)),
                   pos=4)
text(1,2,  sprintf("total base (before mapping): %s",
                   format(sum(t$read_length), big.mark=",", scientific=FALSE)),
                  pos=4)
text(1,1,  sprintf("mean read length: %0.2f bp", mean(t$read_length)), pos=4)
text(1,0,  sprintf("median read length: %f bp", median(t$read_length)), pos=4)

# Distribution plot
hist(t$read_length, main="minION read lengths distribution", xlab="read length (nt)",
     ylab="num. reads")

##
## minION Speed
##

hexbinplot(t$read_length ~ t$duration,
          main="minION sequencing speed",
          xlab="duration (seconds from read-start to read-end)",
          ylab="read length (nt)"
          )

##
## read-lengths over time
##

hexbinplot( t$read_length ~ (t$unix_timestamp_end-t$exp_starttime), aspect=0.5,
            main="minION read length over time",
            xlab="arrival time (seconds since experiment start)",
            ylab="read length (nt)")


##
## Channel Activity: num reads
##

# Summarizes num-reads-per-channel
d = table(t$channel)
# conver to vector of 512 elements
v = rep(0,512)
v[as.integer(rownames(d))] = as.integer(d)
# Reshape as a matrix of 16x32
m = matrix(v, nrow=16)

# create column labels hinting at channel number (e.g. "1-16" for first column)
lbl = sprintf("%d-%d", seq(0,31)*16+1,seq(1,32)*16)

heatmap.2(m,Rowv = NA, Colv=NA, dendogram='none', labRow=seq(0,15),labCol=lbl,trace='none',
          main="minION flowcell channel activity\n (num reads/channel)", col=brewer.pal(9,"Blues"),
          key.title="", key.xlab="num. reads per channel", key.ylab="channel count")

##
## Channel Activity: num bases
##

# Summarizes num-bases-per-channel
d = tapply(t$read_length, t$channel, FUN=sum)
# conver to vector of 512 elements
v = rep(0,512)
v[as.integer(rownames(d))] = as.integer(d)
# Reshape as a matrix of 16x32
m = matrix(v, nrow=16)

heatmap.2(m,Rowv = NA, Colv=NA, dendogram='none', labRow=seq(0,15),labCol=lbl,trace='none',
          main="minION flowcell channel activity\n (num bases/channel)", col=brewer.pal(9,"PuRd"),
          key.title="", key.xlab="num. bases per channel", key.ylab="channel count")


dummy = dev.off()
