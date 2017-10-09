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
suppressPackageStartupMessages(library(naturalsort,quietly=TRUE))
library(optparse,quietly=TRUE)

parser <- OptionParser(
  usage = paste("%prog FILE.snps",
                "This program plots information about a called-SNPS from a minION run",
                "based on SNPS file from 'generate-snp-list.py' output.",
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

snpsfile = file
outfile = paste(file,".pdf",sep="")
data = read.table(snpsfile,comment.char="",sep="\t",quote="", header=T)

pdf(outfile)

# Takes a number, returns a string with the nubmer and thousand separators
# tsep(1000) => "1,000"
tsep=function(x) { format(x, big.mark=",", scientific=FALSE) }

##
## Front Page
##

# Silly hack:
# use an empty plot with text strings to draw some information in the output PDF
plot(0:5, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(1,4,  sprintf("minION run: %s", snpsfile), pos=4)
text(1,3,  sprintf("total SNPs reads: %s", tsep(nrow(data))),pos=4)

##
## 
##
valid_snps = data [ data$af > 0, ]
t = c(
       sprintf("valid SNPs: %s  of which:", tsep( nrow(valid_snps) ) ),
       sprintf("  match REF: %s", tsep(   sum(valid_snps$snp_ref_base==valid_snps$base) ) ),
       sprintf("  match ALT: %s", tsep(   sum(valid_snps$snp_alt_base==valid_snps$base) ) ),
       sprintf("invalid SNPs: %s", tsep(  sum(data$af == 0) ) )
     )
hist(data$af, main="Distribution of Allele-Frequency for called SNPS\n0 = invalid SNP (doesn't match ref or alt)",
      xlab="allele-frequency (from dbSNP-138)",
      ylab="number of SNPS")
legend("topleft",paste(t,sep="\n"),bty='n',border='none')



##
## Valid SNP accumulation over time
##

v = table(valid_snps$arv_time_rel)
x = as.integer(rownames(v))
y = cumsum(as.vector(v))
plot(x,y,type='l',xaxt='n', main="number of valid SNPs over time",
     xlab="time (since experiment-start)",
     ylab="accum. num. of valid SNPs") 

# Hack: draw x-axis as hours:minutes:seconds
f = axTicks(1)
k = data.frame(epoch=f, hours=trunc(f/60/60), minutes=trunc((f%%3600)/60), seconds=(f%%3600)%%60)
l = apply(k, 1, function(x){ sprintf("%02d:%02d:%02d\n(%ds)", x[2],x[3],x[4],x[1]) })
axis(1,at=f, labels=l)


t = valid_snps [ valid_snps$snp_ref_base == valid_snps$base, ]
v = table(t$arv_time_rel)
x = as.integer(rownames(v))
y = cumsum(as.vector(v))
lines(x,y,col="red")

t = valid_snps [ valid_snps$snp_alt_base == valid_snps$base, ]
v = table(t$arv_time_rel)
x = as.integer(rownames(v))
y = cumsum(as.vector(v))
lines(x,y,col="blue")

legend ( 'topleft', 
         c("all valid SNPs,of which:",
           "  match REF",
           "  match ALT"
         ),
         bty="n", lty = 1,
         col = c("black","red","blue"))

##
## Avg. SNPS per read
##
snps_per_read = table(valid_snps$seq_id)
hist(snps_per_read,breaks = 20,
     main="Average valid SNPs per read",
     xlab="number of valid SNPs",
     ylab="number of mapped minION reads")




##
## Invalid SNPs break-down
##
invalid_snps = data [ data$af == 0, ]

# NOTE:
# Ideally the diagonal should be zero
# i.e. A=>A, because if the sequenced nucleotide matched the reference,
# why did we matk it as invalid with AF=0 in 'generate-snp-list.py' ??
# BUT,
# dbSNP contains SNPs which have allele-frequencies of 1.0/0.0, so these are
# not actually invalid allelse, just counted as such.
# Example: rs7237188
#
m = table ( invalid_snps$snp_ref_base, invalid_snps$base )
heatmap.2(m,Rowv = NA, Colv=NA, dendrogram = 'none',trace='none', col=brewer.pal(9,"PuRd"),
          main = "Invalid SNPs cross-tabulation\n(SNPs not matching either ALT/REF)",
          xlab = "base called from minION",
          ylab = "dbSNP ref base")

