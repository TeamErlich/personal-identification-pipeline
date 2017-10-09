#!/usr/bin/env Rscript

# Personal Identification Pipeline - Plotting Scripts
# https://github.com/TeamErlich/personal-identification-pipeline
#
# Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)
# All Rights Reserved.
# This program is licensed under GPL version 3.
# See LICENSE file for full details.


##
## TODO:
## Change the path of the '.matches' file you'd like to plot
## (The output of 'calc-prob-matches.py' script).
file="NO-SUCH-FILE"

data = read.table(file, header=T,sep="\t",comment.char = "",quote = "")

sample_id = basename(file)

# Save resutls to a PDF file
pdf(paste(sample_id,".pdf",sep=""))

# Create an empty plot with pre-set axes
# x => 0 to 3600 seconds (=60 minutes)
# y = 0 to 1 (probability)
sample_id = basename(file)

plot(0,0,pch=0,col="white", xlim=c(0,3600),ylim=c(0,1),
      main=paste("minION Sample",sample_id,"vs DNA.Land users"),
      xlab="minION read arrival time (minutes since experiment-start)",
      ylab="Match probability", xaxt="n")
# Add labels every 10 minutes
axis(1, at=c(0,10,20,30,40,50,60)*60,
        labels=c(0,10,20,30,40,50,60))

## Get list of sample_ids  in the file
items = levels(data$ref_id)

## Iterate each sample and plot it
for (i in items) {
  if (i=="yaniv_erlich") {
    col = "red"
  } else {
    col = "black"
  }

  # Extract sample's row (i.e. the sample vs Yaniv's DNA.Land-23-and-Me-file)
  d = data[ data$ref_id==i, ]

  # Keep the first 60 minutes
  d = d[ d$arv_time_rel <= 60*60, ]

  # Find the max probability for each second
  # (multiple reads,SNPs may arrive in a given second)
  a = tapply(d$posterior, d$arv_time_rel, FUN=max)

  x = as.integer(rownames(a))
  y = a

  # Draw the data for this sample
  lines(x,y, type='l',col=col)
}

## Add a legend (last - to be drawn on top of the lines)
legend("bottomright",
       c("Yaniv Erlich's DNA.Land user","DNA.Land users"),
       col=c("red","black"),
       lty=1
)

