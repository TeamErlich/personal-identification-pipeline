#!/usr/bin/env python

version_info=\
"""
calc-match-prob.py: part of the Personal Identification Pipeline
https://github.com/TeamErlich/personal-identification-pipeline

Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)
All Rights Reserved.
This program is licensed under GPL version 3.
See LICENSE file for full details.

Version 0.3
"""

"""
MinION Personal ID Original Pipeline Step 3 calcProb([file/path
name of snplist for sequence data] [file/path name of target candidate
data]

- Run on snplist & file of target candidate
- returns bayesian-derived probability value of the sample belonging to that candidate

"""

from scipy.stats import binom
from collections import defaultdict
from os.path import basename
import numpy as np
import sys, os, argparse

class ProcessingException(Exception):
    """Base class for exceptions in this module."""
    pass


def print_header(debug_genotypes):
    """Print the heaer line.

    debug_genotypes: if True, prints additional columns
                     (see output format for details).
    """
    # Print header
    header = [
        "ref_id",             #1
        "read_id",            #2
        "snp_id",             #3
        "arv_time_abs",       #4
        "arv_time_rel",       #5
        "arv_num",            #6

        "count_snp_num",            #7
        "count_snp_left",           #8
        "count_snp_processed",      #9
        "count_snp_invalid_allele", #10
        "count_snp_not_in_ref",     #11

        "posterior",            #12
        "fam_posterior",        #13

        "prior",                #14
        "fam_prior",            #15

        "p_obs",                #16
        "fam_p_obs",            #17
        "p_obs_given_match",    #18
        "p_obs_given_mismatch", #19
        "p_first_deg",          #20

        "homo_hit",             #21
        "homo_miss",            #22
        "hetero_hit",           #23
        "hetero_miss",          #24
        "hetero_common",        #25
        "hetero_rare",          #26
        "sex_hits",             #27

        "err_rate",             #28
        ]

    if debug_genotypes:
        header.extend( [
            "dbsnp_ref_base",  #29
            "dnsnp_alt_base",  #30
            "dnsnp_ref_af",    #31
            "dnsnp_alt_af",    #32
            "sample_snp_base", #33
            "sample_snp_af",   #34
            "ref_23andMe_gt",   #35
            "ref_23andMe_gt_desc",  #36
            "snp_match_type",  #37
            ])

    print '\t'.join(map(str,header))

class HeaderOnlyAction(argparse.Action):
    """An argparse-compatible class which prints
    the header line then exits immedaitely.
    (this prevents argparse from complaining about missing 'required' parameters
    which aren't required when only printing the header line).
    """
    def __init__(self, option_strings, dest, **kwargs):
        super(HeaderOnlyAction, self).__init__(option_strings, '', nargs=0, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        print_header(False)
        sys.exit(0)



class OutputFormatHelpAction(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super(OutputFormatHelpAction, self).__init__(option_strings, '', nargs=0, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        text=version_info + """
Matched Probabilities Output Format
===================================

This program reads a list of called SNPs (from minION sample),
and compares then against multiple 23-and-Me reference files.

For each 23-and-Me file, it outputs accumulated probabilities
of a match between the sample and the reference 23-and-Me file
(accumulated over the input SNPs from the sample).

Example
If the input sample file has 10 lines (representing 10 SNPs),
The output will contain 10 lines, representing the probability of a
match against the 23-and-Me file using only the first SNP, then first 2 SNPs,
then the first 3 SNPs, etc.

The output will have 18 tab-separated fields:

1.  ref_id           - the 23-and-Me reference file
2.  read_id          - the minION read ID that contained this called SNPs.
3.  snp_id           - the SNP id (e.g. rs12345) from the input file.

4.  arv_time_abs     - Unix timestamp of the arrival time of this minION read.
5.  arv_time_rel     - minION read arrival time, in seconds-since-experiment-start.
6.  arv_num          - sequencial number of minION read arrival order (1=first read)

7.  count_snp_num    - sequential number of SNP processed (1=first in input
                       file). There will be gaps if SNPs are skipped.
8. count_snp_left   - number of SNPs left to process
9. count_snp_processed - number of SNPs processed so far (not skipped)
10. count_snp_invalid_allele - number of skipped SNPs so far due to invalid
                              allele (i.e. differ from dbSNP REF/ALT base)
11. count_snp_not_in_ref - number of skipped SNPs so far due to missing
                           rsIDs in current 23-and-Me file.

12. posterior        - The updated match probability after this SNP is included
                       in the model.
13. fam_posterior    - The updated familial match probability after this SNP is
                       included in the model.

14. prior             - internal parameters from the last round of the model.
15. fam_prior         - internal parameters from the last round of the model.
16. p_obs             - internal parameters from the last round of the model.
17. fam_p_obs         - internal parameters from the last round of the model.
18. p_obs_given_match - internal parameters from the last round of the model.
19. p_obs_given_mismatch - internal parameters from the last round of the model.
20. p_first_deg       - internal parameters from the last round of the model.

21. hom_hits         - number of homoygous hits (i.e. a match between the input
                       sample SNP, when the corresponding 23-and-Me genotype is
                       homozygous).

22. hom_miss         - number of homoygous misses (i.e. 23-and-Me genotype is
                       homozygous, and the input SNP doesn't match it).

23. het_hits         - number of heterozygous hits.

24. het_miss         - number of heterozygous misses.

25. het_common       - TODO

26. het_rare         - TODO

27. sex_hits         - NOT USED, always zero.

28. error_rate       - The error_rate value used for this SNP

if --debug-genotypes is used, the following fields are added:
29. dbsnp_ref_base     -
30. dnsnp_alt_base     -
31. dnsnp_ref_af       -
32. dnsnp_alt_af       -
33. sample_snp_base    -
34. sample_snp_af      - equal to dbsnp_ref_af, dbsnp_alt_af, or other if
                         there's an a error/bug (because invalid alleles should
                         not be counted)
35. ref_23andMe_gt     - The genotype value from the 23-and-Me file
36. ref_23andMe_gt_desc - "hom_ref", "het", "hom_alt", "other" (based on previous field)
37. snp_match_type     - "hom_hit", "hom_miss", "het_hit", "het_miss"
"""




        print (text)
        sys.exit(0)


class ErrorRateMode():
    """poorman's enum"""
    fixed         = 1  # Fixed value
    nucq          = 2  # based on qual-score of sequenced nucleotide
    read_err_rate = 3  # based on match/mismatch error rate of the read

class ProbabilityMatrix():
    """poorman's enum"""
    old           = 1
    new           = 2

def parse_command_line():
    # Define parameters
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="calcProb - minION Personal-Identification pipeline, step 3",
        version=version_info,
        epilog="""

This script reads a called-SNP-list
        (from a minION sample, generated by generate-snp-list.py),
and compares it to given set of 23-and-Me files.

    $ %(prog)s SAMPLE.txt [23andMe-1.txt] [23andMe-2.txt] [...]

Input file should contain 10 fields (see generate-snp-list.py --help).

Matches are written to STDOUT.

See --output-format-help for details about output format.

NOTE:
genotypes/alleles are compared by SNP-ID (e.g. rs12345).
chromosome/positions are ignored.
  """)


    # Option parameters
    parser.add_argument("--debug", help="print debug information to STDERR",
                        action="store_true")

    parser.set_defaults(error_mode=ErrorRateMode.fixed)
    parser.add_argument("--fixed-error-rate", type=float,metavar="E",
                        help="use fixed error rate (default = 0.15)")

    parser.add_argument("--nucq-error-rate", dest="error_mode",
                        action='store_const', const=ErrorRateMode.nucq,
                        help="error rate based on nucledite's " \
                             "Phred quality-score")

    parser.add_argument("--read-mm-error-rate", dest="error_mode",
                        action='store_const', const=ErrorRateMode.read_err_rate,
                        help="error rate based on read's match/mismatch ratio")


    parser.set_defaults(prob_matrix=ProbabilityMatrix.new)
    parser.add_argument("--prob-matrix-new", dest="prob_matrix",
                        action='store_const',const=ProbabilityMatrix.new,
                        help="use new probability matrix")
    parser.add_argument("--prob-matrix-old", dest="prob_matrix",
                        action='store_const',const=ProbabilityMatrix.old,
                        help="use old probability matrix")


    parser.add_argument("--initial-prior", type=float, default=0.00001,
                        help="initial prior probability. default: %(default)s")

    parser.add_argument("--output-basename", action="store_true",
                        help="print only the basename (without directory)" \
                              "of the 23-and-Me data files.")

    parser.add_argument('--output-format-help', action=OutputFormatHelpAction,
                        help="information about output format")

    parser.add_argument("--debug-genotypes", action="store_true",
                        help="print additional output fields (exposing " \
                        "genotype information about 23-and-Me refs")

    parser.add_argument("--cutoff-time", metavar="SECONDS", type=int,
                        default=3600,
                        help="stop processing upon SNP with " \
                        "arrival_time_relative bigger than this value. " \
                        "default: %(default)s")

    parser.add_argument("--cutoff-low-probability", metavar="P", type=float,
                        default=0.000000001,
                        help="stop processing when posterior probablity dips " \
                        "below this value. Default: %(default)g")

    parser.add_argument("--no-header", action="store_true",
                        help="do not print the header line")

    parser.add_argument('--only-header', action=HeaderOnlyAction,
                        help="Print the header and exit (without "\
                        "--debug-genotype fields)")

    # Positional parameter
    parser.add_argument('sample',
                        help='Sample file to compare (results of calling SNPS' \
                             'on minION fast5 files using XXX');

    # one-or-more reference 23-and-Me files to compare against
    parser.add_argument('datafiles',metavar='23-and-Me', nargs="+",
                        help='23-and-Me files to match against');

    args = parser.parse_args()

    # Validate error-rate parameters and combinations
    if args.error_mode == ErrorRateMode.fixed:
        if args.fixed_error_rate is None:
            args.fixed_error_rate = 0.15 # set default if not specified by user

        # validate values
        if args.fixed_error_rate < 0 or args.fixed_error_rate > 1:
            sys.exit("error: invalid fixed error rate '%s'" \
                     % args.fixed_error_rate)
    else:
        if not (args.fixed_error_rate is None):
            sys.exit("error: don't mix error-rate options (use only one)")

    return args


class SNPMatchBayesianModel:
    """Class Holding the intermediate (and final) results
    of repetitive bayesian model updates.
    """

    def __init__(self,initial_prior, probability_matrix):
        """Create a new model, set prior to 'initial-prior' value"""
        self.prior = initial_prior
        self.fam_prior = initial_prior
        self.posterior = initial_prior
        self.fam_posterior = initial_prior

        self.homo_hit = 0
        self.homo_miss = 0
        self.hetero_hit = 0
        self.hetero_miss = 0

        self.hetero_common = 0
        self.hetero_rare = 0

        self.p_obs = None
        self.fam_p_obs = None
        self.p_obs_given_match = None
        self.p_obs_given_mismatch = None
        self.p_first_deg = None

        self.last_snp_match = None

        self.probability_matrix = probability_matrix


    def update(self, new_base, ref_bases, new_base_af, error_rate):
        """Update the bayesian model.

        new_base  = the new SNP (a single nucleotide) (from the minION sample)
        ref_bases = a string of one or two SNPs (from the 23-and-Me ref file)
        new_base_af  = base's allele-frequency (=population probability)
        error_rate = the error rate attached to new_base
        """

        if self.probability_matrix == ProbabilityMatrix.new:
            f_a = new_base_af  #the freq of allele OBSERVED allele in seq data.
            f_b = 1 - f_a      #the freq of allele UNobserved allele in seq data
            e = error_rate

            #Assuming the observation came from a random person:
            self.p_obs_given_mismatch = f_a * (1 - e) + f_b * e
            m= np.matrix(   [[ 1 - e],  #Hom match
                             [ 0.5  ],  #Het
                             [ e    ]]) #Hom mismatch

            #prob matrix for each mioesis
            tran = np.matrix( [[f_a,   f_b,  0    ],
                               [f_a/2, 0.5,  f_b/2],
                               [ 0,    f_a,  f_b  ]])

            tran_m = tran * m
        else:
            # Use old method with fixed values
            self.p_obs_given_mismatch = new_base_af
            self.p_first_deg = self.p_obs_given_mismatch
            e = error_rate
            m= np.matrix(   [[ 1 - e],  #Hom match
                             [ 0.5  ],  #Het
                             [ e    ]]) #Hom mismatch
            tran_m = m

        # Update prior values: from posterior of previous round
        self.prior = self.posterior
        self.fam_prior = self.fam_posterior


        # Determine 'p_obs_given_match'
        if (len(ref_bases)==1 or ref_bases[0]==ref_bases[1]): # homozygous case
            if new_base==ref_bases[0]: # match
                self.p_obs_given_match = m[0,0]
                self.p_first_deg = tran_m[0,0]

                self.homo_hit += 1
                self.last_snp_match = "hom_hit"
            else:
                self.p_obs_given_match = m[2,0]
                self.p_first_deg = tran_m[2,0]

                self.homo_miss += 1
                self.last_snp_match = "hom_miss"

        else: # heterozygous
            if new_base in ref_bases:
                self.p_obs_given_match = m[1,0]
                self.p_first_deg = tran_m[1,0]

                self.hetero_hit += 1
                self.last_snp_match = "het_hit"
            else:
                # Assume this is invalid case, yet it COULD happen.
                # The 23-and-Me file can have alleles which are not in dbSNP
                # ref/alt. Issue a warning, and don't update posteriors.
                msg = "warning: encountered hetero_miss, new_base='%s' "\
                      "ref_bases='%s'. Skipping." % (new_base, ref_bases)
                print >>sys.stderr,msg
                self.hetero_miss += 1
                self.last_snp_match = "het_miss"

                self.p_obs_given_match = -1
                self.p_obs_given_mismatch = -1
                self.p_first_deg = -1
                self.p_obs = -1
                self.fam_p_obs = -1
                return

            if new_base_af>.5:
                self.hetero_common += 1
            else:
                self.hetero_rare += 1

        # Update model for self-match
        self.p_obs = self.prior * self.p_obs_given_match + \
                     self.fam_prior * self.p_first_deg + \
                     ( 1 - self.prior - self.fam_prior) * self.p_obs_given_mismatch

        self.posterior = self.prior * self.p_obs_given_match / self.p_obs

        # Update model for familial match
        if self.probability_matrix == ProbabilityMatrix.old and \
           (new_base_af < 0.1 or new_base_af > 0.9):
            # Old method when Allele-frequency was extreme
            self.fam_posterior = self.fam_prior
        else:
            # New method, or Old method when 0.1<=AF<=0.9
            self.fam_posterior = self.fam_prior * self.p_first_deg / self.p_obs



class CalledSNPRead:
    """Class representing one SNP/Read record
    (the output of generate-snp-list.py).

    see generate-snp-list.py --help for details.

 1.  dbSNP chromosome
 2.  dbSNP Position (1-based)
 3.  dbSNP ID (e.g. 'rsID12345')
 4.  dbSNP ref-base
 5.  dbSNP alt-base
 6.  dbSNP ref-allele-freq
 7.  dbSNP alt-allele-freq
 8.  seq/read ID (from BEDSeq,FASTQ files)
 9.  base (nucleotide found in the read for this position)
 10. numeric sanger quality-score (Phred) for this nucleotide, 0-93
     (0 = ASCII 33 '@' in FASTQ file)
 11. Allele-Frequency for this genotype
     either same as field 6 / 7, or 0 if field 10 != field 6/7.
 12. Arrival-time-Absolute: (unix-time): arrival time of read from minION
 13. Arrival-time-relative: (seconds):   seconds-since-experiment-start
 14. Arrival-number (1 = first to arrive. Gaps possible due to skipped-reads)
 15. read match/mismatch error rate (from nt_mm_err_rate in BESeq)

    """
    fields = ['chrom', 'pos', 'snp_id',
              'snp_ref_base', 'snp_alt_base',
              'snp_ref_af', 'snp_alt_af',

              'read_id', 'nuc', 'qual',
              'af',
              'arv_time_abs','arv_time_rel','arv_num',
              'mm_err_rate'
             ]
    def __init__(self, *args):
        for i in zip(self.fields,args):
            if i[0] in ['pos', 'arv_time_rel']:
                setattr(self, i[0], int(i[1]))
            elif i[0] in ['af','qual','mm_err_rate']:
                setattr(self, i[0], float(i[1]))
            else:
                setattr(self, i[0], i[1])


def load_called_snps_file(filename):
    """Loads a list of Called SNPs/Reads
    (the output of 'generate-snp-list.py').

    returns a list-of-CalledSNPRead objects.
    """
    try:
        snps = []

        f = open(filename,'r')
        for linenum,line in enumerate(f):
            err = "input error in '%s' line %d: " % (filename, linenum+1)
            line = line.strip()
            fields = line.split('\t')
            if len(fields) != 15:
                sys.exit(err + "expecting 15 fields, found %d" % (len(fields)))

            if linenum==0:
                if fields[0] != "chrom":
                    sys.exit(err + "expecting header line (first word: " \
                             "'chrom') - is this a snps file?")
                continue

            try:
                s = CalledSNPRead(*fields)
            except ValueError as e:
                sys.exit(err + "invalid value (%s)" % (str(e)))

            snps.append(s)

        return snps

    except IOError as e:
        sys.exit("I/O error while reading SNP/reads file '%s': %s" \
                 % (filename, str(e)))


def get_ref_gt_description(ref,alt, gt):
    """Given a reference base, alternate base,
    and one (or two) genotype bases,
    returns the appropriate description
    (hom_ref/het/hom_alt/other).

    Example:
       get_ref_gt_description("A","T","A")  => "hom_ref"
       get_ref_gt_description("A","T","T")  => "hom_alt"
       get_ref_gt_description("A","T","AA") => "hom_ref"
       get_ref_gt_description("A","T","TA") => "het"
       get_ref_gt_description("A","T","AT") => "het"
       get_ref_gt_description("A","T","TT") => "hom_alt"
       get_ref_gt_description("A","T","TC") => "other"
    """
    if len(gt)==1 or \
       (len(gt)==2 and gt[0] == gt[1]):
        if gt[0]==ref:
            return "hom_ref"
        elif gt[0]==alt:
            return "hom_alt"
        else:
            return "other"

    if len(gt)!=2:
        # TODO: handle it in a friendlier way
        raise Exception("input error: len(gt)!=2 (gt='%s')" % gt)

    # if we got here, genotype is het
    if (gt[0]==ref and gt[1]==alt) or \
       (gt[0]==alt and gt[1]==ref):
        return "het"

    return "other"



def load_23andMe_snps(filename):
    """
    Loads a 23-and-Me file,
    returns a dictionary of snps->genotypes.

    throws 'ProcessingException' upon errors.
    """
    seqData={}
    try:

        f=open(filename,'r')
        for linenum,line in enumerate(f):
            err = "input error in %s:%s: " % (filename, linenum+1)

            if line[0]=='#':
                continue

            line = line.strip()

            fields = line.split('\t')
            if len(fields)!=4:
                msg = "%s expecting 4 fields, got %d" % (err, len(fields))
                raise ProcessingException(msg)


            [snp, ch, pos, val] = fields

            if not (len(val)==1 or len(val)==2):
                msg = "%s invalid genotype value %s" % (err, val)
                raise ProcessingException(msg)

            # XXX TODO: validate content of val, preferably not with regex

            if snp in seqData:
                msg = "%s duplicated SNP-ID %s" % (err, snp)
                raise ProcessingException(msg)

            seqData[snp]=val

        f.close

    except IOError as e:
        msg = "I/O error while reading '%s': %s" % (filename, str(e))
        raise ProcessingException(msg)

    return seqData




def calcProb(sample_snps,
             ref_snps,
             ref_id,
             initPrior=None,
             error_rate_mode=None,
             fixed_error_rate=None,
             probability_matrix=None,
             debug_genotypes=False,
             cutoff_seconds=None,
             cutoff_low_prob=None):
    """
    Calculate probabilities of a genetic match,
    comparing the SNPs from the sample (sample_snps),
    against the SNPs from a reference 23-and-Me file (ref_snps).

    NOTE:
    This assumes 'sample_snps' are sorted by arrival-time
    (or other criterial, e.g. arrival-number) -
    as the probabilities are accumulated, and printed per sample-snp.
    """
    ## XXX TODO: if these are not used, remove them.
    heterohits=0
    sexhits=0
    misses=0

    count_snp_num = 0
    count_snps_left = len(sample_snps)
    count_snp_invalid_allele = 0
    count_snp_not_in_ref = 0
    count_snp_processed = 0

    heterorare=0
    heterocommon=0

    #p = [initPrior, initPrior, initPrior]
    #famP = [initPrior, initPrior, initPrior]
    #notes= [defaultdict(lambda:0), defaultdict(lambda:0), defaultdict(lambda:0)]

    model = SNPMatchBayesianModel(initPrior, probability_matrix)

    for snp in sample_snps:
        count_snp_num += 1
        count_snps_left -= 1

        # Terminate early based on SNP timing
        # (e.g. don't process SNPs arriving after an hour from
        #  experiment-start-time)
        if (cutoff_seconds is not None) and \
           snp.arv_time_rel > cutoff_seconds:
            break

        # AF (Allele-Frequency) will be zero if the nucleotide
        # coming from minION does not match any of the alleles in the
        # SNP database (e.g. rs12345's alleles are 'A/T' and the sequenced
        # nucleotide was 'G'.
        if snp.af==0:
            count_snp_invalid_allele +=1
            continue

        # SNP-ID from the SNP database, but not mentioned in this 23-and-Me file.
        if snp.snp_id not in ref_snps:
            count_snp_not_in_ref += 1
            #print "snp not found in ref: %s" % snp
            continue

        count_snp_processed += 1

        if error_rate_mode == ErrorRateMode.fixed:
            error_rate = fixed_error_rate
        elif error_rate_mode == ErrorRateMode.nucq:
            error_rate = pow(10,(snp.qual/-10))
        elif error_rate_mode == ErrorRateMode.read_err_rate:
             error_rate = snp.mm_err_rate

        model.update(snp.nuc, ref_snps[snp.snp_id], snp.af, error_rate)

        # Terminate after dipping below a cutoff probability
        if (cutoff_low_prob is not None) and \
           model.posterior < cutoff_low_prob:
            break

        # Output fields (explicitly convert floating points to
        # non-scientific notation strings)
        results = [
            ref_id,                     #1
            snp.read_id,                #2
            snp.snp_id,                 #3
            snp.arv_time_abs,           #4
            snp.arv_time_rel,           #5
            snp.arv_num,                #6

            count_snp_num,              #7
            count_snps_left,            #8
            count_snp_processed,        #9
            count_snp_invalid_allele,   #10
            count_snp_not_in_ref,       #11

            model.posterior,            #12
            model.fam_posterior,        #13

            model.prior,                #14
            model.fam_prior,            #15

            model.p_obs,                #16
            model.fam_p_obs,            #17
            model.p_obs_given_match,    #18
            model.p_obs_given_mismatch, #19
            model.p_first_deg,          #20

            model.homo_hit,             #21
            model.homo_miss,            #22
            model.hetero_hit,           #23
            model.hetero_miss,          #24
            model.hetero_common,        #25
            model.hetero_rare,          #26
            0,                          #27 (sex_hits = zero)

            error_rate                  #28
        ]

        if debug_genotypes:
            # returns 'hom_ref/het/hom_alt/other'
            ref_gt_desc = get_ref_gt_description(snp.snp_ref_base,
                                   snp.snp_alt_base,
                                   ref_snps[snp.snp_id])

            results.extend([
                snp.snp_ref_base,      #29
                snp.snp_alt_base,      #30
                snp.snp_ref_af,        #31
                snp.snp_alt_af,        #32
                snp.nuc,               #33
                snp.af,                #34
                ref_snps[snp.snp_id],  #35
                ref_gt_desc,           #36
                model.last_snp_match   #37
                ])

        print '\t'.join(map(str,results))





if __name__ == "__main__":
    args = parse_command_line()

    if args.debug:
        print >>sys.stderr, "Sample file:", args.sample
        print >>sys.stderr, "Initial prior:", args.initial_prior
        print >>sys.stderr, "Cutoff Seconds: ", args.cutoff_time
        print >>sys.stderr, "Cutoff low probability: ", args.cutoff_low_probability
        print >>sys.stderr, len(args.datafiles), "file(s) to process"

    # Load the Called SNPs from the sample file
    sample_snps = load_called_snps_file(args.sample)

    if not args.no_header:
        print_header(args.debug_genotypes)

    ok = True
    for f in args.datafiles:
        if args.debug:
            print >>sys.stderr,"Processing 23andMe data file: ", f

        try:
            ref_snps = load_23andMe_snps(f)
        except ProcessingException as e:
            print >>sys.stderr, "failed to load 23-and-Me file data '%s': %s" \
                % (f, str(e))
            ok = False

            # Skip this file
            continue

        # The ID to identify the sample
        # (will be printed as the first column, e.g. openSNP/DNA.Land sample_id)
        ref_id = f
        if args.output_basename:
            ref_id = basename(ref_id)

        calcProb(sample_snps, ref_snps, ref_id,
                 initPrior=args.initial_prior,
                 error_rate_mode=args.error_mode,
                 fixed_error_rate=args.fixed_error_rate,
                 probability_matrix=args.prob_matrix,
                 debug_genotypes=args.debug_genotypes,
                 cutoff_seconds=args.cutoff_time,
                 cutoff_low_prob=args.cutoff_low_probability)

    if not ok:
        sys.exit("error: processing one or more files failed")
