#!/usr/bin/env python

# TODO: find a more accurate name,

# TODO: document output format.

version_info=\
"""
generate-snp-list.py: part of the Personal Identification Pipeline
https://github.com/TeamErlich/personal-identification-pipeline

Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)
All Rights Reserved.
This software is restricted to educational, research, not-for-profit purposes.
See LICENSE file for full details.

Version 0.3
"""

import sys, argparse, os
import pysam
from pprint import pprint


def parse_command_line():
    # Define parameters
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Personal-ID Pipeline: SNP-Caller",
        version=version_info,
        epilog="""
MinION Personal ID Pipeline - SNP "Caller"

This script generates a list of SNPs (with alleles and probabilities)
from a (pre-processed) minION sequencing run.

The SNP database must be preprocessed with bgzip/tabix.
See 'setup-snp138common.sh' for details.

See sam-to-bedseq.py for BEDSeq information.

Output Format:
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
        """)

    parser.add_argument('--verbose', dest='verbose', action='store_true',
                        help="show progress information on STDERR")

    parser.add_argument('bedseq',
                        metavar='BEDSeq',
                        help='input BEDSeq file');

    parser.add_argument('SNPdb', metavar='SNPdb',
                        help='file with SNP database (in UCSC GB ' \
                             'format of snp138common)');

    args = parser.parse_args()

    # check for tabix index
    if not os.path.exists(args.SNPdb) or \
       not os.path.exists(args.SNPdb+".tbi"):
        sys.exit("SNP database (%s) or tabix index (%s.tbi) not found. " \
                 "please see 'setup-snp138common.sh' for building instructions." \
                 % (args.SNPdb,args.SNPdb))

    return args


class BEDSeqRead:
    """Class representing one read from a BEDSeq file.

    BEDSeq fields are (see sam-to-bedseq.py):
    1. chrom
    2. start-pos (0-based)
    3. end-pos   (1-based)
    4. read-id
    5. dummy (always zero)
    6. strand    (+ or -)
    7. absolute arrival time       (seconds)
    8. relative arrival time       (seconds)
    9. arrival number              (1 = first read to arrive from minION)
    10. original sequence length   (in FASTQ/timing file, as came of the minION)
    11. mapped sequence length     (before normalization, as mapped by BWA)
    12. normalized sequence length (= end-pos - start-pos)
    13. sam_nm                     (number of unmatching bases, from SAM's NM tag)
    14. nt_align                   (#nuc that aligned to ref, match or mismatch)
    15. nt_match                   (#nuc that align+match the reference)
    16. nt_mismatch                (#nuc that align+mismatch, e.g snps)
    17. nt_ins                     (#nuc that are insertions)
    18. nt_del                     (#nuc that are deletions)
    19. nt_sclip                   (#nuc that are soft-clipped)
    20. nt_hclip                   (#nuc that are hard-clipped)
    21. nt_mm_err_rate             ( 1 - (nt_match/nt_align) )
    22. CIGAR                      (CIGAR value from SAM file)
    23. MD                         (MD tag from SAM file)
    24. normalized sequence        (ALWAYS forward strand)
    25. normalized quality scores
    """
    fields = ['chr', 'start', 'end', 'id', 'dummy',  'sense',
              'arv_time_abs', 'arv_time_rel', 'arv_num',
              'orig_seq_len', 'mapped_seq_len','norm_seq_len',
              'sam_nm',
              'nt_align', 'nt_match','nt_mismatch',
              'nt_ins','nt_del','nt_sclip','nt_hclip',
              'nt_mm_err_rate',
              'cigar','sam_md',
              'sequence','readqual']
    def __init__(self, *args):
        for i in zip(self.fields,args):
            if i[0] in ['start', 'end']:
                # NOTE: There are other numeric fields, but they are
                # not currently used numerically.
                setattr(self, i[0], int(i[1]))
            elif i[0] in ['nt_mm_err_rate']:
                setattr(self, i[0], float(i[1]))
            else:
                setattr(self, i[0], i[1])

    def __repr__(self):
        return "<BEDSeqRead(id='%s' %s:%d-%d arv_num=%s>" % \
            (self.id, self.chr, self.start, self.end, self.arv_num)


class SNP:
    """Class representing one SNP record.
    The fields match the data file from the UCSC Genome Browser,
    track 'snp138common.txt',
    see: http://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&c=chr2&g=snp138Common
    """
    fields = ['bin', 'chrom', 'chromStart', 'chromEnd', 'name', 'score',
              'strand', 'refNCBI', 'refUCSC', 'observed', 'molType', 'type',
              'valid', 'avHet', 'avHetSE', 'func', 'locType', 'weight',
              'exceptions', 'submitterCount', 'submitters', 'alleleFreqCount',
              'alleles', 'alleleNs', 'alleleFreqs', 'bitfields']
    def __init__(self, *args):
        for i in zip(self.fields,args):
            if i[0] in ['chromStart', 'chromEnd', 'score', 'alleleFreqCount']:
                setattr(self, i[0], int(i[1]))
            else:
                setattr(self, i[0], i[1])


class SNPTabix:
    """A Class to encapsulate fetching of SNP records based on given reads
    (chrom/start/pos)
    """

    def __init__(self, filename):
        self.tabix_file_name = filename
        # TODO: catch TABIX exceptions
        self.tabix = pysam.TabixFile(filename)


    def fetch_snps(self,chrom,start,end):
        """Returns list-of-SNP() given chrom/start/end.

        Equivalent to:
            tabix FILE chrom:start-end

        Also validates each returned item is a valid SNP record.

        Aborts on input errors.

        TODO: consolidate error messages.
        """
        snps = []
        for row in self.tabix.fetch(chrom, start, end):
            flds = row.strip().split()
            if len(flds) < 26:
                print >>sys.stderr, \
"""tabix/snps error: expecting 26 field, got %s
when fetching region %s:%s-%s .
bad line:
  %s

To troubleshoot manually, try:
  tabix %s %s:%s-%s

""" % ( len(flds), read.chr, read.start, read.end, row,
        self.tabix_file_name, read.chr, read.start, read.end)
                for i,v in enumerate(flds):
                    print >>sys.stderr, "field[%s] = '%s'" % (i+1, v)
                sys.exit("bad SNPs database. aborting.")

            try:
                snp = SNP(*flds)
            except ValueError as e:
                print >>sys.stderr, \
"""tabix/snps error: one of the snp-database fields is invalid
when fetching region %s:%s-%s .
error was:
  %s
bad line:
  %s

To troubleshoot manually, try:
  tabix %s %s:%s-%s

""" % ( read.chr, read.start, read.end, str(e), row,
        self.tabix_file_name, read.chr, read.start, read.end)
                for i,v in enumerate(flds):
                    print >>sys.stderr, "field[%s] = '%s'" % (i+1, v)
                sys.exit("bad SNPs database. aborting.")

            if snp.strand != '+':
                continue


            if snp.alleleFreqCount != 2:
                msg = "warning: SNP '%s' has allelFreqCount of '%s' " \
                      "(only '2' supported)" % (snp.name, snp.allelFreqCount)
                print >>sys.stderr,msg
                continue

            # Build allele-frequency dictionary
            # e.g.
            #   alleles = "A,C,"
            #   alleleFreqs = "0.1,0.9,"
            #   => { 'A':0.1, 'C':0.9 }
            snp_bases = snp.alleles.split(',')
            snp_bases = [x for x in snp_bases if x]
            snp_freqs = snp.alleleFreqs.split(',')
            snp_freqs = [x for x in snp_freqs if x]
            if len(snp_bases) != 2 or len(snp_freqs) != 2:
                msg = "warning: SNP '%s': alleles/alleleFreqss counts!=2"\
                      % (snp.name)
                print >>sys.stderr,msg
                continue

            d = dict(zip(snp_bases, snp_freqs))
            setattr(snp, 'allele_freq_dict', d)

            # Find the ALT base, and save it
            alt = [x for x in d if x != snp.refUCSC]
            if len(alt)!=1:
                # Something is terribly wrong:
                # expected just two items, and after removing the ref base,
                # only ALT nucleotide left
                msg = "warning: SNP '%s' error calculating ALT base/af" % \
                      (snp.name)
                print >>sys.stderr,msg
                continue
            alt=alt[0]
            setattr(snp,'alt_nuc',alt)


            snps.append(snp)


        return snps



def BEDSeq_snp_caller(bedseq, snp_db_tabix,progress=False):
    """
    Reads a BEDSeq file (BED positions + normalize sequences),
    then compares each nucleotide against overlapping SNPs.

    For each match, writes the SNP information + the minION read it matched.
    """
    try:
        # print header
        out = [
            "chrom",            #1
            "pos",              #2
            "snp_id",           #3
            "snp_ref_base",     #4
            "snp_alt_base",     #5
            "snp_ref_af",       #6
            "snp_alt_af",       #7

            "seq_id",           #8
            "base",             #9
            "qual",             #10
            "af",               #11
            "arv_time_abs",     #12
            "arv_time_rel",     #13
            "arv_num",          #14
            "mm_err_rate",      #15
          ]
        print '\t'.join(out)

        snptabix = SNPTabix(snp_db_tabix)

        inf = open(bedseq,'r')
        for linenum,line in enumerate(inf):
            err = "input error in '%s' line %d: " % (bedseq, linenum+1)

            ##
            ## Read & Parse a line from the BEDSeq file
            ##
            fields = line.strip().split('\t')
            if len(fields) != 25:
                sys.exit(err + "expecting 24 fields, found %d" % (len(fields)))

            if linenum==0:
                # header line
                if fields[0] != "chrom":
                    sys.exit(err + "invaid header line, expecting 'chrom' " \
                             "as first field, found '%s'" % (fields[0]))
                continue

            try:
                read = BEDSeqRead(*fields)
            except ValueError as e:
                # Not terribly helpful, but better than nothing
                sys.exit(err + "invalid values: %s" % (str(e)))

            if progress:
                print >>sys.stderr, "Calling SNPs on read %s (%s:%s-%s)" % \
                                    (read.id,read.chr,read.start,read.end)

            ##
            ## Find overlapping SNP
            ## (using tabix from the pre-processed snps db file)
            snps = snptabix.fetch_snps(read.chr, read.start, read.end)

            ##
            ## Compare each nucleotide in the read against the SNPs
            ##
            snp_names = []
            for snp in snps:
                if read.start > snp.chromStart or read.end <= snp.chromEnd:
                    continue

                nuc = read.sequence[snp.chromStart-read.start]
                # '|' is a marker in BEDSeq for a nucleotide that was 'normalized'
                # and did not exist in the original sequence (e.g. a deletion).
                if nuc == '|':
                    continue

                rq = read.readqual[snp.chromStart-read.start]
                rq = ord(rq)-33

                # ref/alt info from dbSNP
                ref_base = snp.refUCSC
                ref_af = snp.allele_freq_dict[ref_base]
                alt_base = snp.alt_nuc
                alt_af = snp.allele_freq_dict[alt_base]

                # Find allele-freq of the sequenced/mapped nucleotide
                # from the BEDSeq input:
                af = snp.allele_freq_dict.get(nuc,0)


                out=[
                    snp.chrom,         #1
                    snp.chromStart+1,  #2  1-base SNPs, like 23-and-Me
                    snp.name,          #3
                    ref_base,          #4
                    alt_base,          #5
                    ref_af,            #6
                    alt_af,            #7

                    read.id,           #8
                    nuc,               #9
                    rq,                #10
                    af,                #11
                    read.arv_time_abs, #12
                    read.arv_time_rel, #13
                    read.arv_num,      #14
                    read.nt_mm_err_rate  #15
                ]
                snp_names.append(snp.name)

                print '\t'.join(map(str,out))

            if progress:
                print >>sys.stderr,"  found %d snps: %s" \
                       % (len(snp_names), ','.join(snp_names))

    except IOError as e:
        sys.exit("error reading file-list from '%s':" % (str(e)))




if __name__ == "__main__":
    args = parse_command_line()

    reads = BEDSeq_snp_caller(args.bedseq, args.SNPdb,progress=args.verbose)
