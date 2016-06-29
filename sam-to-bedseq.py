#!/usr/bin/env python
from __future__ import print_function


version_info=\
"""
Sam-To-BEDSeq: part of the Personal Identification Pipeline
https://github.com/TeamErlich/personal-identification-pipeline

Copyright (C) 2016 Yaniv Erlich (yaniv@cs.columbia.edu)
All Rights Reserved.
This software is restricted to educational, research, not-for-profit purposes.
See LICENSE file for full details.

Version 0.3
"""

import sys,argparse, re
from pprint import pprint

class BEDSeqHelpAction(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super(BEDSeqHelpAction, self).__init__(option_strings, '', nargs=0, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        text=version_info + """

BEDSeq output Format
====================

This script reads a SAM file, and writes a BED file
with additional columns of normalized sequence+quality scores.

The output columns are:
   1. chrom
   2. start-pos (0-based)
   3. end-pos   (1-based)
   4. read-id
   5. dummy (always zero)
   6. strand    (+ or -)
   7. absolute arrival time      (seconds)
   8. relative arrival time      (seconds)
   9. arrival number             (1 = first read to arrive from minION)
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

Arrival Time, Number, original sequence lengths (columns 7, 8, 9,10):
  Will be calculated based on matching the sequence IDs to the
  timing file. If no timing file is given, these will be -1 .

  arrival_time_abs: taken from the 'unix_timestamp_end' value of
                    the read from the timing file.

  arrival_time_rel: delta between 'unix_timestamp_end' and 'exp_starttime'
                    from the timing file.

  arrival_num:      sequencial order of arrival of the read.
                    NOTE: gaps are possible, as this script counts
                    ALL reads from the timing file, but writes
                    only MAPPED reads from the SAM file.
                    e.g. if reads num 1 and 3 are in the output BEDSeq
                    file, it means the second read to arrive of the minION
                    did not map and was skipped.

  orig_seq_len:     The length of the read in the FASTQ file off the minION.
                    If the read was hard-clipped by the genome mapping program,
                    the mapped length will be shorter.

Sequence Normalization (columns 24, 25):
  Insertions in the SAM file are discarded (as they don't exist in the
  reference genome), and deletions are padded with '|' (a character
  that is invalid for both nucleotides and quality scores).

  Given the following SAM input:
    sequence:  AACCGGAT
    quality:   %%%%%%%%
    CIGAR:     2M2D6M

  The normalized sequence will be 'AA||CCGGAT',
  The normalized quality  will be '%%||%%%%%%'.

  This ensures that the coordinate of all nucleotides match the
  reference genome (e.g. the position of 'T' is 8th in the SAM sequence,
  but 10th in the reference genome. It will be 10th in the
  output FASTQ file as well.

The following should hold true (if not - it's a bug):
  nt_align   = nt_match + nt_mismatch
  sam_nm     = nt_mismatch + nt_ins + nt_del
"""
        print (text)
        sys.exit(0)


def parse_command_line():
    # Define parameters
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="SAM-to-BEDSeq script",
        version=version_info,
        epilog="""
This script reads a SAM file, and writes a BEDSeq file
(a BED file with additional columns of minION timing and
 normalized sequence+qualityscores).

example:

    # Generate FASTQ,TIMES by poretools.
    $ poretools fastq /server/data/myruns/test1/ > test1.fastq
    $ poretools times /server/data/myruns/test1/ > test1.times

    # Map reads
    $ bwa mem test1.fastq hg19.fa > test1.sam

    # Generate a BEDSeq file, combining timing information
    # and normalizing mapped sequences:
    $ %(prog)s --times test1.times test1.sam > test1.bedseq

NOTE:
1. non-mappers are ignored.

See --bedseq-help for more details about the output format.

        """)

    parser.add_argument('--bedseq-help', action=BEDSeqHelpAction,
                        help="information about BEDSeq output format")
    parser.add_argument('--times', metavar='FILE',
                        help='poretools timing file');
    parser.add_argument('--debug-times', action="store_true",
                        help='pring loaded timing information and exit')

    # Positional parameter
    parser.add_argument('sam',   metavar='SAM',   help='SAM file to read');

    args = parser.parse_args()
    return args


class ReadTiming:
    """Class representing Timing of one read, based on information from
    a 'poretools time' file.
    """
    fields = ['filename','abs_time','delta_time','seq_len']
    def __init__(self, *args):
        """
        parameters:
         filename:   filename (fast5) of the read
         delta-time: seconds since 'experiment-start-time'
         abs-time:   seconds since epoch, sequence arrival-end time
         seq-len:    length of sequence

        NOTE:
        when the ReadTiming instance is created, the 'index' is yet
        unknown (because we don't assume the input timing file is sorted)
        """
        setattr(self, 'index','None')
        for i in zip(self.fields,args):
            if i[0] == 'filename':
                setattr(self, i[0], i[1])
            else:
                setattr(self, i[0], int(i[1]))

    def set_index(self,index):
        setattr(self,'index',int(index))

    def __repr__(self):
        return "<ReadTiming(index=%s delta=%d abs=%d seq_len=%d filename=%s)>" \
               % (str(self.index), self.delta_time, self.abs_time, \
                  self.seq_len, self.filename)




def load_timing_file(filename):
    """
    Read the 'poretools times' file
    as generated by:
          poretools times DIR/ > file.times

    returns a dictionary of filename=>ReadTimes() .

    NOTE: we don't assume the input file is sorted by end_unixtime.
    """
    timing = []

    first_experiment_start_time = None

    try:
        f = open(filename,'r')
        for linenum,line in enumerate(f):
            err = "input error in '%s' line %d: " % (filename, linenum+1)
            line = line.strip()
            fields = line.split('\t')

            if len(fields)<7:
                sys.exit(err + "expecting 7 fields, found %d" % (len(fields)))

            if linenum==0:
                if fields[0] != "channel":
                    sys.exit(err + "expecting header line (first word: " \
                             "'channel') - is this a 'poretools times' file?")
                continue

            # Get the experiment start time
            try:
                experiment_start_time = int(fields[3])
            except ValueError as e:
                sys.exit(err + "invalid exp_starttime value in field 4: "\
                         "'%s'" % (fields[3]))

            if first_experiment_start_time is None:
                first_experiment_start_time = experiment_start_time
            if first_experiment_start_time != experiment_start_time:
                sys.exit(err + "timing time contains exp-start-time from " \
                         "several different runs. aborting")

            filename = fields[1]
            try:
                unix_timestamp_end = int(fields[6])
            except ValueError as e:
                sys.exit(err + "invalid unit timestamp in field 7: '%s'" \
                         % (fields[6]))
            try:
                seq_len = int(fields[2])
            except ValueError as e:
                sys.exit(err + "invalid read-length in field 3: '%s'" \
                         % (fields[2]))

            if unix_timestamp_end < experiment_start_time:
                sys.exit(err + "invalid unix timestamps: exp_starttime>end")


            timing.append( ReadTiming(
                             filename,
                             unix_timestamp_end,
                             unix_timestamp_end - experiment_start_time,
                             seq_len
                           ) )

        # Add sequencing index number to each read, based on sorted unix time
        timing = sorted(timing, key=lambda x:x.delta_time)

        # Add index (based on sorting, above) and create a dictionary.
        td = dict()
        for i,t in enumerate(timing):
            t.set_index(i)
            td[t.filename] = t

        return td

    except IOError as e:
        sys.exit("I/O error while reading times file '%s': %s" \
                 % (filename, str(e)))


def get_cigar_counts(cigar):
    """
    Given a cigar string, returns tuple counts of:
     (matches, insertions, deletions, soft_clip,hard_clips)

    NOTE:
    'matches' refers to aligned nucleotides: each one can be a match or a
    mismatch compared to the reference (see 'MD' tag to differentiate those).

    (match,ins,del,sclip,hclips) = get_cigar_counts("59M2D16M")
       => match == 65
          ins = 0
          del = 2
          sclip = 0
          hclip = 0
    """
    # initialize zero counts for select cigar tags
    counts = dict( [ (x,0) for x in 'IDMSH' ] )

    # Split CIGAR string into op-codes (pairs of number+character)
    cigarpieces = re.split('([0-9]+[MIDNSHPX=])', cigar)
    # discard empty pieces (due to the way 're.split' works
    cigarpieces = [x for x in cigarpieces if x]

    for piece in cigarpieces:
        n=int(piece[:-1])
        op=piece[-1]
        if not (op in counts):
            raise ValueError("implementation error: CIGAR op '%s' not handled "\
                             "(cigar = '%s')" % (cigar))
        counts[op] = counts[op] + n

    return ( counts['M'],
             counts['I'],
             counts['D'],
             counts['S'],
             counts['H'] )


def get_md_counts(md):
    """
    Given an MD string (from a SAM file 'MD:Z:x' tag),
    returns the number of nucleotides matched and mismatched
    (deletions/insertions are ignored).

    Example:
      (match,mismatch) = get_md_counts("59^TT0G15)
       => match => 64
          mismatch => 1

    The '^TT' is deletion and is ignored.
    """
    parts = re.split('([0-9]+)', md)
    parts = [x for x in parts if x]

    matches = 0
    mismatches = 0
    for x in parts:
        try:
            v = int(x)
            # Numeric? it's number of matched nucleotides
            matches = matches + v
        except ValueError:
            # Not numeric:
            # either deletions (with '^' prefix) or mismatches
            if not x.startswith('^'):
                mismatches = mismatches + len(x)
    return (matches, mismatches)



def normalize_cigar_seq(cigar, seq, qual):
    """Normalizes a nucleotide sequence based on CIGAR string.

    Insertions in the SEQ are discarded (as they don't exist in the
    reference genome), and deletions are padded with '|' (a character
    that is invalid for both nucleotides and quality scores).

    Raises ValueError exception on invalid input.

    Example:

      sequence:  AACCGGAT
      quality:   %%%%%%%%
      CIGAR:     2M2D6M
      The normalized sequence will be 'AA||CCGGAT'.
      The normalized quality  will be '%%||%%%%%%'.

      This ensures that the coordinate of all nucleotides match the
      reference genome (e.g. the position of 'T' is 8th in the SAM sequence,
      but 10th in the reference genome. It will be 10th in the
      output FASTQ file as well.

    Usage:

      (norm_seq,norm_qual) = normalize_cigar_seq('2M2D6M','AACCGGAT','%%%%%%%%')
    """
    if len(seq) != len(qual):
        raise ValueError("length of sequence,quality-scores differ")

    if cigar=="*":
        # non-mapper, return as-is
        return (seq,qual)

    # Split CIGAR string into op-codes (pairs of number+character)
    cigarpieces = re.split('([0-9]+[MIDNSHPX=])', cigar)

    # discard empty pieces (due to the way 're.split' works
    cigarpieces = [x for x in cigarpieces if x]

    outseq = ""
    outqual = ""
    i=0
    for piece in cigarpieces:
        n=int(piece[:-1])
        op=piece[-1]

        if op=='S':
            # Soft-Clipping - skip nucleotides from the input string
            i += n

        elif op=='I':
            # Insertions - skip nucleotides from the input string

            #Optionally: save the skipped nucleotide
            #for c in seq[i:(i+n)]:
            #    ins[c]+=1
            i += n

        elif op=='M':
            # Matches/Mismatches: add matched nucleotides

            # Python will silently ignore string slicing longer than the string
            # itself, so explicitly check for that.
            if (i+n) > len(seq):
                raise ValueError("cigar string '%s' at op '%s%s' expands " \
                                 "to longer than sequence (len=%d) " \
                                 % (cigar,n,op, len(seq)))

            outseq +=  seq[i:(i+n)]
            outqual += qual[i:(i+n)]
            i+=n
        elif op=='D':
            # Deletions: add padding
            pad = '|'*n
            outseq += pad
            outqual += pad
        elif op=='H':
            # Hard-clipping: silently ignore
            # (in hard-clipping, the nucleotides have been removed from
            # the sequence string before being written to the SAM/BAM file)
            pass
        elif op=='N' or op=='P' or op=='X':
            # TODO: implement these if we ever encounter them,
            # (or silently ignore)
            raise ValueError("implementation  error: CIGAR op '%s' not handled")
        else:
            # Should never happen (unless the cigar regex above is wrong)
            # raise generic 'exception' so it wouldn't be caught as ValueError
            raise Exception("programmer is stupid (op='%s' cigar='%s')" %
                            (op, cigar))

    return (outseq, outqual)




cigar_re = re.compile("^\*|([0-9]+[MIDNSHPX=])+$")
seq_re = re.compile("^\*|[A-Za-z=.]+$")
def validate_sam_fields(flds):
    """Given fields from a SAM file (after splitting by tab),
    raises ValueError (with detailed information)
    if something isn't valid.

    returns silently if fields are valid.

    NOTE:
    Not all fields are validated at the moment.
    """
    if len(flds)<11:
        raise ValueError("expecting at least 11 fields, found %d " \
                         % (len(flds)))

    try:
        flags = int(flds[1])
    except ValueError as e:
        raise ValueError("invalid flags '%s' - expecting numeric value" \
                         % flds[1])

    # if it's a mapper, POS must be positive integer value
    if flags & 4 == 0:
        try:
            pos = int(flds[3])
            if pos<1:
                raise ValueError("invalid position '%s' - expecting value>=1"\
                                 % flds[3])
        except ValueError as e:
            raise ValueError("invalid position '%s' - expecting numeric value" \
                             % flds[3])

    cigar=flds[5]
    if not cigar_re.search(cigar):
        raise ValueError("invalid CIGAR string '%s'" % cigar)

    seq=flds[9]
    if not seq_re.search(seq):
        raise ValueError("invalid SEQ string '%s' " % seq)

    qual=flds[10]
    if len(seq) != len(qual):
        raise ValueError("SEQ and QUAL fields differ in length " \
                         "(%d vs %d)" % (len(seq), len(qual)))


def build_sam_tags(flds):
    """
    Given a list of fields from a SAM file
    (all fields, including the first 11 fixed fields),
    returns a dictionary with the SAM tags (e.g. 'MD', 'NM').

    Tags with type 'i' are converted to integers.
    Tags with type 'f' are converted to floats.

    Example:
      tags = build_sam_tags( ["NM:i:0","MD:Z:77"])
        => tags = { 'NM':0, 'MD':'77' }
    """

    # Split tags into tuples of (name,type,value)
    # e.g. ["NM:i:0","MD:Z:77"] => [('NM', 'i', '0'), ('MD', 'Z', '77')]
    in_tags = [ tuple(x.split(':')) for x in flds[11:]]
    out_tags = {}
    for n,t,v in in_tags:
        if t=="i":
            v = int(v)
        elif t=='f':
            v = float(v)
        out_tags[n]=v
    return out_tags


def sam_to_bedseq(sam_filename, timing_dict=None):
    """Converts a SAM file to a BEDSeq file,
    after expanding the sequence/quality fields.
    """
    try:
        # Print Header
        out = [
                    "chrom",                  #1
                    "start_pos",
                    "end_pos",
                    "seq_id",
                    "dummy",
                    "strand",                 #6

                    "arrival_time_abs",       #7
                    "arrival_time_rel",       #8
                    "arrival_num",            #9

                    "orig_seq_len",          #10
                    "mapped_seq_len",        #11
                    "norm_seq_len",          #12

                    "sam_nm",                #13
                    "nt_align",              #14
                    "nt_match",              #15
                    "nt_mismatch",           #16
                    "nt_ins",                #17
                    "nt_del",                #18
                    "nt_sclip",              #19
                    "nt_hclip",              #20
                    "nt_mm_err_rate",        #21

                    "cigar",                 #22
                    "sam_md",                #23

                    "norm_seq",              #24
                    "norm_qual"              #25

              ]
        print ('\t'.join(out))

        sam=file(sam_filename,'r')

        for linenum,line in enumerate(sam):
            err = "input error in '%s' line %d: " % (sam_filename, linenum+1)
            line = line.strip()

            if line[:1]=='@':
                continue

            flds = line.split('\t')
            try:
                validate_sam_fields (flds)
                flags = int(flds[1])

                if (flags & 4)==4:
                    # non-mapper: skip it
                    continue

                cigar = flds[5]
                seq   = flds[9]
                qual  = flds[10]
                (norm_seq, norm_qual) = normalize_cigar_seq(cigar,seq,qual)

                (nt_align,nt_ins,nt_dels,nt_sclips,\
                                         nt_hclips) = get_cigar_counts(cigar)

                sam_tags = build_sam_tags(flds)

                if not ('MD' in sam_tags):
                    raise ValueError("MD tag not found in this read")
                sam_md = sam_tags['MD']

                sam_nm = sam_tags.get('NM',0)

                (nt_match,nt_mismatch) = get_md_counts(sam_md)

                nt_mm_err_rate = 1.0 - ((1.0 * nt_match) / nt_align)
            except ValueError as e:
                sys.exit(err + str(e))

            # Construct BED Output
            chrom =  flds[2]
            seq_id = flds[0]

            start_pos = int(flds[3]) - 1 # 1-based SAM to 0-based BED
            end_pos = start_pos + len(norm_seq)+1

            positive_orientation = (flags & 16) == 0
            strand = "+"
            if not positive_orientation:
                strand = "-"

            # Get arrival/timing information for this read
            arrival_time_abs = -1
            arrival_time_rel = -1
            arrival_num      = -1
            input_seq_len    = -1
            if timing_dict:
                if not seq_id in timing_dict:
                    sys.exit(err + "sequence ID '%s' not found in timing file" \
                             % (seq_id))

                readtiming = timing_dict[seq_id]
                arrival_time_rel = readtiming.delta_time
                arrival_time_abs = readtiming.abs_time
                arrival_num  = readtiming.index
                input_seq_len = readtiming.seq_len


            out = [
                    chrom,                  #1
                    str(start_pos),
                    str(end_pos),
                    seq_id,
                    '0',
                    strand,                 #6

                    str(arrival_time_abs),  #7
                    str(arrival_time_rel),  #8
                    str(arrival_num),       #9

                    str(input_seq_len),     #10
                    str(len(seq)),          #11
                    str(len(norm_seq)),     #12

                    str(sam_nm),            #13
                    str(nt_align),          #14
                    str(nt_match),          #15
                    str(nt_mismatch),       #16
                    str(nt_ins),            #17
                    str(nt_dels),           #18
                    str(nt_sclips),         #19
                    str(nt_hclips),         #20
                    str(nt_mm_err_rate),    #21

                    cigar,                  #22
                    sam_md,                 #23

                    norm_seq,               #24
                    norm_qual               #25
                  ]

            print ('\t'.join(out))

            # Sanity checks
            # (after printing, so that if there's an error, the last printed
            #  line contains the offending data)
            if nt_align != (nt_match + nt_mismatch):
                msg = "Sanity-check failed: nt_align(%s) != nt_match(%s) " \
                      " + nt_mismatch(%s)" % (nt_align,nt_match,nt_mismatch)
                sys.exit(err + msg)
            if sam_nm != (nt_ins + nt_dels + nt_mismatch):
                msg = "Sanity-check failed: sam_nm(%s) != nt_ins(%s) " \
                      " + nt_dels(%s) + nt_mismatch(%s)" \
                      % (sam_nm,nt_ins,nt_dels,nt_mismatch)
                sys.exit(err + msg)

    except IOError as e:
        if e.filename == sam_filename:
            sys.exit("I/O error while reading SAM file: %s" % (str(e)))

        # perhaps output I/O error?
        sys.exit("I/O error: %s" % (str(e)))




if __name__ == "__main__":
    args = parse_command_line()

    timing = None
    if args.times:
        timing = load_timing_file(args.times)
        if args.debug_times:
            pprint(timing)
            sys.exit()

    sam_to_bedseq(args.sam,timing)
