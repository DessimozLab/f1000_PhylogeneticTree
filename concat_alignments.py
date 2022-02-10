#!/usr/bin/env python

import logging
import sys
import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq, UnknownSeq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict


logger = logging.getLogger("concat")


def load_alignments(alignmentfiles, format):
    alignments = []
    for file in alignmentfiles:
        try:
            for alignment in AlignIO.parse(file, format=format):
                logger.debug("loaded alignment of length {} from {}".format(len(alignment), file))
                alignments.append(alignment)
        except ValueError as e:
            logger.error("Cannot parse input file {}: {}".format(file, e))
            raise
    logger.info("Successfully loaded {} alignments from {} input files"
                .format(len(alignments), len(alignmentfiles)))
    return alignments


def concatenate(alignments):
    # Get the full set of labels (i.e. sequence ids) for all the alignments
    all_labels = set(seq.id for aln in alignments for seq in aln)
    logger.debug("extracted {} different labels in all alignments: {}"
                 .format(len(all_labels), all_labels))

    # Make a dictionary to store info as we go along
    # (defaultdict is convenient -- asking for a missing key gives back an empty list)
    concat_buf = defaultdict(list)

    for aln in alignments:
        length = aln.get_alignment_length()

        # check if any labels are missing in the current alignment
        these_labels = set(rec.id for rec in aln)
        missing = all_labels - these_labels
        logger.debug("alignment of length {} with {} sequences, {} missing ({})"
                     .format(length, len(these_labels), len(missing), missing))

        # if any are missing, create unknown data of the right length,
        # stuff the string representation into the concat_buf dict
        for label in missing:
            new_seq = UnknownSeq(length, character='X')
            concat_buf[label].append(str(new_seq))

        # else stuff the string representation into the concat_buf dict
        for rec in aln:
            concat_buf[rec.id].append(str(rec.seq))

    # Stitch all the substrings together using join (most efficient way),
    # and build the Biopython data structures Seq, SeqRecord and MultipleSeqAlignment
    msa = MultipleSeqAlignment(SeqRecord(Seq(''.join(seq_arr)), id=label)
                                for (label, seq_arr) in concat_buf.items())
    logger.info("concatenated MSA of {} taxa and total length {} created"
                .format(len(msa), len(msa[0])))
    return msa


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate alignments",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--format-in', default='fasta',
                        help="input format of the alignments. Any format that is understood"
                             "by Biopython's AlignIO module is possible.")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                        help="Path to the output file where the concatenated multiple "
                             "sequence alignment will be written")
    parser.add_argument('-u', '--format-output', default='phylip-relaxed',
                        help="output format of the concatenated multiple sequence alignment")
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help="Produce some output and status reports")
    parser.add_argument('-d', '--debug', action="store_true", default=False,
                        help="Be more verbose for debugging purposes")
    parser.add_argument('alignment', nargs='+', type=str,
                        help="Path to the alignment files. Use shell expansion to pass many files "
                             "in a simple way, e.g. \"/path/to/folder/*.fa\".")
    conf = parser.parse_args()

    level = logging.WARNING
    if conf.verbose:
        level = logging.INFO
    if conf.debug:
        level = logging.DEBUG
    logging.basicConfig(level=level, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logger.debug("Concatenate alignments: arguments: {}".format(conf))

    alignments = load_alignments(conf.alignment, conf.format_in.lower())
    msa = concatenate(alignments)
    AlignIO.write(msa, conf.output, conf.format_output.lower())
