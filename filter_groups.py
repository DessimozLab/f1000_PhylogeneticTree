#!/usr/bin/env python

import math
import multiprocessing
import argparse
import os
import Bio.SeqIO
import re
import logging
logger = logging.getLogger('filter-group')
RE_SPECIES = re.compile("\[(?P<species>[^]]*)\].*")


def get_species_from_header(header):
    m = RE_SPECIES.search(header)
    if m is not None:
        return m.group('species')
    raise InvalidHeaderFormat(header)


def convert_fasta_file(infile, outfile, min_nr_species):

    def convert_id_to_species_name(record):
        record.id = get_species_from_header(record.description)
        if record.id in seen:
            raise SeveralSequencePerSpeciesException("Species \"{}\" seen more than once in {}"
                                                     .format(record.id, infile))
        seen.add(record.id)
        return record

    seen = set([])
    records = [convert_id_to_species_name(z) for z in Bio.SeqIO.parse(infile, 'fasta')]
    if len(records) >= min_nr_species:
        cnts = Bio.SeqIO.write(records, outfile, 'fasta')
        logger.debug("\"{}\" contains {} sequences. IDs have been converted and file rewritten to {}"
                     .format(infile, cnts, outfile))
    else:
        logger.debug("skipping \"{}\" with {} sequences".format(infile, len(records)))
        cnts = 0
    return cnts


def get_matching_files(input_dir, pattern):
    with os.scandir(input_dir) as dir_iter:
        for f in dir_iter:
            if f.name.endswith(pattern):
                yield f


def convert_files(input_dir, output_dir, min_nr_species, pattern='.fa', nr_processes=None):
    os.makedirs(output_dir, exist_ok=True)
    args_iter = ((x.path, os.path.join(output_dir, x.name), min_nr_species) for x in get_matching_files(input_dir, pattern))
    with multiprocessing.Pool(processes=nr_processes) as mp:
        logger.info("converting files in {} with {} processes in parallel"
                    .format(input_dir, nr_processes if nr_processes else os.cpu_count()))
        result = mp.starmap(convert_fasta_file, args_iter)
        converted = [x for x in result if x > 0]
    logger.info("{} (of {}) files filtered and converted. Converted files "
                "contain {:.2f} species on average (sd {:.4f})"
                .format(len(converted), len(result), sum(converted)/len(converted),
                        math.sqrt((sum(x**2 for x in converted) - (sum(converted)**2)/len(converted))/(len(converted)-1))
                        ))


class InvalidHeaderFormat(Exception):
    pass


class SeveralSequencePerSpeciesException(Exception):
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Filter marker genes based on how complete they cover the "
                                     "species set")
    parser.add_argument('-n', '--min-nr-species', type=int, required=True,
                        help="Specify the minimum number of species that should be "
                             "covered by each marker gene.")
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Path to a directory containing the input files to be filtered')
    parser.add_argument('-o', '--output', required=True,
                        help="Path to the output directory where the filtered fasta files will "
                             "be stored. This directory does not need to exist. If it does, existing "
                             "files will not be removed, but might be overwritten.")
    parser.add_argument('-e', '--ext', type=str, default=".fa",
                        help="Extension of files to be considered for conversion. (default %(default)s)")
    parser.add_argument('-#', '--threads', type=int,
                        help="Nr of threads to use to filter input files in parallel. "
                             "(defaults to the number of available cores)")
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help="Produce some output and status reports")
    parser.add_argument('-d', '--debug', action="store_true", default=False,
                        help="Be more verbose for debugging purposes")

    conf = parser.parse_args()
    level = logging.WARNING
    if conf.verbose:
        level = logging.INFO
    if conf.debug:
        level = logging.DEBUG
    logging.basicConfig(level=level, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logger.debug("Concatenate alignments: arguments: {}".format(conf))

    convert_files(conf.input, conf.output, conf.min_nr_species, pattern=conf.ext, nr_processes=conf.threads)
