#!/usr/bin/env python
from __future__ import print_function
import sys
import argparse
import csv
import os.path
import logging
from collections import defaultdict
from signal import signal, SIGPIPE, SIGINT, SIG_DFL
signal(SIGPIPE, SIG_DFL)
signal(SIGINT, SIG_DFL)


def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def init_logger(args):
    logger = logging.getLogger(__name__)
    ch = logging.StreamHandler()
    if args.debug:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)
        ch.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
        ch.setLevel(logging.WARNING)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return(logger)


def run_argparse():
    parser = argparse.ArgumentParser(
        description='Extract spacers that target the same protospacer'
                    ' sequence.')
    parser.add_argument(
        'infile', help='Input crisprhit file from crisprhit.py.',
        type=extant_file)
    parser.add_argument(
        '--debug', help='Turn on debugging messages.', action='store_true')
    parser.add_argument(
        '--verbose', help='Print verbose progress messages.',
        action='store_true')
    args = parser.parse_args()
    args.logger = init_logger(args)
    return args


def parse_crisprhit(args):
    spacers_revcom = []  # sequence of spacer_revcom for each entry
    data = {}  # dict with data about spacer
    positions = {}  # k = spacer_revcom; v = tuple of start, end pos
    red_spacers = defaultdict(list)
    with open(args.infile, 'rb') as tsv:
        reader = csv.DictReader(tsv, delimiter='\t')
        for i, line in enumerate(reader):
            s_revcom = line['spacer_revcom']
            red_spacers[s_revcom].append(line['spacer'])
            if s_revcom not in spacers_revcom:
                spacers_revcom.append(s_revcom)
                data[s_revcom] = line
                positions[s_revcom] = (line['start'], line['end'])
    return(spacers_revcom, data, positions, red_spacers)


def find_overlap(args, spacers_revcom, positions):
    for i, seq in enumerate(spacers_revcom):
        pos1 = positions[seq]


def main():
    args = run_argparse()
    spacers_revcom, data, positions, red_spacers = parse_crisprhit(args)
    find_overlap(args, spacers_revcom, positions)


if __name__ == '__main__':
    main()
