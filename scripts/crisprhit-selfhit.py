#!/usr/bin/env python
from __future__ import print_function
import sys
import argparse
import csv
import os.path
import logging
#import operator
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
    data = {}  # dict with data about spacer
    positions = {}  # k = spacer_revcom; v = tuple of start, end pos
    red_spacers = defaultdict(list)
    with open(args.infile, 'rb') as tsv:
        reader = csv.DictReader(tsv, delimiter='\t')
        for line in reader:
            # print(line)
            s_revcom = line['spacer_revcom']
            red_spacers[s_revcom].append(line['spacer'])
            if s_revcom not in data:
                data[s_revcom] = line
                positions[s_revcom] = (int(line['start']), int(line['end']))
    return(data, positions, red_spacers, reader.fieldnames)


def find_overlap(args, positions):
    overlaps = []
    grouplist = []
    pos_sorted = sorted(positions.values())
    refcoord = pos_sorted[0]
    pos_sorted.remove(refcoord)
    grouplist.append(refcoord)
    for i, coord in enumerate(pos_sorted):
        r = range(min(refcoord), max(refcoord))
        if coord[0] in r or coord[1] in r:
            if coord not in grouplist:
                grouplist.append(coord)
        else:
            overlaps.append(grouplist)
            grouplist = [coord]
        if i == len(pos_sorted) - 1:
            overlaps.append(grouplist)
        refcoord = coord
    ol_check = {}
    for i, ol in enumerate(overlaps):
        overlap = -1
        if len(ol) > 1:
            overlap = i
        for coord in ol:
            ol_check[coord] = overlap
    return(ol_check)


def print_overlap(args, data, ol_check, red_spacers, header):
    output = []
    for spacer in data:
        coord = (int(data[spacer]['start']), int(data[spacer]['end']))
        if ol_check[coord] > -1:
            line = [ol_check[coord]] + \
                    [data[spacer][x] for x in header]
            output.append(line)
    header = ['group'] + header
    output.sort()
    for line in output:
        line[0] = 'group{}'.format(line[0])
    writer = csv.writer(sys.stdout, delimiter='\t')
    writer.writerow(header)
    writer.writerows(output)


def main():
    args = run_argparse()
    data, positions, red_spacers, header = parse_crisprhit(args)
    ol_check = find_overlap(args, positions)
    print_overlap(args, data, ol_check, red_spacers, header)


if __name__ == '__main__':
    main()
