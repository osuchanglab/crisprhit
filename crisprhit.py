#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import sys
import argparse
import os.path
import logging
# import warnings
from signal import signal, SIGPIPE, SIGINT, SIG_DFL
from Bio import SeqIO
from collections import defaultdict
signal(SIGPIPE, SIG_DFL)
signal(SIGINT, SIG_DFL)


__version__ = '0.9.0'
date = 'June 1, 2017'
positions = {}
positions['seed'] = range(1, 6) + [7, 8]
positions['interference'] = [6, 12, 18, 24, 30]
positions['priming'] = [28]
positions['stable'] = [10, 11, 12, 18, 22, 24, 25, 29, 30, 31]
positions['guide'] = 'SSSSSISS XXY     Y   X YX  *XYX'
# warnings.filterwarnings("error")


def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


class PrintVersion(argparse.Action):
    def __init__(self, nargs=0, **kwargs):
        if nargs != 0:
            raise ValueError('nargs for PrintVersion must be 0; '
                             'it is just a flag.')
        super(PrintVersion, self).__init__(nargs=nargs, **kwargs)

    def __call__(self, parser, values, namespace, option_string=None):
        fn = os.path.basename(__file__)
        sys.stderr.write('{} v. {}\n'.format(fn, __version__))
        sys.stderr.write('Generated {}\n'.format(date))
        parser.exit()


def init_logger(args):
    logger = logging.getLogger(__name__)
    ch = logging.StreamHandler()
    if args.debug:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)
        ch.setLevel(logging.INFO)
    elif args.quiet:
        logger.setLevel(logging.CRITICAL)
        ch.setLevel(logging.CRITICAL)
    else:
        logger.setLevel(logging.WARNING)
        ch.setLevel(logging.WARNING)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return(logger)


def run_argparse():
    parser = argparse.ArgumentParser(
            description='Extract sequence based on BLAST matching from a '
                        'FASTA or Genbank file.')
    parser.add_argument(
            'infile', help='Sequence file (FASTA or Genbank) of interest.',
            type=extant_file)
    parser.add_argument(
            'spacerfile', help='FASTA formatted spacers used as query.',
            type=extant_file)
    parser.add_argument(
            'listfiles', help='Tab-delimited BLAST data (outfmt 6 or 7) '
                              'containing output from search using CRISPR '
                              'spacer sequences as the query. "-" for STDIN.',
            type=str, nargs='+')
    parser.add_argument(
            '--filetype', help='File type provided (default = auto).',
            default='auto', choices=['FASTA', 'GENBANK', 'auto'])
    parser.add_argument(
            '--PAM', help='PAM sequence (default = TT)', default='TT',
            type=str)
    parser.add_argument(
            '--outtype', help='Output sequence type (default = all). Valid '
                              'for FASTA output only.',
            default='all', choices=['all', 'PAM'])
    parser.add_argument(
            '--outfmt', help='Output file format. (default = table)',
            default='table', choices=['fasta', 'table', 'basepair'])
    parser.add_argument(
            '--length', help='Matching protospacer lengths '
                             '(default = fill).',
            choices=['full', 'seed', 'fill'], default='fill')
    parser.add_argument(
            '--restrict', help='Restrict output length to this.', type=int,
            default=0)
    parser.add_argument(
            '--plength', help='PAM length to include. (default = 2)', type=int,
            default=2)
    parser.add_argument(
            '--width', help='Output column width. (default = 60)', type=int,
            default=60)
    parser.add_argument(
            '--match', help='Force PAM match. (default = False)',
            choices=['True', 'False', 'Partial'], default='False')
    parser.add_argument(
            '--filter', help='Filter output to a particular spacer category.',
            choices=['all', 'perfect', 'hq', 'priming', 'partial', 'other'],
            action='append', default=['all'])
    parser.add_argument(
            '--mmlimit', help='Set mismatch limit. (default = 15)', default=15,
            type=int)
    parser.add_argument(
            '--prefix', help='Prefix for arbitrary spacer names. '
                             '(default = ps)',
            type=str, default='ps')
    parser.add_argument(
            '--hits', help='Amount of hits to report for a spacer. '
                           '(default = all)',
            default='all', choices=['all', 'top'])
    parser.add_argument(
            '-v', '--verbose', help='Print progress messages.',
            action='store_true')
    parser.add_argument(
            '-q', '--quiet', help='Hide warning messages.',
            action='store_true')
    parser.add_argument(
            '--debug', help='Print debugging messages.', action='store_true')
    parser.add_argument(
            '-V', '--version', help='Print version message and quit.',
            action=PrintVersion)
    args = parser.parse_args()
    args.logger = init_logger(args)
    args.base = os.path.basename(args.infile)
    args.warning = 0
    if len(args.filter) > 1 and 'all' in args.filter:
        args.filter.remove('all')
    return args


def parse_tag(args, tag, tags):
    if tag.startswith('#') or not tag:
        return(args, tags)
    tag = tag.rstrip('\n')
    tagdata = tag.split('\t')
    if len(tagdata) < 12:
        msg = 'Cannot properly split BLAST data for {}. Skipping.'
        args.logger.warning(msg.format(tag))
        args.warning += 1
        return(args, tags)
    else:
        spacer = tagdata[0]
        target = tagdata[1]
        alnlen = int(tagdata[3])
        mismatch = int(tagdata[4])
        spacerstart = int(tagdata[6])
        spacerend = int(tagdata[7])
        start = int(tagdata[8])
        end = int(tagdata[9])
        data = (start, end, alnlen, mismatch, spacer, spacerstart, spacerend)
        tags[target].append(data)
    return(args, tags)


def parse_opts(args, tags):
    if args.listfiles[0] == '-':
        msg = 'Reading BLAST data from STDIN.'
        args.logger.debug(msg)
        for line in sys.stdin.read().split():
            args, tags = parse_tag(args, line, tags)
    else:
        msg = 'Reading coordinates from files:'
        args.logger.debug(msg)
        for listfile in args.listfiles:
            args.logger.debug(' + {}'.format(listfile))
            with open(listfile) as l:
                for line in l:
                    args, tags = parse_tag(args, line, tags)
    msg = 'Finished reading coordinates.'
    args.logger.debug(msg)
    if args.filetype == 'auto':
        msg = 'Automatically detecting file type.'
        args.logger.debug(msg)
        with open(args.infile, 'r') as f:
            line = f.readline().strip()
            if line[0] == '>':
                args.logger.debug('FASTA format detected.')
                args.filetype = 'FASTA'
            elif 'LOCUS' in line:
                args.logger.debug('Genbank format detected.')
                args.filetype = 'GENBANK'
            else:
                msg = 'Unable to determine filetype of {}.'
                args.logger.critical(msg.format(args.infile))
                msg = 'Provide filetype of sequence data with --filetype.'
                args.logger.critical(msg)
                sys.exit()
    if len(args.PAM) < args.plength:
        msg = 'Filling PAM with N to match given plength.'
        args.logger.info(msg)
        nfill = args.plength - len(args.PAM)
        args.PAM += 'N' * nfill
    elif len(args.PAM) > args.plength:
        msg = 'Setting plength to given PAM length.'
        args.logger.info(msg)
        args.plength = len(args.PAM)
    positions['guide'] = 'P' * args.plength + positions['guide']
    return(args, tags)


def parse_spacers(args):
    spacers = {}
    msg = 'Starting parse of spacer file {}.'
    args.logger.debug(msg.format(os.path.basename(args.spacerfile)))
    fasta = SeqIO.parse(args.spacerfile, 'fasta')
    for record in fasta:
        spacers[record.id] = (record.seq, record.reverse_complement().seq)
    return spacers


def compare_spacer(args, spacer, pspacer, PAMhits):
    hit = 0
    seedhit = 0
    matchstick = ''
    # mm = mismatch counting
    # hq, seed, prime, stable, total
    mm = defaultdict(int)
    hit += PAMhits
    mm['pam'] = len(args.PAM) - PAMhits
    mm['total'] += mm['pam']
    # Need to reverse sequences for basepairing
    rspacer = spacer[::-1]
    rpspacer = pspacer[::-1]
    for i, c in enumerate(rspacer):
        pos = i + 1
        try:
            if c == str(rpspacer[i]):
                matchstick += '|'
                hit += 1
                if pos in positions['seed']:
                    seedhit += 1
            else:
                matchstick += ' '
                mm['total'] += 1
                if pos in positions['interference']:
                    mm['hq'] += 1
                if pos in positions['seed']:
                    mm['seed'] += 1
                if pos in positions['priming']:
                    mm['priming'] += 1
                if pos in positions['stable']:
                    mm['stable'] += 1
        except IndexError:
            msg = 'Length of sequence ({}) does not match spacer ({})'
            args.logger.critical(msg.format(len(pspacer), len(spacer)))
            sys.exit()
    quality = 'other'
    if not mm['total']:
        quality = 'perfect'
    elif mm['hq'] <= 2 and mm['seed'] == 0 and mm['pam'] == 0 and \
            mm['priming'] == 0 and mm['total'] <= 3:
        quality = 'hq'
    elif mm['stable'] <= 1 and mm['total'] <= 8 and mm['priming'] == 1:
        quality = 'priming'
    elif mm['stable'] <= 2 and mm['pam'] <= 1 and mm['total'] <= 7:
        quality = 'priming'
    elif mm['stable'] >= 3:
        quality = 'stable'
    elif mm['seed'] <= 1 and mm['pam'] <= 1:
        quality = 'partial'
    elif mm['seed'] <= 2 and mm['pam'] == 0:
        quality = 'partial'
    elif mm['seed'] == 0:
        quality = 'partial'
    return(hit, mm['total'], seedhit, mm['priming'], quality, matchstick)


def parse_contigs(args, tags, spacers):
    msg = 'Starting parse of file {}.'
    args.logger.debug(msg.format(args.base))
    infile = SeqIO.parse(args.infile, args.filetype.lower())
    output = []
    pamcheck = defaultdict(lambda: defaultdict(int))
    seen = set()
    # ps = protospacer count
    ps = 1
    prefix = 'ps'
    for record in infile:
        if record.id in tags:
            k = 1
            for line in tags[record.id]:
                # data = (start, end, alnlen,
                #         mismatch, target, spacerstart, spacerend)
                # p for proto, need reverse strand from hit, so we swap them
                pstart = line[1]
                pend = line[0]
                skip = False
                spacername = line[4]
                if args.hits == 'top' and spacername in seen:
                    msg = 'Skipping {} hit {} by {} due to multiple hits.'
                    args.logger.debug(msg.format(record.id, k, spacername))
                    k += 1
                    continue
                else:
                    seen.add(spacername)
                slen = len(spacers[spacername][0])

                if args.length == 'full':
                    if line[2] != slen:
                        skip = True
                elif args.length == 'seed':
                    if line[2] != line[6]:
                        skip = True
                elif args.length == 'fill':
                    lfill = slen - line[6]
                    rfill = line[5] - 1
                    if pstart < pend:
                        # pstart smaller, subtract from start, add to end
                        pstart = pstart - lfill
                        pend = pend + rfill
                    else:
                        # pstart larger, add to start, subtract from end
                        pstart = pstart + lfill
                        pend = pend - rfill
                if skip:
                    msg = 'Skipping {} hit {} due to length mismatch.'
                    args.logger.debug(msg.format(record.id, k))
                    k += 1
                    continue
                PAM = ''
                seq = ''
                strand = ''
                # plength is PAM length
                plength = args.plength
                if pstart < pend:
                    strand = 'plus'
                    start = pstart
                    end = pend
                    seq = record[start - 1:end].seq
                    PAM = record[end:end + plength].seq.lower()
                else:
                    strand = 'minus'
                    start = pend
                    end = pstart
                    plength += 1
                    seq = record[start - 1:end].reverse_complement().seq
                    PAM = record[start - plength:start - 1]\
                        .reverse_complement().seq.lower()
                skip = False
                PAMhits = 0
                for i, base in enumerate(args.PAM):
                    if base.lower() == PAM[i] or PAM[i] == 'N':
                        PAMhits += 1
                if args.match == 'Partial' and not PAMhits:
                    skip = True
                elif args.match == 'True' and PAMhits != len(args.PAM):
                    skip = True
                if skip:
                    msg = 'Skipping {} hit {} due to PAM mismatch.'
                    args.logger.debug(msg.format(record.id, k))
                    k += 1
                    continue
                outline = {}
                # Must compare to reverse complement of spacer
                spacer = spacers[spacername][1]
                # Fill sequence for seed matches
                if slen != len(seq):
                    fill = 'N' * (slen - len(seq))
                    seq = fill + seq
                hits, mismatches, seedhits, primingmm, qual, matchstick = \
                    compare_spacer(args, spacer, seq, PAMhits)
                if mismatches > args.mmlimit:
                    msg = 'Skipping {} hit {} due to hitting mismatch limit.'
                    args.logger.debug(msg.format(record.id, k))
                    k += 1
                    continue
                outline['id'] = record.id
                outline['start'] = str(pstart)
                outline['end'] = str(pend)
                outline['spacer'] = spacername
                outline['seq'] = seq
                outline['PAM'] = PAM
                outline['PAMhits'] = PAMhits
                outline['strand'] = strand
                outline['length'] = len(seq)
                outline['hits'] = hits
                outline['mismatches'] = mismatches
                outline['seedhits'] = seedhits
                outline['primingmm'] = primingmm
                outline['matchstick'] = matchstick
                outline['quality'] = qual
                outline['name'] = prefix + str(ps).zfill(5)
                PAMn = list(str(PAM))
                if args.PAM.find('N') > -1:
                    for i in [p for p, c in enumerate(args.PAM) if c == 'N']:
                        PAMn[i] = 'n'
                PAMn = ''.join(PAMn)
                pamcheck[PAMn]['matches'] += 1
                pamcheck[PAMn]['hits'] += hits
                pamcheck['total']['matches'] += 1
                pamcheck['total']['hits'] += hits
                output.append(outline)
                k += 1
                ps += 1
    msg = 'Finished searching file {}.'
    args.logger.debug(msg.format(args.base))
    return(output, pamcheck)


def print_output(args, output, spacers):
    msg = 'Printing output of file {}.'
    args.logger.debug(msg.format(args.base))
    if args.outfmt == 'table':
        header = ['#name', 'id', 'seq', 'PAM', 'start', 'end', 'strand',
                  'spacer', 'hits', 'misses', 'PAM_hits', 'seed_hits',
                  'priming_mm', 'quality']
        print('\t'.join(header))
    for line in output:
        if 'all' in args.filter:
            pass
        elif line['quality'] not in args.filter:
            msg = 'Skipping {} due to quality filter'
            args.logger.info(msg.format(line['name']))
            continue
        if args.outfmt == 'fasta':
            out = ''
            if line['strand'] == 'plus':
                out = '>' + line['id'] + ':' + line['start'] + '..' + \
                       line['end'] + '_' + line['spacer']
            elif line['strand'] == 'minus':
                out = '>' + line['id'] + ':c' + line['end'] + '..' + \
                      line['start'] + '_' + line['spacer']
            out += ' ' + 'hits:{} '.format(line['hits']) + \
                   'mismatch:{} '.format(line['mismatches']) + \
                   'PAM_hits:{}'.format(line['PAMhits']) + \
                   'seed_hits:{}'.format(line['seedhits']) + \
                   'priming_mm:{}'.format(line['primingmm']) + \
                   'quality:{}'.format(line['quality'])
            print(out)
            if args.outtype == 'all':
                out = line['seq'] + line['PAM']
            elif args.outtype == 'PAM':
                out = line['PAM']
            if args.restrict:
                print(out[-args.restrict:])
            else:
                print(out)
        elif args.outfmt == 'table':
            out = [line['name'], line['id'], line['seq'], line['PAM'],
                   line['start'], line['end'], line['strand'],
                   line['spacer'], line['hits'], line['mismatches'],
                   line['PAMhits'], line['seedhits'], line['primingmm'],
                   line['quality']]
            print('\t'.join(map(str, out)))
        elif args.outfmt == 'basepair':
            # ps = protospacer
            # s = spacer
            lfill = 10
            slen = len(spacers[line['spacer']][0])
            out = [line['name'], line['id'], line['start'] + ':' + line['end'],
                   line['strand'], line['spacer'], line['quality']]
            print('#' + ' '.join(map(str, out)))
            out = positions['guide']
            print(' ' * lfill + out)
            print(' ' * (lfill + args.plength) + spacers[line['spacer']][0])
            seq = line['seq']
            out = line['matchstick']
            print(' ' * (lfill + args.plength) + out)
            rpad = args.width - slen - args.plength - lfill
            print(line['start'].ljust(lfill) + line['PAM'][::-1] +
                  seq[::-1] + line['end'].rjust(rpad))


def compare_PAM(args, pamcheck):
    givenpam = args.PAM.lower()
    PAMscores = []
    top10 = []
    for pamseq in pamcheck:
        score = (pamcheck[pamseq]['hits']/pamcheck[pamseq]['matches']) * \
                (pamcheck[pamseq]['matches']/pamcheck['total']['matches'])
        PAMscores.append((pamseq, score))
    top10 = sorted(PAMscores, key=lambda x: x[1], reverse=True)[1:11]
    matched = 0
    for i, v in enumerate(top10):
        args.logger.debug('gPAM={}\tePAM={}'.format(givenpam, v[0]))
        if givenpam == v[0]:
            matched = i + 1
            continue
    if matched:
        if matched == 1:
            msg = 'Given PAM ({}) is consistent with empirical data.'
            args.logger.info(msg.format(givenpam.upper()))
            return True
        else:
            msg = 'Given PAM ({}) is in the top 10 PAM sequences ({}).'
            args.logger.info(msg.format(givenpam.upper(), matched))
    msg = 'Given PAM sequence might not be accurate.'
    args.logger.warning(msg)
    if not args.verbose and not args.debug:
        msg = 'Use --verbose to see possible PAM sequences.'
        args.logger.warning(msg)
    msg = 'These are possible PAM sequences for your data:\n'
    msg += ' ' * 10 + '\t'.join(['Rank', 'Seq', 'Score', 'Count'])
    for i, v in enumerate(top10):
        pamhits = str(pamcheck[v[0]]['matches'])
        item = '\n'
        item += ' ' * 10
        item += '{}\t{}\t{:.2f}\t{}'.format(i+1, v[0].upper(), v[1], pamhits)
        msg += item
    args.logger.info(msg)
    return False


def main():
    args = run_argparse()

    tags = defaultdict(list)

    args, tags = parse_opts(args, tags)

    spacers = parse_spacers(args)

    output, pamcheck = parse_contigs(args, tags, spacers)

    print_output(args, output, spacers)

    compare_PAM(args, pamcheck)


if __name__ == '__main__':
    main()
