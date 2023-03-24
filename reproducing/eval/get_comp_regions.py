#!/usr/bin/env python3

import argparse
import gzip
import csv
from parascopy.inner import common
from parascopy.inner.genome import Genome, Interval


def create_reader(inp):
    for line in inp:
        if line.startswith('##'):
            continue
        assert line.startswith('#')
        fields = line.lstrip('#').strip().split('\t')
        return csv.DictReader(inp, delimiter='\t', fieldnames=fields)


def load_chrom_cn(inp, genome):
    chrom_cn = {}
    for line in inp:
        if line.startswith('#'):
            continue
        chrom, start, end, samples, cn = line.strip().split('\t')
        assert start == '0' and end == 'inf'
        chrom_id = genome.chrom_id(chrom)
        for sample in samples.split(','):
            chrom_cn[(sample, chrom_id)] = int(cn)
    return chrom_cn


def get_ref_agcn(regions, sample, genome, chrom_cn):
    total_cn = 0
    for region in regions:
        chrom_id = region.chrom_id
        chrom_cn.get((sample, chrom_id))
        cn = chrom_cn.get((sample, chrom_id))
        if cn is None:
            cn = chrom_cn.get(('*', chrom_id))
        if cn is None:
            cn = 2
        total_cn += cn
    return total_cn


def get_regions(row, genome):
    regions = [Interval(genome.chrom_id(row['chrom']), int(row['start']), int(row['end']))]
    if row['homologous_regions'] == '*':
        return regions
    for s in row['homologous_regions'].split(','):
        regions.append(Interval.parse(s, genome))
    return regions


def process(reader, genome, chrom_cn, args):
    qual_thresh = args.qual
    filters_cond = args.filters
    cn_cond = args.cn
    min_refcn, max_refcn = args.cn_bounds
    pooled = args.pooled

    for row in reader:
        regions = get_regions(row, genome)
        region1 = regions[0]
        sample = row['sample']
        ref_agcn = get_ref_agcn(regions, sample, genome, chrom_cn)

        if not (min_refcn <= ref_agcn <= max_refcn):
            continue

        if pooled:
            min_qual = float(row['qual'])
            filters_pass = row['filter'] == 'PASS'
            agcn_match = ref_agcn == int(row['CN'])
            pscn_match = True

        else:
            min_qual = min(float(row['qual']), float(row['agCN_qual']))
            filters_pass = row['filter'] == 'PASS' and row['agCN_filter'] == 'PASS'
            ref_pscn = get_ref_agcn((region1,), sample, genome, chrom_cn)
            agcn_match = ref_agcn == int(row['agCN'])
            pscn_match = ref_pscn == int(row['CN'])

        if (qual_thresh >= 0 and min_qual < qual_thresh) or (qual_thresh < 0 and min_qual >= -qual_thresh):
            continue
        if (filters_cond == 'only-pass' and not filters_pass) or (filters_cond == 'not-pass' and filters_pass):
            continue
        cn_match = agcn_match and pscn_match
        if (cn_cond == 'only-ref' and not cn_match) or (cn_cond == 'non-ref' and cn_match):
            continue
        if cn_cond == 'pscn-ref' and not pscn_match:
            continue

        yield region1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, metavar='<file>',
        help='Input bed.gz file.')
    parser.add_argument('-f', '--fasta-ref', required=True, metavar='<file>',
        help='Fasta reference file.')
    parser.add_argument('-o', '--output', required=True, metavar='<file>',
        help='Output BED file.')
    parser.add_argument('-r', '--ref-cn', required=True, metavar='<file>',
        help='Input BED file with reference copy numbers (see --assume-cn in "parascopy call").')

    parser.add_argument('-q', '--qual', default=20, type=float, metavar='<float>',
        help='Copy number quality threshold [default: %(default)s]. '
            'If negative, only output regions where aggregate or paralog-specific CN quality '
            'is LESS than the absolute threshold value.')
    parser.add_argument('-F', '--filters', choices=('only-pass', 'any', 'not-pass'), default='only-pass',
        help='Output regions with the corresponding filters [default: %(default)s].')
    parser.add_argument('-n', '--cn', choices=('only-ref', 'any', 'non-ref', 'pscn-ref'), default='only-ref',
        help='Output regions with the matching (non-matching) copy number and '
            'the reference copy number [default: %(default)s].')
    parser.add_argument('-c', '--cn-bounds', default=(3, 8), nargs=2, type=int, metavar='<int>',
        help='Minimal and maximal reference aggr. copy number values [default: 3 8].')
    parser.add_argument('-p', '--pooled', action='store_true',
        help='Create pooled calling regions.')
    args = parser.parse_args()

    with gzip.open(args.input, 'rt') as inp, open(args.output, 'w') as out, Genome(args.fasta_ref) as genome:
        with common.open_possible_gzip(args.ref_cn) as ref_cn_inp:
            chrom_cn = load_chrom_cn(ref_cn_inp, genome)

        reader = create_reader(inp)
        regions = list(process(reader, genome, chrom_cn, args))
        regions.sort()
        regions = Interval.combine_overlapping(regions)

        out.write('# {}\n'.format(common.command_to_str()))
        for region in regions:
            out.write(region.to_bed(genome) + '\n')


if __name__ == '__main__':
    main()
