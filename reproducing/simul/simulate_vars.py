#!/usr/bin/env python3

import argparse
import pysam
import random
import sys
from tqdm import tqdm
import scipy.stats
import numpy as np
from numpy.random import MT19937, RandomState, SeedSequence

import parascopy.inner.common as common
import parascopy.inner.itree as itree
from parascopy.inner.genome import Genome, Interval


def get_exclude_variants(args, genome):
    if args.exclude is None:
        return None

    excl_tree = itree.MultiChromTree()
    with pysam.VariantFile(args.exclude) as inp:
        for record in inp:
            rec_interval = Interval(genome.chrom_id(record.chrom), record.start, record.start + len(record.ref))
            excl_tree.add(rec_interval, None)
    return excl_tree


def create_header(genome, sample, chroms):
    vcf_header = pysam.VariantHeader()
    vcf_header.add_line('##command="{}"'.format(' '.join(sys.argv)))

    for name, length in zip(genome.chrom_names, genome.chrom_lengths):
        if chroms is None or name in chroms:
            vcf_header.add_line('##contig=<ID={},length={}>'.format(name, length))
    vcf_header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    vcf_header.add_sample(sample)
    return vcf_header


def create_substitution(ref, rs):
    alt = ref
    while alt == ref:
        alt = 'ACGT'[rs.randint(4)]
    return alt


def create_insertion(ref, rs):
    insert_length = rs.randint(1, 4)
    alt = ref + ''.join('ACGT'[rs.randint(4)] for _ in range(insert_length))
    return alt


def create_deletion(seq, i, rs):
    deletion_length = rs.randint(1, 4)
    ref = seq[i : i + 1 + deletion_length]
    alt = seq[i]
    assert len(ref) > 1
    return ref, alt


def generate_chrom(genome, chrom_id, exclude_tree, mut_rates, hetero_rate, out_vcf, rs):
    chrom = genome.chrom_name(chrom_id)
    chrom_len = genome.chrom_len(chrom_id)
    seq = genome.fetch_interval(Interval(chrom_id, 0, chrom_len))
    mut_rates = np.cumsum(mut_rates)
    mut_prob = mut_rates[-1]
    assert mut_prob <= 1.0

    with tqdm(total=chrom_len) as pbar:
        i = 0
        while True:
            dist = scipy.stats.geom.rvs(mut_prob, random_state=rs)
            i += dist
            if i >= chrom_len:
                break
            pbar.update(dist)

            ref = seq[i]
            if ref == 'N':
                i += 1
                pbar.update()
                continue

            r = rs.random_sample() * mut_prob
            if r <= mut_rates[0]:
                alt = create_substitution(ref, rs)
            elif r <= mut_rates[1]:
                alt = create_insertion(ref, rs)
            else:
                ref, alt = create_deletion(seq, i, rs)
                if 'N' in ref:
                    i += 1
                    pbar.update()
                    continue

            if exclude_tree is not None:
                region = Interval(chrom_id, i - 1, i + len(ref) + 1)
                if exclude_tree.overlap_size(region) > 0:
                    i += 1
                    pbar.update()
                    continue

            r2 = rs.random_sample()
            if r2 > hetero_rate:
                gt = (1, 1)
            elif r2 > hetero_rate / 2:
                gt = (0, 1)
            else:
                gt = (1, 0)

            record = out_vcf.new_record()
            record.chrom = chrom
            record.start = i
            record.alleles = (ref, alt)
            record.samples[0]['GT'] = gt
            record.samples[0].phased = True
            out_vcf.write(record)

            i_upd = len(ref) + 1
            i += i_upd
            pbar.update(i_upd)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta-ref', required=True,
        help='Fasta reference file.')
    parser.add_argument('-m', '--mutation-rate', required=True, nargs=3, type=float,
        help='Mutation rate (substitution, insertion, deletion).')
    parser.add_argument('-H', '--hetero-rate', required=True, type=float,
        help='Rate of the heterozygous variants.')
    parser.add_argument('-o', '--output', required=True,
        help='Output VCF file.')

    parser.add_argument('-e', '--exclude', required=False,
        help='Optional: VCF file, do not generate variants that overlap exclude variants.')
    parser.add_argument('-s', '--sample', default='sim',
        help='Sample name [default: "%(default)s"].')
    parser.add_argument('-S', '--seed', type=int,
        help='Optional seed.')
    parser.add_argument('-c', '--chroms', nargs='+',
        help='Optional: Only generate for a set of chromosoms.')
    args = parser.parse_args()

    genome = Genome(args.fasta_ref)
    excl_tree = get_exclude_variants(args, genome)
    chroms = set(args.chroms) if args.chroms else None
    header = create_header(genome, args.sample, chroms)

    rs = RandomState(MT19937(SeedSequence(args.seed)))
    with pysam.VariantFile(args.output, 'wz' if args.output.endswith('.gz') else 'w', header=header) as out_vcf:
        for chrom_id, chrom in enumerate(genome.chrom_names):
            if chroms is None or chrom in chroms:
                sys.stderr.write('Chromosome {}\n'.format(chrom))
                generate_chrom(genome, chrom_id, excl_tree, args.mutation_rate, args.hetero_rate, out_vcf, rs)

    if args.output.endswith('.gz'):
        common.Process(['tabix', '-p', 'vcf', args.output]).finish()


if __name__ == '__main__':
    main()
