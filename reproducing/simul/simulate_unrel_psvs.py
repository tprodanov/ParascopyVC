#!/usr/bin/env python3

import argparse
import pysam
import numpy as np
import random
import sys

import parascopy.inner.common as common


def create_header(in_header, sample):
    out_header = in_header.copy()
    out_header.add_line('##command="{}"'.format(' '.join(sys.argv)))
    out_header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    out_header.add_sample(sample)
    return out_header


def simulate_psvs(in_vcf, out_vcf, rate, max_len):
    last_chrom = None
    last_end = None

    for psv in in_vcf:
        if psv.chrom == last_chrom and psv.start <= last_end:
            continue
        if max(map(len, psv.alleles)) > max_len:
            continue
        fval = psv.info['fval']
        if np.isnan(fval):
            continue

        q = np.clip((1 - fval) * rate, 0, 1)
        p = 1 - q
        ref_hom = p * p
        alt_hom = q * q
        r = random.random()
        alt_allele = random.randrange(1, len(psv.alleles))
        if r <= ref_hom:
            continue
        elif r <= ref_hom + alt_hom:
            gt = (alt_allele, alt_allele)
        else:
            gt = (alt_allele, 0) if random.randrange(2) else (0, alt_allele)

        last_chrom = psv.chrom
        last_end = psv.start + len(psv.ref)
        var = out_vcf.new_record()
        var.chrom = psv.chrom
        var.start = psv.start
        var.alleles = psv.alleles
        for key in psv.info:
            var.info[key] = psv.info[key]

        var.samples[0]['GT'] = gt
        var.samples[0].phased = True
        out_vcf.write(var)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True,
        help='Input VCF file with expanded PSVs.')
    parser.add_argument('-r', '--rate', type=float, default=1,
        help='Unreliable PSV rate. When rate = 1, generates PSVs with probabilities '
            'according to f-values. Smaller rate leads to less unreliable PSVs, and higher '
            'rate leadds to more unreliable PSVs. [default: %(default)s].')
    parser.add_argument('-o', '--output', required=True,
        help='Output VCF file.')

    parser.add_argument('-s', '--sample', default='sim',
        help='Sample name [default: "%(default)s"].')
    parser.add_argument('-S', '--seed', type=int,
        help='Optional seed.')
    parser.add_argument('-m', '--max-len', type=int, default=30,
        help='Maximum PSV length [default: %(default)s].')
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    with pysam.VariantFile(args.input) as in_vcf:
        header = create_header(in_vcf.header, args.sample)
        with pysam.VariantFile(args.output, 'wz' if args.output.endswith('.gz') else 'w', header=header) as out_vcf:
            simulate_psvs(in_vcf, out_vcf, args.rate, args.max_len)

    if args.output.endswith('.gz'):
        common.Process(['tabix', '-p', 'vcf', args.output]).finish()


if __name__ == '__main__':
    main()
