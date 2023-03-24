#!/usr/bin/env python3

import sys
import gzip
import os
from collections import defaultdict
import numpy as np
import argparse


def float_or_nan(s):
    if s == 'None':
        return np.nan
    return float(s)


def write_summary(filename, thresholds, only_thresholds):
    if os.path.isfile(filename):
        dirname = os.path.dirname(filename)
    else:
        dirname = filename
        filename = os.path.join(dirname, 'weighted_roc.tsv.gz')

    print(f'====  {dirname}  ====')
    if filename.endswith('.gz'):
        inp = gzip.open(filename, 'rt')
    else:
        inp = open(filename)

    n_baseline = None
    n_calls = None
    roc_matrix = []
    for line in inp:
        line = line.strip().split()
        if line[0].startswith('#'):
            if line[1] == 'baseline':
                n_baseline = int(line[-1])
                print(f'Total baseline:      {n_baseline}')
            elif line[1] == 'call':
                n_calls = int(line[-1])
                print(f'Total call variants: {n_calls}')
            elif line[1] == 'true_positives_baseline':
                print('Threshold    TP(b)    TP(c)       FP       FN  Precision  Recall  F1-score')
                print('--------------------------------------------------------------------------')
        else:
            roc_matrix.append(list(map(float_or_nan, line)))
    inp.close()

    if not n_baseline or not n_calls:
        print('Cannot continue\n')
        return

    roc_matrix = np.array(roc_matrix)
    # Swap two columns.
    roc_matrix[:, [2, 3]] = roc_matrix[:, [3, 2]]

    ixs = defaultdict(list)
    for thresh in thresholds:
        curr_ixs = np.where(roc_matrix[:, 0] >= thresh)[0]
        if len(curr_ixs) > 0:
            curr_i = curr_ixs[-1]
            ixs[curr_i].append('≥ {:.0f}'.format(thresh))

    true_pos_rate = roc_matrix[:, 1] / n_baseline
    if not only_thresholds:
        for col, flag in [(5, 'best precision'), (6, 'best recall'), (7, 'best F1')]:
            curr_i = np.argmax(roc_matrix[:, col] + true_pos_rate * 0.0001)
            ixs[curr_i].append(flag)

    for i, flags in sorted(ixs.items()):
        s = '{:9.3f}  {:7.0f}  {:7.0f}  {:7.0f}  {:7.0f}  {:9.4f}  {:6.4f}  {:8.4f}'.format(*roc_matrix[i])
        if flags:
            s += '  ({})'.format(', '.join(flags))
        print(s)
    print()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', nargs='+', help='Input directories')
    parser.add_argument('-t', '--thresholds', nargs='+', default=(0, 10, 20, 50, 100), type=float,
        help='Quality thresholds [default %(default)s].')
    parser.add_argument('-T', '--only-thresholds', action='store_true',
        help='Print only accuracy values for quality thresholds.')
    args = parser.parse_args()

    for path in args.input:
        write_summary(path, args.thresholds, args.only_thresholds)


if __name__ == '__main__':
    main()
