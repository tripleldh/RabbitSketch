#!/usr/bin/env python3
"""
subsample_fastq.py – randomly subsample a FASTQ/FASTQ.gz file.
Usage:  python subsample_fastq.py <input.fastq[.gz]> <fraction> <out.fastq.gz> [seed=42]
"""
import sys, gzip, random

def open_any(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)

def subsample(inp, frac, out, seed=42):
    rng = random.Random(seed)
    written = total = 0
    with open_any(inp) as fi, gzip.open(out, 'wt') as fo:
        while True:
            lines = [fi.readline() for _ in range(4)]
            if not lines[0]:
                break
            total += 1
            if rng.random() < frac:
                fo.writelines(lines)
                written += 1
    print(f"  {inp}: kept {written}/{total} reads ({written/total*100:.1f}%)", flush=True)

if __name__ == '__main__':
    if len(sys.argv) < 4:
        sys.exit("Usage: subsample_fastq.py <input> <fraction> <output> [seed]")
    inp   = sys.argv[1]
    frac  = float(sys.argv[2])
    out   = sys.argv[3]
    seed  = int(sys.argv[4]) if len(sys.argv) > 4 else 42
    subsample(inp, frac, out, seed)
