#!/usr/bin/env python3
"""Dump a .hic intra-chromosomal map to a 3-column sparse matrix for SpectralTAD."""

from __future__ import annotations

import argparse
from pathlib import Path

import hicstraw


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--hic", required=True)
    p.add_argument("--chrom", required=True, help="Chromosome name as stored in the .hic (e.g. 22)")
    p.add_argument("--binsize", type=int, default=10000)
    p.add_argument("--norm", default="NONE", help="NONE/VC/VC_SQRT/KR (SpectralTAD recommends NONE)")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    hic = hicstraw.HiCFile(args.hic)
    chrom_names = {c.name for c in hic.getChromosomes()}
    if args.chrom not in chrom_names:
        raise SystemExit(f"Chromosome {args.chrom!r} not in file. Available: {sorted(chrom_names)}")

    records = hicstraw.straw(
        "observed", args.norm, args.hic, args.chrom, args.chrom, "BP", int(args.binsize)
    )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with out.open("w") as fh:
        for rec in records:
            if hasattr(rec, "binX"):
                x, y, v = int(rec.binX), int(rec.binY), float(rec.counts)
            else:
                x, y, v = rec
            if not (v == v) or v == 0:
                continue
            fh.write(f"{int(x)}\t{int(y)}\t{v}\n")
            n += 1
    print(f"Wrote {n} contacts to {out}")


if __name__ == "__main__":
    main()
