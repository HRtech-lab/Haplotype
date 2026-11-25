"""Read-backed phasing of diploid, biallelic single-nucleotide variants.

Use --help for the BAM/VCF interface. Positions in reports are one-based.
SAM and VCF specifications: https://samtools.github.io/hts-specs/
Read-based phasing: https://whatshap.readthedocs.io/en/latest/guide.html
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from contextlib import contextmanager
import os
from pathlib import Path
from tempfile import NamedTemporaryFile
import pysam


def load_sites(path, sample=None):
    sites = {}
    with pysam.VariantFile(str(path)) as source:
        names = list(source.header.samples)
        if sample is None and len(names) == 1:
            sample = names[0]
        if sample not in names:
            raise ValueError(
                "Select one VCF sample with --sample. Available: " + ", ".join(names)
            )
        for r in source:
            if set(r.filter.keys()) - {"PASS", "."}:
                continue
            alleles = r.alleles or ()
            if (
                len(alleles) != 2
                or any(base not in {"A", "C", "G", "T"} for base in alleles)
                or alleles[0] == alleles[1]
            ):
                continue
            if r.samples[sample].get("GT") not in ((0, 1), (1, 0)):
                continue
            pos = sites.setdefault(r.contig, {})
            if r.pos in pos and pos[r.pos] != alleles:
                raise ValueError(f"Inconsistent variant records at {r.contig}:{r.pos}")
            pos[r.pos] = alleles
    return sample, sites


def fragment_patterns(path, sites, sample, min_mapq=20, min_baseq=20):
    """Merge observed alleles once per fragment. Sets retain conflicts between mates."""
    fragments = defaultdict(lambda: defaultdict(set))
    mode = "r" if str(path).lower().endswith(".sam") else "rb"
    with pysam.AlignmentFile(str(path), mode) as reads:
        rg_to_sm = {d["ID"]: d.get("SM") for d in reads.header.to_dict().get("RG", [])}
        known = {v for v in rg_to_sm.values() if v}
        if known and sample not in known:
            raise ValueError(f"Sample {sample!r} is absent from BAM read-group samples")
        if sites and not set(sites).intersection(reads.references):
            raise ValueError("VCF and alignments have no shared contig names")
        for r in reads.fetch(until_eof=True):
            if (
                r.is_unmapped
                or r.is_secondary
                or r.is_supplementary
                or r.is_duplicate
                or r.is_qcfail
            ):
                continue
            if r.mapping_quality == 255 or r.mapping_quality < min_mapq:
                continue
            if r.is_paired and (
                not r.is_proper_pair
                or r.mate_is_unmapped
                or r.next_reference_id != r.reference_id
            ):
                continue
            chrom = r.reference_name
            variants = sites.get(chrom, {})
            if not variants or not r.query_sequence or not r.query_name:
                continue
            group = r.get_tag("RG") if r.has_tag("RG") else ""
            if known and rg_to_sm.get(group) != sample:
                continue
            qualities = r.query_qualities
            if qualities is None and min_baseq > 0:
                continue
            key = (chrom, group, r.query_name)
            for qpos, rpos in r.get_aligned_pairs(matches_only=True):
                pos = rpos + 1
                if pos not in variants or (
                    qualities is not None and qualities[qpos] < min_baseq
                ):
                    continue
                base = r.query_sequence[qpos]
                if base in variants[pos]:
                    fragments[key][pos].add(variants[pos].index(base))
    counts = defaultdict(Counter)
    for (chrom, _, _), calls in fragments.items():
        pattern = tuple(
            sorted(
                (pos, next(iter(bases)))
                for pos, bases in calls.items()
                if len(bases) == 1
            )
        )
        if pattern:
            counts[chrom][pattern] += 1
    return counts


class PhaseForest:
    """Union with a parity bit from each node to its parent."""

    def __init__(self, positions):
        self.parent = {p: p for p in positions}
        self.bit = dict.fromkeys(positions, 0)
        self.size = dict.fromkeys(positions, 1)

    def find(self, p):
        if self.parent[p] != p:
            root, up = self.find(self.parent[p])
            self.bit[p] ^= up
            self.parent[p] = root
        return self.parent[p], self.bit[p]

    def join(self, a, b, parity):
        ra, pa = self.find(a)
        rb, pb = self.find(b)
        if ra == rb:
            return "consistent_cycle" if pa ^ pb == parity else "conflicting_cycle"
        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.bit[rb] = pa ^ pb ^ parity
        self.size[ra] += self.size[rb]
        return "joined"


def phase(sites, patterns, min_support=2, min_margin=1):
    rows = []
    edges = []
    for chrom, variants in sites.items():
        votes = defaultdict(lambda: [0, 0])
        coverage = Counter()
        for pattern, count in patterns.get(chrom, {}).items():
            for pos, allele in pattern:
                coverage[pos] += count
            for (a, x), (b, y) in zip(pattern, pattern[1:]):
                votes[a, b][x ^ y] += count
        forest = PhaseForest(variants)
        ordered = sorted(
            votes.items(), key=lambda v: (-abs(v[1][0] - v[1][1]), -max(v[1]), v[0])
        )
        for (a, b), (same, opposite) in ordered:
            if same == opposite:
                state = "tie"
            elif max(same, opposite) < min_support or abs(same - opposite) < min_margin:
                state = "low_support"
            else:
                state = forest.join(a, b, int(opposite > same))
            edges.append((chrom, a, b, same, opposite, state))
        groups = defaultdict(list)
        for pos in sorted(variants):
            root, parity = forest.find(pos)
            groups[root].append((pos, parity))
        for block in sorted(groups.values(), key=lambda z: z[0][0]):
            first, anchor = block[0]
            linked = len(block) > 1
            for pos, b in block:
                ref, alt = variants[pos]
                gt = b ^ anchor
                rows.append(
                    (
                        chrom,
                        pos,
                        ref,
                        alt,
                        first if linked else ".",
                        gt if linked else ".",
                        1 - gt if linked else ".",
                        coverage[pos],
                    )
                )
    ranks = {chrom: i for i, chrom in enumerate(sites)}
    return sorted(rows, key=lambda r: (ranks[r[0]], r[1])), edges


@contextmanager
def atomic_stream(path):
    path = Path(path)
    with NamedTemporaryFile("w", dir=path.parent, encoding="utf-8", delete=False) as f:
        tmp = f.name
        try:
            yield f
            f.flush()
            os.fsync(f.fileno())
            f.close()
            os.replace(tmp, path)
        finally:
            if os.path.exists(tmp):
                os.unlink(tmp)


def table(path, header, rows):
    with atomic_stream(path) as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(map(str, row)) + "\n")


def legacy_outputs(directory, sites, patterns, rows):
    grouped = defaultdict(list)
    for r in rows:
        if r[4] != ".":
            grouped[r[0], r[4]].append(r)
    with atomic_stream(directory / "Test_first_merging_outputs.txt") as out:
        for chrom, group in patterns.items():
            for pattern, n in sorted(group.items()):
                labels = [str(p) + sites[chrom][p][b] for p, b in pattern]
                print(f"[{chrom};{labels};{n}]", file=out)
    with atomic_stream(directory / "test_Second_merging_HASSAN.txt") as out:
        idx = 1
        for (_, start), block in grouped.items():
            print(
                f"BLOCK: offset: {idx} len: {len(block)} phased: {len(block)} SPAN: {block[-1][1]-start}",
                file=out,
            )
            for chrom, pos, ref, alt, _, h1, h2, _ in block:
                print(idx, h1, h2, chrom, pos, ref, alt, "0/1", sep="\t", file=out)
                idx += 1
            print("********", file=out)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", default="napus.vcf")
    p.add_argument(
        "--bam", default="napus_sorted.bam", help="BAM or SAM; an index is optional"
    )
    p.add_argument("--sample")
    p.add_argument("--out-dir", type=Path, default=Path("."))
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--min-baseq", type=int, default=20)
    p.add_argument(
        "--min-support",
        type=int,
        default=2,
        help="Minimum number of fragments supporting the winning parity",
    )
    p.add_argument("--min-margin", type=int, default=1)
    p.add_argument(
        "--legacy",
        action="store_true",
        help="Also write the original two report filenames",
    )
    import sys

    args = list(sys.argv[1:] if argv is None else argv)
    a = p.parse_args(args)
    if min(a.min_mapq, a.min_baseq) < 0 or min(a.min_support, a.min_margin) < 1:
        p.error(
            "Quality thresholds must be nonnegative and evidence thresholds positive"
        )
    try:
        sample, sites = load_sites(a.vcf, a.sample)
        patterns = fragment_patterns(a.bam, sites, sample, a.min_mapq, a.min_baseq)
        rows, edges = phase(sites, patterns, a.min_support, a.min_margin)
        a.out_dir.mkdir(parents=True, exist_ok=True)
        names = ["phased.tsv", "links.tsv"]
        if a.legacy or not args:
            names += [
                "Test_first_merging_outputs.txt",
                "test_Second_merging_HASSAN.txt",
            ]
        if {Path(a.bam).resolve(), Path(a.vcf).resolve()} & {
            (a.out_dir / n).resolve() for n in names
        }:
            raise ValueError("An output would replace an input")
        table(
            a.out_dir / "phased.tsv",
            ("contig", "pos", "ref", "alt", "block_start", "hap1", "hap2", "fragments"),
            rows,
        )
        table(
            a.out_dir / "links.tsv",
            ("contig", "left", "right", "same", "opposite", "decision"),
            edges,
        )
        if a.legacy or not args:
            legacy_outputs(a.out_dir, sites, patterns, rows)
    except (OSError, ValueError) as error:
        p.exit(2, f"{p.prog}: {error}\n")
    print(
        f"{sample}: {len(rows)} heterozygous SNVs; {sum(r[4]!='.' for r in rows)} in linked blocks"
    )


if __name__ == "__main__":
    main()
