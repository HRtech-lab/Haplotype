"""Make a small artificial diploid marker dataset for the phasing example."""

from pathlib import Path
import pysam


def make_fixture(directory):
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    vcf = root / "markers.vcf"
    bam = root / "reads.bam"
    variants = [(5, "A", "G"), (10, "C", "T"), (17, "G", "A"), (27, "C", "T")]
    with vcf.open("w") as f:
        f.write(
            '##fileformat=VCFv4.2\n##contig=<ID=markers,length=48>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        )
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tplant\n")
        for pos, ref, alt in variants:
            f.write(f"markers\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t1/0\n")
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": "markers", "LN": 48}],
        "RG": [{"ID": "set1", "SM": "plant"}],
    }
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for i, alleles in enumerate([(0, 1, 0), (1, 0, 1)] * 3):
            read = pysam.AlignedSegment(out.header)
            read.query_name = f"read{i}"
            read.reference_id = 0
            read.reference_start = 0
            read.flag = 0
            read.mapping_quality = 55
            read.cigarstring = "2S7M1I5M1D7M"
            bases = list("C" * 22)
            for q, r in read.get_aligned_pairs(matches_only=True):
                for (pos, ref, alt), index in zip(variants[:3], alleles):
                    if r + 1 == pos:
                        bases[q] = (ref, alt)[index]
            read.query_sequence = "".join(bases)
            read.query_qualities = [35] * len(bases)
            read.set_tag("RG", "set1")
            out.write(read)
    return vcf, bam
