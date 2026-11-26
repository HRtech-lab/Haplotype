# Haplotype

A **haplotype** is the combination of alleles carried together on one chromosome copy. A diploid genotype records the two alleles at a site. Phasing determines how alleles at different sites are arranged across the two copies. Two sites can both have genotype `0/1`, yet the underlying pair of haplotypes may be `00 / 11` or `01 / 10`. Reads spanning the sites can distinguish these possibilities. [NHGRI glossary](https://www.genome.gov/genetics-glossary/haplotype)

The program phases diploid, biallelic heterozygous SNVs from a single sample using existing short-read alignments. It produces evidence counts and separate phase blocks. The worked example has known answers and includes an insertion, a deletion, soft clipping and an unlinked fourth marker. [Notebook and figures](Worked%20Examples/01%20Read-backed%20phasing/README.md)

![Observed alleles and coverage](Worked%20Examples/01%20Read-backed%20phasing/figure.svg)

## Try the included example

From the repository root, run the worked example on its generated fixture. Its output shows three markers in one phase block and a fourth marker whose phase is unresolved.

```sh
python -m pip install -r requirements.txt
python "Worked Examples/01 Read-backed phasing/read_backed_phasing.py"
```

## Run on your own alignments

```sh
python haplotype.py --vcf calls.vcf.gz --bam reads.bam --sample plant --out-dir results
python -m unittest -v
```

BAM or SAM and VCF or BGZF-compressed VCF are accepted, with matching contig names. Input indexes are optional. The sample is inferred from a single-sample VCF. Where BAM read groups have sample names, only groups assigned to the selected sample contribute. For an alignment file without sample tags, provide reads from that sample alone.

## Following the evidence

1. The VCF provides PASS or unfiltered, diploid `0/1`, `1/0` or equivalent phased calls at single-base substitutions. Input phase is ignored. VCF positions are one-based. Aligned query/reference pairs from pysam map observations through the full CIGAR string.
2. A fragment contributes once at each marker. Primary unpaired reads and proper, same-contig pairs can contribute. Duplicate, failed-QC, secondary and supplementary alignments are skipped. Overlapping mates that disagree leave that site unresolved for that fragment.
3. Consecutive observed markers vote for equal or opposite allele indices. Links are considered by vote margin, then winning support and position. A parity forest joins consistent blocks and records contradictory cycles. Ties and weak links remain explicit in the link report.
4. The lowest coordinate in each connected block anchors haplotype 1 to reference. A block has no parental orientation without independent evidence. Unlinked sites remain in the output with missing phase.

The default mapping and base quality thresholds are both 20. The default link requires two fragments supporting its winning orientation and a margin of at least one vote. `--help` lists the controls. The forest is a transparent heuristic. Its edge counts do not constitute calibrated confidence probabilities. The program retains informative fragments in memory and is intended for small datasets or preselected regions. Diploid phase calls do not establish the phase of a polyploid genome. [Pysam aligned pairs](https://pysam.readthedocs.io/en/latest/api.html) · [SAM/VCF specifications](https://samtools.github.io/hts-specs/) · [WhatsHap methods and usage](https://whatshap.readthedocs.io/en/latest/guide.html)

## Reports

| File | Meaning |
|---|---|
| `phased.tsv` | One row per retained variant. Allele indices 0 and 1 refer to REF and ALT. A dot marks unresolved phase. |
| `links.tsv` | Same/opposite fragment counts and the disposition of each link. |

With no command-line arguments, the script reads the earlier `napus.vcf` and `napus_sorted.bam` filenames and also emits both legacy-named reports. Re-running replaces each report atomically. The modern tables preserve explicit reference and alternate columns. The historical report style is available with `--legacy`; its genotype field is normalised to `0/1`.

## Sources and reuse

The example marker data are generated in `demo_data.py`. New code and teaching materials in this revision use the [MIT licence](LICENSE). The definitions, file formats and phasing literature are linked above.
