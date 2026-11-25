"""Regression checks for aligned coordinates, overlapping mates and evidence ties."""

import tempfile, unittest
from pathlib import Path
from collections import Counter
import pysam
from demo_data import make_fixture
from haplotype import load_sites, fragment_patterns, phase, main


class PhaseTests(unittest.TestCase):
    def test_indel_softclip_final_read_and_unlinked_site(self):
        with tempfile.TemporaryDirectory() as d:
            v, b = make_fixture(d)
            sample, sites = load_sites(v)
            patterns = fragment_patterns(b, sites, sample)
            self.assertEqual(sum(patterns["markers"].values()), 6)
            rows, edges = phase(sites, patterns)
            self.assertEqual(
                [(r[1], r[5], r[6]) for r in rows],
                [(5, 0, 1), (10, 1, 0), (17, 0, 1), (27, ".", ".")],
            )
            main(["--bam", str(b), "--vcf", str(v), "--out-dir", d, "--legacy"])
            original = (Path(d) / "phased.tsv").read_bytes()
            main(["--bam", str(b), "--vcf", str(v), "--out-dir", d])
            self.assertEqual(original, (Path(d) / "phased.tsv").read_bytes())

    def test_pairs_overlap_conflict_and_sample_filter(self):
        with tempfile.TemporaryDirectory() as d:
            v, _ = make_fixture(d)
            sample, sites = load_sites(v)
            bam = Path(d) / "mates.bam"
            h = {
                "SQ": [{"SN": "markers", "LN": 48}],
                "RG": [{"ID": "a", "SM": "plant"}, {"ID": "b", "SM": "other"}],
            }
            with pysam.AlignmentFile(str(bam), "wb", header=h) as out:
                for flag, base, group, name in [
                    (99, "A", "a", "pair"),
                    (147, "G", "a", "pair"),
                    (0, "G", "b", "other"),
                    (1024, "G", "a", "duplicate"),
                ]:
                    r = pysam.AlignedSegment(out.header)
                    r.query_name = name
                    r.flag = flag
                    r.reference_id = 0
                    r.reference_start = 0
                    r.next_reference_id = 0
                    r.next_reference_start = 0
                    r.cigarstring = "12M"
                    r.mapping_quality = 50
                    r.query_sequence = "CCCC" + base + "CCCC" + "T" + "CC"
                    r.query_qualities = [35] * 12
                    r.set_tag("RG", group)
                    out.write(r)
            result = fragment_patterns(bam, sites, sample)
            self.assertEqual(result["markers"], {((10, 1),): 1})

    def test_cycle_tie_and_empty(self):
        sites = {"ctg": {i: ("A", "G") for i in [1, 2, 3, 4]}}
        pairs = Counter(
            {
                ((1, 0), (2, 0)): 6,
                ((2, 0), (3, 1)): 5,
                ((1, 0), (3, 0)): 3,
                ((3, 0), (4, 0)): 2,
                ((3, 0), (4, 1)): 2,
            }
        )
        rows, links = phase(sites, {"ctg": pairs})
        self.assertEqual([r[5] for r in rows], [0, 0, 1, "."])
        self.assertEqual({e[-1] for e in links}, {"joined", "conflicting_cycle", "tie"})
        self.assertEqual(phase({}, {}), ([], []))


if __name__ == "__main__":
    unittest.main()
