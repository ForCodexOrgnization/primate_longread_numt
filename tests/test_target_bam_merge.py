#!/usr/bin/env python3
"""Unit checks for safe sample-platform targeted BAM planning (no BAM needed)."""
import csv
import tempfile
import unittest
from pathlib import Path

from scripts.merge_sample_platform_numt_bam import merge_intervals, gather


def table(path, rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader()
        w.writerows(rows)


class TargetMergeTests(unittest.TestCase):
    def test_intervals_union_only_within_contig(self):
        loci=[("chr2",10,20),("chr1",200,230),("chr1",100,150),
              ("chr1",145,180),("chr2",19,23)]
        self.assertEqual(merge_intervals(loci),
                         [("chr1",100,180),("chr1",200,230),("chr2",10,23)])

    def test_nonready_runs_fail_closed(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            r=dict(sample="S1",species="Species_one",run="RUN_A",
                   platform="ONT",reference="/ref.fa",besthit="/besthit.tsv")
            with self.assertRaises(RuntimeError):
                gather([r],root,1000,"chrM",False)
            with self.assertRaises(RuntimeError):
                gather([r],root,1000,"chrM",True)

    def test_input_candidate_intervals_and_partial(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            def row(name):
                return dict(sample="S1",species="Species_one",run=name,
                            platform="ONT",reference="/ref.fa",besthit="/besthit.tsv")
            base=root/"S1"/"RUN_A"
            base.mkdir(parents=True)
            (base/"complete.flag").touch()
            for suffix in ("aligned.bam","aligned.bam.bai"):
                (base/suffix).write_bytes(b"not_a_real_bam")
            table(base/"chrm_candidate_validation"/"chrM_candidate_locus_validation.tsv",
                  [dict(locus_id="L1",nuclear_chrom="chr1",
                        nuclear_start0="100",nuclear_end0="300")])
            with self.assertRaises(RuntimeError):
                gather([row("RUN_A"),row("RUN_B")],root,1000,"chrM",False)
            available,missing,intervals=gather(
                [row("RUN_A"),row("RUN_B")],root,1000,"chrM",True)
            self.assertEqual(len(available),1)
            self.assertEqual(len(missing),1)
            self.assertEqual(intervals,[("chr1",100,300)])


if __name__=="__main__":
    unittest.main()
