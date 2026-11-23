#!/usr/bin/env python3
"""Pure-stdlib smoke test for sample-level read aggregation (no real data)."""
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "aggregate_sample_level_validation.py"


def write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def read(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))


class AggregationTest(unittest.TestCase):
    def test_one_strict_read_per_run_combines_to_two(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            result=root/"results"
            manifest=root/"manifest.tsv"
            runs=[("RUN_A","ONT"),("RUN_B","PACBIO"),("RUN_C","ONT")]
            write(manifest,[dict(sample="S1",species="Species_one",run=run,
                                 platform=platform,reference="/ref.fa",
                                 besthit="/candidates.tsv") for run,platform in runs])
            for run,_ in runs[:2]:
                base=result/"S1"/run
                nuc=base/"chrm_candidate_validation"
                chrm=base/"longread_chrM_homology"
                loc=dict(locus_id="L1", nuclear_chrom="chr1",
                         nuclear_start0="100",nuclear_end0="1200",
                         chrM_hit="chrM",chrM_start="500",chrM_end="1200",
                         group="BLAST_QC_PASS",n_strict_sequence_support="1")
                write(nuc/"chrM_candidate_locus_validation.tsv",[loc])
                write(nuc/"chrM_candidate_read_evidence.tsv",
                      [dict(locus_id="L1",read="reused_read_name",mapq="60",
                            locus_identity_to_nuclear_ref="0.98",
                            locus_fraction="0.97",left_aligned="100",
                            right_aligned="100",max_internal_deletion="0",
                            strict_sequence_support="1")])
                write(chrm/"per_locus_chrM_homology.tsv",
                      [dict(locus_id="L1",chrM_status="NOT_SUPPORTED")])
                write(chrm/"per_read_chrM_homology.tsv",
                      [dict(locus_id="L1",read="reused_read_name",
                            expected_chrM_homology="1")])
            out=root/"summary"
            cmd=[sys.executable,str(SCRIPT),"--manifest",str(manifest),
                 "--results-root",str(result),"--outdir",str(out)]
            subprocess.run(cmd,check=True,capture_output=True,text=True)
            loc=read(out/"sample_locus_validation.tsv")
            self.assertEqual(len(loc),1)
            self.assertEqual(loc[0]["n_nuclear_strict_reads"],"2")
            self.assertEqual(loc[0]["n_expected_chrM_reads"],"2")
            self.assertEqual(loc[0]["final_status"],"VALIDATED")
            groups=read(out/"sample_platform_locus_validation.tsv")
            self.assertEqual(len(groups),2)
            self.assertTrue(all(g["final_status"]=="NOT_VALIDATED" for g in groups))
            sample=read(out/"sample_summary.tsv")[0]
            self.assertEqual(sample["completeness"],"PARTIAL")
            self.assertEqual(sample["n_completed_runs"],"2")
            self.assertEqual(len(read(out/"run_status.tsv")),3)
            # Fail-closed when a final publication requires all runs complete.
            p=subprocess.run(cmd+["--require-complete"],capture_output=True,text=True)
            self.assertNotEqual(p.returncode,0)
            self.assertIn("pending",p.stderr.lower()+p.stdout.lower())


if __name__=="__main__":
    unittest.main()
