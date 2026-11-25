#!/usr/bin/env python3
"""Offline tests for cohort manager metadata, status, and retry classification."""
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts import cohort_manager as m


class CohortManagerTests(unittest.TestCase):
    def row(self, run="SRR1", preset="map-ont"):
        return dict(sample="SRS1",species="Species_one",run=run,
                    platform="ONT",preset=preset,fastq_ftp="ftp.sra.ebi.ac.uk/test.fastq.gz",
                    reference="/ref.fa",besthit="/besthit.tsv")

    def test_manifest_deduplicates_runs_and_rejects_unsafe_presets(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"manifest.tsv"
            import csv
            def save(rows):
                with f.open("w",newline="") as out:
                    w=csv.DictWriter(out,fieldnames=m.FIELDS,delimiter="\t")
                    w.writeheader()
                    w.writerows(rows)
            save([self.row(),self.row()])
            self.assertEqual(len(m.load_manifest([f])),1)
            save([self.row(preset="REVIEW")])
            with self.assertRaises(ValueError):
                m.load_manifest([f])

    def test_slurm_array_job_active_and_failed(self):
        completed = lambda rc=0, out="", err="": subprocess.CompletedProcess(
            args=["squeue"], returncode=rc, stdout=out, stderr=err)
        with patch.object(m.subprocess, "run", return_value=completed(out="RUNNING\n")) as mock:
            self.assertEqual(m.query_slurm("123_1"),("active","RUNNING",None))
            mock.assert_called_once()
        with patch.object(m.subprocess, "run", return_value=completed()), \
             patch.object(m, "cmd", return_value="123|FAILED|1:0\n123_1|FAILED|75:0\n"):
            self.assertEqual(m.query_slurm("123_1"),("terminal","FAILED",75))
        # Bouchet: finished array task is removed from squeue with nonzero
        # exit code and "Invalid job id specified". Must fall back to sacct.
        with patch.object(m.subprocess, "run", return_value=completed(
                1, err="slurm_load_jobs error: Invalid job id specified")), \
             patch.object(m, "cmd", return_value="123|COMPLETED|0:0\n123_1|FAILED|75:0\n") as acct:
            self.assertEqual(m.query_slurm("123_1"),("terminal","FAILED",75))
            self.assertEqual(acct.call_args.args[0][4],"123")
        with patch.object(m.subprocess, "run", return_value=completed(
                1, err="slurm_load_jobs error: Invalid job id specified")), \
             patch.object(m, "cmd", return_value=""):
            self.assertEqual(m.query_slurm("123_1"),("unknown","NO_ACCOUNTING",None))
        with patch.object(m.subprocess, "run", return_value=completed(
                1, err="slurm_load_jobs error: Unauthorized")):
            with self.assertRaises(RuntimeError):
                m.query_slurm("123_1")

    def test_parent_only_sacct_record(self):
        """Bouchet may omit 123_1 while reporting parent 123 FAILED 75."""
        import subprocess
        queue=subprocess.CompletedProcess(args=["squeue"],returncode=1,
            stdout="",stderr="slurm_load_jobs error: Invalid job id specified")
        with patch.object(m.subprocess, "run", return_value=queue), \
             patch.object(m, "cmd", return_value=(
                 "123|FAILED|75:0\n"
                 "123.batch|FAILED|75:0\n"
                 "123.extern|COMPLETED|0:0\n")):
            self.assertEqual(m.query_slurm("123_1"),("terminal","FAILED",75))

    def test_http_403_detection_is_job_specific(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            (root/"logs").mkdir()
            (root/"logs"/"lr_stream_123_1.log").write_text(
                "Download aborted; status=403; partial=0 bytes\n")
            (root/"logs"/"lr_stream_124_1.log").write_text(
                "Download aborted; protocol error\n")
            self.assertTrue(m.is_remote_403_failure("123_1",root))
            self.assertFalse(m.is_remote_403_failure("124_1",root))
            self.assertFalse(m.is_remote_403_failure("125_1",root))

    def test_completed_outputs_are_detected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            row=self.row()
            base=root/"SRS1"/"SRR1"
            (base/"chrm_candidate_validation").mkdir(parents=True)
            (base/"longread_chrM_homology").mkdir(parents=True)
            for rel in ("chrm_candidate_validation/chrM_candidate_locus_validation.tsv",
                        "chrm_candidate_validation/chrM_candidate_read_evidence.tsv",
                        "longread_chrM_homology/per_locus_chrM_homology.tsv",
                        "longread_chrM_homology/per_read_chrM_homology.tsv"):
                target = base / rel
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text("header\n")
            self.assertTrue(m.output_complete(row,"validate",root))
            self.assertFalse(m.output_complete(row,"align",root))

    def test_merge_complete_requires_all_runs_and_verified_counts(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);target=root/"SRS1"/"ONT"/"numt_target_merged.bam"
            target.parent.mkdir(parents=True)
            target.write_bytes(b"test")
            Path(str(target)+".bai").write_bytes(b"test")
            audit=target.parent/"numt_target_merged.audit.json"
            payload=dict(sample="SRS1",platform="ONT",partial=False,
                         included_runs=["SRR1"],n_included_runs=1,
                         target_record_counts_verified=True,
                         per_run_target_records={"SRR1":10})
            audit.write_text(json.dumps(payload))
            self.assertTrue(m.merge_complete(("SRS1","ONT"),[self.row()],root))
            self.assertFalse(m.merge_complete(("SRS1","ONT"),
                                              [self.row(),self.row("SRR2")],root))
            payload["partial"]=True
            audit.write_text(json.dumps(payload))
            self.assertFalse(m.merge_complete(("SRS1","ONT"),[self.row()],root))


if __name__=="__main__":
    unittest.main()
