from pathlib import Path
import tempfile
import shutil
from unittest import TestCase
from types import SimpleNamespace
from urllib.error import HTTPError

import pandas as pd

from pyhiv.loading.reference_update import (
    CrfReference,
    GenBankRecordPartial,
    GenBankRecordUnavailable,
    LanlAlignment,
    LanlReferenceRecord,
    ReconciliationError,
    accessions_from_text,
    build_desired_reference_map,
    extract_lanl_alignment_content,
    load_lanl_crf_references,
    parse_reference_name,
    rebuild_sequences_with_locations,
    update_lanl_reference_alignment,
    update_reference_dataset,
    validate_hiv1_complete_or_near_complete,
)


class TestReferenceUpdate(TestCase):

    def setUp(self):
        self.tmp_dir = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: shutil.rmtree(self.tmp_dir, ignore_errors=True))
        self.references_dir = self.tmp_dir / "reference_fastas"
        self.references_dir.mkdir()
        (self.tmp_dir / "HXB2_fasta").mkdir()
        self.sequences_path = self.tmp_dir / "sequences_with_locations.tsv"

    def write_fasta(self, name, sequence="ACGTACGT"):
        path = self.references_dir / name
        path.write_text(f">{Path(name).stem}\n{sequence}\n", encoding="utf-8")
        return path

    def write_table(self, rows):
        pd.DataFrame(
            rows,
            columns=[
                "subtype",
                "country",
                "year",
                "isolate",
                "accession",
                "align sequence",
                "sequence",
                "features",
                "group",
            ],
        ).to_csv(self.sequences_path, sep="\t", index=False)

    def load_crf_html(self, html):
        class FakeResponse:
            def __enter__(self):
                return self

            def __exit__(self, *args):
                return None

            def read(self):
                return html.encode("utf-8")

        from unittest.mock import patch
        with patch("pyhiv.loading.reference_update.urlopen", return_value=FakeResponse()):
            return load_lanl_crf_references()

    def test_parse_reference_name_keeps_accession_and_subtype(self):
        parsed = parse_reference_name("PP056116-159_01103")

        self.assertEqual(parsed["accession"], "PP056116")
        self.assertEqual(parsed["subtype"], "159_01103")

    def test_rebuild_adds_reference_folder_entries_with_none_features(self):
        self.write_fasta("AB123456-A1.fasta")
        self.write_table([])

        result = rebuild_sequences_with_locations(
            self.references_dir,
            self.sequences_path,
        )

        table = pd.read_csv(self.sequences_path, sep="\t", dtype=str, keep_default_na=False)
        self.assertEqual(result["added_rows"], 1)
        self.assertEqual(table.loc[0, "accession"], "AB123456")
        self.assertEqual(table.loc[0, "subtype"], "A1")
        self.assertEqual(table.loc[0, "features"], "None")

    def test_update_reports_up_to_date(self):
        self.write_fasta("AB123456-A1.fasta")
        self.write_table([
            {
                "subtype": "A1",
                "country": "PT",
                "year": "2024",
                "isolate": "sample",
                "accession": "AB123456",
                "align sequence": "ACGTACGT",
                "sequence": "ACGTACGT",
                "features": "{'gag': (1, 4)}",
                "group": "M",
            }
        ])

        result = update_reference_dataset(
            base_dir=self.tmp_dir,
            update_lanl_alignment=False,
            crf_loader=lambda: {
                "A1": CrfReference("A1", "CRFA1", "sample", "AB123456")
            },
        )

        self.assertTrue(result.up_to_date)
        self.assertEqual(result.total_references, 1)

    def test_update_adds_missing_crf_sequence_and_tsv_row(self):
        self.write_table([])

        result = update_reference_dataset(
            base_dir=self.tmp_dir,
            email="test@example.com",
            update_lanl_alignment=False,
            crf_loader=lambda: {
                "205_0107": CrfReference("205_0107", "CRF205_0107", "20HZ2132", "PX262425")
            },
            genbank_fetcher=lambda accession: ("ACGTACGT", None),
            sleep_seconds=0,
        )

        table = pd.read_csv(self.sequences_path, sep="\t", dtype=str, keep_default_na=False)
        self.assertFalse(result.up_to_date)
        self.assertEqual(result.added_sequences, 1)
        self.assertEqual(result.added_rows, 1)
        self.assertTrue((self.references_dir / "PX262425-205_0107.fasta").exists())
        self.assertEqual(table.loc[0, "accession"], "PX262425")
        self.assertEqual(table.loc[0, "features"], "None")

    def test_update_skips_unavailable_genbank_accession(self):
        self.write_table([])

        def failing_fetcher(accession):
            raise HTTPError(
                url=f"https://example.test/{accession}",
                code=400,
                msg="Bad Request",
                hdrs=None,
                fp=None,
            )

        result = update_reference_dataset(
            base_dir=self.tmp_dir,
            email="test@example.com",
            update_lanl_alignment=False,
            crf_loader=lambda: {
                "151_0107": CrfReference("151_0107", "CRF151_0107", "220565", "PP982937")
            },
            genbank_fetcher=failing_fetcher,
            sleep_seconds=0,
        )

        self.assertEqual(result.added_sequences, 0)
        self.assertEqual(result.failed_sequences, 0)
        self.assertEqual(result.failed_accessions, ())
        self.assertEqual(result.pending_ncbi_references, 1)
        self.assertEqual(
            result.pending_ncbi_accessions,
            ("151_0107:PP982937:authoritative_reference_pending_ncbi",),
        )
        self.assertFalse((self.references_dir / "PP982937-151_0107.fasta").exists())

    def test_update_skips_empty_genbank_record(self):
        self.write_table([])

        result = update_reference_dataset(
            base_dir=self.tmp_dir,
            email="test@example.com",
            update_lanl_alignment=False,
            crf_loader=lambda: {
                "151_0107": CrfReference("151_0107", "CRF151_0107", "220565", "PP982937")
            },
            genbank_fetcher=lambda accession: (_ for _ in ()).throw(
                GenBankRecordUnavailable(accession)
            ),
            sleep_seconds=0,
        )

        self.assertEqual(result.added_sequences, 0)
        self.assertEqual(result.failed_sequences, 0)
        self.assertEqual(result.failed_accessions, ())
        self.assertEqual(result.pending_ncbi_references, 1)
        self.assertEqual(
            result.pending_ncbi_accessions,
            ("151_0107:PP982937:authoritative_reference_pending_ncbi",),
        )
        self.assertFalse((self.references_dir / "PP982937-151_0107.fasta").exists())

    def test_pending_ncbi_is_idempotent_and_not_a_failed_sequence(self):
        self.write_table([])
        crfs = {
            "151_0107": (
                CrfReference("151_0107", "CRF151_0107", "220564", "PP982936"),
                CrfReference("151_0107", "CRF151_0107", "220565", "PP982937"),
                CrfReference("151_0107", "CRF151_0107", "ZLQ04681", "PP982938"),
            )
        }

        def pending_fetcher(accession):
            raise GenBankRecordUnavailable(accession)

        first = update_reference_dataset(
            base_dir=self.tmp_dir,
            email="test@example.com",
            update_lanl_alignment=False,
            crf_loader=lambda: crfs,
            genbank_fetcher=pending_fetcher,
            sleep_seconds=0,
        )
        second = update_reference_dataset(
            base_dir=self.tmp_dir,
            email="test@example.com",
            update_lanl_alignment=False,
            crf_loader=lambda: crfs,
            genbank_fetcher=pending_fetcher,
            sleep_seconds=0,
        )

        self.assertEqual(first.added_sequences, 0)
        self.assertEqual(second.added_sequences, 0)
        self.assertEqual(second.added_rows, 0)
        self.assertEqual(second.updated_rows, 0)
        self.assertEqual(second.reclassified_references, 0)
        self.assertEqual(second.lanl_fasta_files_removed, 0)
        self.assertEqual(second.failed_sequences, 0)
        self.assertEqual(second.pending_ncbi_references, 3)
        self.assertFalse(second.up_to_date)
        self.assertEqual(
            second.pending_ncbi_accessions,
            (
                "151_0107:PP982936:authoritative_reference_pending_ncbi",
                "151_0107:PP982937:authoritative_reference_pending_ncbi",
                "151_0107:PP982938:authoritative_reference_pending_ncbi",
            ),
        )

    def test_update_skips_partial_genbank_record(self):
        self.write_table([])

        result = update_reference_dataset(
            base_dir=self.tmp_dir,
            email="test@example.com",
            update_lanl_alignment=False,
            crf_loader=lambda: {
                "151_0107": CrfReference("151_0107", "CRF151_0107", "220565", "PP982937")
            },
            genbank_fetcher=lambda accession: (_ for _ in ()).throw(
                GenBankRecordPartial(accession)
            ),
            sleep_seconds=0,
        )

        self.assertEqual(result.added_sequences, 0)
        self.assertEqual(result.failed_sequences, 1)
        self.assertEqual(
            result.failed_accessions,
            ("151_0107:PP982937:not_complete_or_near_complete",),
        )
        self.assertFalse((self.references_dir / "PP982937-151_0107.fasta").exists())

    def test_genbank_validation_rejects_short_partial_sequence(self):
        record = SimpleNamespace(
            id="AB123456",
            name="AB123456",
            description="HIV-1 isolate sample partial genome",
        )

        with self.assertRaises(GenBankRecordPartial):
            validate_hiv1_complete_or_near_complete(record, "A" * 1200)

    def test_genbank_validation_accepts_partial_genome_when_near_full_length(self):
        record = SimpleNamespace(
            id="OP921950",
            name="OP921950",
            description="HIV-1 isolate sample partial genome",
        )

        validate_hiv1_complete_or_near_complete(record, "A" * 8974)

    def test_genbank_validation_accepts_long_sequence_with_gene_annotation(self):
        record = SimpleNamespace(
            id="OR260538",
            name="OR260538",
            description="HIV-1 isolate sample pol gene partial CDS",
        )

        validate_hiv1_complete_or_near_complete(record, "A" * 9000)

    def test_genbank_validation_accepts_nearly_complete_sequence(self):
        record = SimpleNamespace(
            id="AB123456",
            name="AB123456",
            description="HIV-1 isolate sample complete genome",
        )

        validate_hiv1_complete_or_near_complete(record, "A" * 9000)

    def test_refresh_features_keeps_none_when_genbank_record_is_empty(self):
        self.write_fasta("AB123456-A1.fasta")
        self.write_table([
            {
                "subtype": "A1",
                "country": "PT",
                "year": "2024",
                "isolate": "sample",
                "accession": "AB123456",
                "align sequence": "ACGTACGT",
                "sequence": "ACGTACGT",
                "features": "None",
                "group": "M",
            }
        ])

        result = rebuild_sequences_with_locations(
            self.references_dir,
            self.sequences_path,
            refresh_features=True,
            fetcher=lambda accession: (_ for _ in ()).throw(
                GenBankRecordUnavailable(accession)
            ),
            sleep_seconds=0,
        )

        table = pd.read_csv(self.sequences_path, sep="\t", dtype=str, keep_default_na=False)
        self.assertEqual(result["updated_feature_rows"], 1)
        self.assertEqual(table.loc[0, "features"], "None")

    def test_update_syncs_lanl_alignment_before_crf_check(self):
        self.write_table([])
        alignment = LanlAlignment(
            year=2025,
            content=(
                "Ref.A1.PT.2024.sample.AB123456,AC-GT\n"
                "Ref.205_0107.CN.2024.20HZ2132.PX262425,ACGT--\n"
            ),
        )

        result = update_reference_dataset(
            base_dir=self.tmp_dir,
            lanl_alignment_loader=lambda: alignment,
            crf_loader=lambda: {
                "205_0107": CrfReference("205_0107", "CRF205_0107", "20HZ2132", "PX262425")
            },
        )

        table = pd.read_csv(self.sequences_path, sep="\t", dtype=str, keep_default_na=False)
        self.assertEqual(result.lanl_alignment_year, 2025)
        self.assertTrue(result.lanl_alignment_updated)
        self.assertEqual(result.lanl_fasta_files_added, 2)
        self.assertEqual(result.added_sequences, 0)
        self.assertTrue((self.tmp_dir / "HIV1_REF_2025_genome_DNA.csv").exists())
        self.assertTrue((self.references_dir / "AB123456-A1.fasta").exists())
        self.assertTrue((self.references_dir / "PX262425-205_0107.fasta").exists())
        self.assertEqual(set(table["accession"]), {"AB123456", "PX262425"})

    def test_lanl_alignment_sync_includes_multiple_valid_crf_references(self):
        self.write_table([])
        alignment = LanlAlignment(
            year=2026,
            content=(
                "Ref.120_0107.CN.2015.CN_SZ_MSM_MSM_LS15083.ON351495,AAAA\n"
                "Ref.120_0107.CN.2015.HES_LS16846.OK662596,CCCC\n"
                "Ref.120_0107.CN.2015.HES_LS17858.OK662597,GGGG\n"
            ),
        )

        result = update_reference_dataset(
            base_dir=self.tmp_dir,
            lanl_alignment_loader=lambda: alignment,
            crf_loader=lambda: {},
        )

        table = pd.read_csv(self.sequences_path, sep="\t", dtype=str, keep_default_na=False)
        self.assertEqual(result.lanl_fasta_files_added, 3)
        self.assertEqual(result.total_references, 3)
        self.assertEqual(
            set(table["accession"]),
            {"ON351495", "OK662596", "OK662597"},
        )
        self.assertTrue((self.references_dir / "ON351495-120_0107.fasta").exists())
        self.assertTrue((self.references_dir / "OK662596-120_0107.fasta").exists())
        self.assertTrue((self.references_dir / "OK662597-120_0107.fasta").exists())

    def test_crf_database_three_near_complete_genomes_adds_three_references(self):
        self.write_table([])
        crfs = {
            "120_0107": (
                CrfReference("120_0107", "CRF120_0107", "a", "OK662596"),
                CrfReference("120_0107", "CRF120_0107", "b", "OK662597"),
                CrfReference("120_0107", "CRF120_0107", "c", "ON351495"),
            )
        }

        result = update_reference_dataset(
            base_dir=self.tmp_dir,
            email="test@example.com",
            update_lanl_alignment=False,
            crf_loader=lambda: crfs,
            genbank_fetcher=lambda accession: ("A" * 9000, None),
            sleep_seconds=0,
        )

        table = pd.read_csv(self.sequences_path, sep="\t", dtype=str, keep_default_na=False)
        self.assertEqual(result.added_sequences, 3)
        self.assertEqual(result.crfs_with_multiple_references, 1)
        self.assertEqual(len(list(self.references_dir.glob("*.fasta"))), 3)
        self.assertEqual(len(table), 3)
        self.assertEqual(table["accession"].nunique(), 3)

    def test_reference_strain_plus_additional_references_are_kept(self):
        self.write_table([])
        html = """
        <table>
          <tr>
            <td>CRF183_0107</td>
            <td>HLJ15116 (PV568686)</td>
            <td>CRF01, CRF07</td>
            <td>author</td>
            <td></td>
            <td></td>
          </tr>
          <tr><td colspan="6">
            CRF183_0107 was described using 3 near-complete genomes:
            HLJ15116 (PV568686), HLJ18014 (PV568687), and HLJ19012 (PV568688).
          </td></tr>
        </table>
        """
        crfs = self.load_crf_html(html)

        result = update_reference_dataset(
            base_dir=self.tmp_dir,
            email="test@example.com",
            update_lanl_alignment=False,
            crf_loader=lambda: crfs,
            genbank_fetcher=lambda accession: ("A" * 9000, None),
            sleep_seconds=0,
        )

        table = pd.read_csv(self.sequences_path, sep="\t", dtype=str, keep_default_na=False)
        self.assertEqual(result.added_sequences, 3)
        self.assertEqual(
            set(table["accession"]),
            {"PV568686", "PV568687", "PV568688"},
        )

    def test_crf_parser_includes_near_full_length_and_excludes_partial(self):
        html = """
        <table>
          <tr><td>CRF181_BC</td><td>20BS178 (PV781618)</td></tr>
          <tr><td>
            CRF181_BC was described based on 2 near-full-length genomes:
            20BS178 (PV781618) and 19BS006 (PV781619).
            A partial genome 19BS007 (PV781620) was also sequenced.
          </td></tr>
        </table>
        """

        crfs = self.load_crf_html(html)

        self.assertEqual(
            [reference.accession for reference in crfs["181_BC"]],
            ["PV781618", "PV781619"],
        )

    def test_crf_parser_keeps_valid_references_before_partial_in_same_sentence(self):
        html = """
        <table>
          <tr><td>CRF181_BC</td><td>20BS178 (PV781618)</td></tr>
          <tr><td>
            CRF181_BC was described using two near-full-length genomes
            PV781618 and PV781619 and one partial genome PV781620.
          </td></tr>
        </table>
        """

        crfs = self.load_crf_html(html)

        self.assertEqual(
            [reference.accession for reference in crfs["181_BC"]],
            ["PV781618", "PV781619"],
        )

    def test_crf_parser_keeps_near_full_length_and_excludes_partial_accession(self):
        html = """
        <table>
          <tr><td>CRF149_01B</td><td>sample (MW110767)</td></tr>
          <tr><td>
            CRF149_01B was described using two near-full-length genomes
            MW110767 and MW110768 and one partial genomic sequence JQ028315.
          </td></tr>
        </table>
        """

        crfs = self.load_crf_html(html)

        self.assertEqual(
            [reference.accession for reference in crfs["149_01B"]],
            ["MW110767", "MW110768"],
        )

    def test_crf_parser_excludes_pol_sequence_reference(self):
        html = """
        <table>
          <tr><td>CRF146_BC</td><td>pol sample (OR260538)</td></tr>
          <tr><td>
            OR260538 is a pol sequence. CRF146_BC was described using
            three near-complete genomes OR260539, OR260540 and OR260541.
          </td></tr>
        </table>
        """

        crfs = self.load_crf_html(html)

        self.assertEqual(
            [reference.accession for reference in crfs["146_BC"]],
            ["OR260539", "OR260540", "OR260541"],
        )

    def test_crf_parser_recognizes_near_complete_sequences_with_word_count(self):
        html = """
        <table>
          <tr><td>CRF190_0708</td><td>sample (OP100001)</td></tr>
          <tr><td>
            CRF190_0708 was described using eight near-complete sequences
            OP100001, OP100002, OP100003, OP100004, OP100005, OP100006,
            OP100007 and OP100008.
          </td></tr>
        </table>
        """

        crfs = self.load_crf_html(html)

        self.assertEqual(len(crfs["190_0708"]), 8)
        self.assertEqual(crfs["190_0708"][0].accession, "OP100001")
        self.assertEqual(crfs["190_0708"][-1].accession, "OP100008")

    def test_crf_parser_expands_accession_ranges_in_valid_context(self):
        html = """
        <table>
          <tr><td>CRF173_63A6</td><td>sample (PQ523366)</td></tr>
          <tr><td>
            CRF173_63A6 was described using 12 complete genomes PQ523366-PQ523377.
          </td></tr>
        </table>
        """

        crfs = self.load_crf_html(html)

        accessions = [reference.accession for reference in crfs["173_63A6"]]
        self.assertEqual(len(accessions), 12)
        self.assertEqual(accessions[0], "PQ523366")
        self.assertEqual(accessions[-1], "PQ523377")

    def test_crf_parser_expands_requested_regression_ranges(self):
        html = """
        <table>
          <tr><td>CRF163_0107</td><td>sample (PQ218317)</td></tr>
          <tr><td>CRF163_0107 was defined using three complete genomes PQ218317-PQ218319.</td></tr>
          <tr><td>CRF168_0107</td><td>sample (PP975254)</td></tr>
          <tr><td>CRF168_0107 was described using four full-length genomes PP975254-PP975257.</td></tr>
          <tr><td>CRF190_0708</td><td>sample (PV745739)</td></tr>
          <tr><td>CRF190_0708 was described using eight near-complete genomes PV745739-PV745746.</td></tr>
          <tr><td>CRF191_0708</td><td>sample (PV745747)</td></tr>
          <tr><td>CRF191_0708 was described using four near-full-length genomes PV745747-PV745750.</td></tr>
          <tr><td>CRF193_cpx</td><td>sample (PV743138)</td></tr>
          <tr><td>CRF193_cpx was defined by three complete genomes PV743138-PV743140.</td></tr>
          <tr><td>CRF197_0107</td><td>sample (PX492343)</td></tr>
          <tr><td>CRF197_0107 was based on three near-complete genome sequences PX492343-PX492345.</td></tr>
        </table>
        """

        crfs = self.load_crf_html(html)

        self.assertEqual(
            [reference.accession for reference in crfs["163_0107"]],
            ["PQ218317", "PQ218318", "PQ218319"],
        )
        self.assertEqual(
            [reference.accession for reference in crfs["168_0107"]],
            ["PP975254", "PP975255", "PP975256", "PP975257"],
        )
        self.assertEqual(
            [reference.accession for reference in crfs["190_0708"]],
            [
                "PV745739",
                "PV745740",
                "PV745741",
                "PV745742",
                "PV745743",
                "PV745744",
                "PV745745",
                "PV745746",
            ],
        )
        self.assertEqual(
            [reference.accession for reference in crfs["191_0708"]],
            ["PV745747", "PV745748", "PV745749", "PV745750"],
        )
        self.assertEqual(
            [reference.accession for reference in crfs["193_cpx"]],
            ["PV743138", "PV743139", "PV743140"],
        )
        self.assertEqual(
            [reference.accession for reference in crfs["197_0107"]],
            ["PX492343", "PX492344", "PX492345"],
        )

    def test_crf_parser_finds_multiple_genomes_without_using_reference_strain_only(self):
        html = """
        <table>
          <tr><td>CRF149_01B</td><td>sample1 (MW110767)</td></tr>
          <tr><td>CRF149_01B was defined using two near-full-length genomes MW110767 and MW110768.</td></tr>
          <tr><td>CRF170_0107</td><td>old_ref (OR514934)</td></tr>
          <tr><td>CRF170_0107 was described using three complete genomes OR514934-OR514936.</td></tr>
          <tr><td>CRF205_0107</td><td>20HZ2132 (PX262425)</td></tr>
          <tr><td>CRF205_0107 was defined by four near-complete genomes PX262425-PX262428 and one near-full-length genome OP921950.</td></tr>
        </table>
        """

        crfs = self.load_crf_html(html)

        self.assertEqual(
            [reference.accession for reference in crfs["149_01B"]],
            ["MW110767", "MW110768"],
        )
        self.assertEqual(
            [reference.accession for reference in crfs["170_0107"]],
            ["OR514934", "OR514935", "OR514936"],
        )
        self.assertEqual(
            [reference.accession for reference in crfs["205_0107"]],
            ["PX262425", "PX262426", "PX262427", "PX262428", "OP921950"],
        )

    def test_crf_parser_uses_following_detail_row_without_full_label(self):
        html = """
        <table>
          <tr>
            <td>CRF170_0107</td>
            <td>YNM132 (OR514935)</td>
            <td>CRF01, CRF07</td>
            <td>Y. Li</td>
            <td>3</td>
          </tr>
          <tr>
            <td colspan="7">
              This circulating recombinant form is based on the near-complete genomes
              of 3 sequences 170_0107.CN.YNM132.OR514935,
              170_0107.CN.YNM46.OR514934 and 170_0107.CN.YNM144.OR514936.
            </td>
          </tr>
          <tr>
            <td>CRF197_0108</td>
            <td>23BN004 (PX492343)</td>
            <td>CRF01, CRF08</td>
            <td>M Chen</td>
            <td>12</td>
          </tr>
          <tr>
            <td colspan="7">
              CRF197_0108 was described based on the near complete genomes of
              3 isolates; 197_0108.CN.2023.23BN004, 197_0108.CN.2023.23BN005
              and 197_0108.CN.2023.23BN051 (PX492343 to PX492345).
            </td>
          </tr>
        </table>
        """

        crfs = self.load_crf_html(html)

        self.assertEqual(
            [reference.accession for reference in crfs["170_0107"]],
            ["OR514935", "OR514934", "OR514936"],
        )
        self.assertEqual(
            [reference.accession for reference in crfs["197_0108"]],
            ["PX492343", "PX492344", "PX492345"],
        )

    def test_crf_parser_uses_reclassified_accessions_from_detail_narrative(self):
        html = """
        <table>
          <tr><td>CRF149_01B</td><td>BL1947-00 (MW110767)</td><td>CRF55, B</td></tr>
          <tr>
            <td colspan="7">
              Some sequences previously characterized as CRF55_01B
              (MW110767, MW110768 and JQ028315) had a different recombination
              pattern and should be classified as CRF149_01B.
            </td>
          </tr>
        </table>
        """

        crfs = self.load_crf_html(html)

        self.assertEqual(
            [reference.accession for reference in crfs["149_01B"]],
            ["MW110767", "MW110768", "JQ028315"],
        )

    def test_crf_parser_does_not_attach_hiv2_detail_to_previous_hiv1_crf(self):
        html = """
        <table>
          <tr><td>CRF205_0107</td><td>20HZ2132 (PX262425)</td><td>CRF01, CRF07</td></tr>
          <tr><td colspan="7">CRF205_0107 was described based on 5 genomes PX262425-PX262428 and OP921950.</td></tr>
          <tr><td colspan="7">The first HIV-2 CRF was defined based on 3 genomes and isolate 7312A (L36874).</td></tr>
        </table>
        """

        crfs = self.load_crf_html(html)

        self.assertEqual(
            [reference.accession for reference in crfs["205_0107"]],
            ["PX262425", "PX262426", "PX262427", "PX262428", "OP921950"],
        )

    def test_crf_parser_excludes_lanl_gene_only_sequence(self):
        html = """
        <table>
          <tr><td>CRF205_0107</td><td>sample (PX262425)</td></tr>
          <tr><td>
            CRF205_0107 was defined by four near-complete genomes PX262425-PX262428.
            OP000001 is a pol gene sequence.
          </td></tr>
        </table>
        """

        crfs = self.load_crf_html(html)

        self.assertEqual(
            [reference.accession for reference in crfs["205_0107"]],
            ["PX262425", "PX262426", "PX262427", "PX262428"],
        )

    def test_isolate_that_looks_like_accession_is_not_parsed_as_accession(self):
        self.assertEqual(
            accessions_from_text("151_0107.ZLQ04681.PP982938"),
            ("PP982938",),
        )
        self.assertEqual(
            accessions_from_text("172_0755.CN.2023.JH20233033.PQ356207"),
            ("PQ356207",),
        )

    def test_crf_parser_keeps_terminal_accession_from_lanl_dotted_name(self):
        html = """
        <table>
          <tr>
            <td>CRF172_0755</td>
            <td>JH20233033 (PQ356207)</td>
            <td>CRF07, CRF55</td>
          </tr>
          <tr>
            <td colspan="7">
              CRF172_0755 was described based on one near-complete genome
              172_0755.CN.2023.JH20233033.PQ356207.
            </td>
          </tr>
        </table>
        """

        crfs = self.load_crf_html(html)

        self.assertEqual(
            [reference.accession for reference in crfs["172_0755"]],
            ["PQ356207"],
        )

    def test_lanl_alignment_sync_reclassifies_existing_reference(self):
        self.write_fasta("AB123456-OLD_CRF.fasta", "AAAA")
        self.write_table([
            {
                "subtype": "OLD_CRF",
                "country": "PT",
                "year": "2024",
                "isolate": "sample",
                "accession": "AB123456",
                "align sequence": "AAAA",
                "sequence": "AAAA",
                "features": "{'gag': (1, 4)}",
                "group": "M",
            }
        ])

        result = update_reference_dataset(
            base_dir=self.tmp_dir,
            lanl_alignment_loader=lambda: LanlAlignment(
                year=2026,
                content="Ref.NEW_CRF.PT.2024.sample.AB123456,CCCC\n",
            ),
            crf_loader=lambda: {},
        )

        table = pd.read_csv(self.sequences_path, sep="\t", dtype=str, keep_default_na=False)
        self.assertEqual(result.lanl_fasta_files_added, 0)
        self.assertEqual(result.lanl_fasta_files_updated, 1)
        self.assertEqual(result.lanl_fasta_files_removed, 1)
        self.assertEqual(result.updated_rows, 1)
        self.assertFalse((self.references_dir / "AB123456-OLD_CRF.fasta").exists())
        self.assertTrue((self.references_dir / "AB123456-NEW_CRF.fasta").exists())
        self.assertEqual(table.loc[0, "subtype"], "NEW_CRF")
        self.assertEqual(table.loc[0, "sequence"], "CCCC")
        self.assertEqual(table.loc[0, "features"], "{'gag': (1, 4)}")

    def test_duplicate_tsv_rows_are_collapsed_to_current_lanl_classification(self):
        self.write_fasta("OR514935-170_0107.fasta", "CCCC")
        self.write_table([
            {
                "subtype": "151_0107",
                "country": "PT",
                "year": "2024",
                "isolate": "old",
                "accession": "OR514935",
                "align sequence": "AAAA",
                "sequence": "AAAA",
                "features": "None",
                "group": "M",
            },
            {
                "subtype": "170_0107",
                "country": "CN",
                "year": "2025",
                "isolate": "current",
                "accession": "OR514935",
                "align sequence": "CCCC",
                "sequence": "CCCC",
                "features": "{'gag': (1, 4)}",
                "group": "M",
            },
        ])

        result = update_reference_dataset(
            base_dir=self.tmp_dir,
            update_lanl_alignment=False,
            crf_loader=lambda: {
                "170_0107": CrfReference("170_0107", "CRF170_0107", "current", "OR514935")
            },
        )

        table = pd.read_csv(self.sequences_path, sep="\t", dtype=str, keep_default_na=False)
        self.assertEqual(result.duplicate_rows_removed, 1)
        self.assertEqual(len(table), 1)
        self.assertEqual(table.loc[0, "accession"], "OR514935")
        self.assertEqual(table.loc[0, "subtype"], "170_0107")
        self.assertEqual(table.loc[0, "isolate"], "current")
        self.assertEqual(table.loc[0, "features"], "{'gag': (1, 4)}")

    def test_crf_additional_reference_survives_second_alignment_sync(self):
        self.write_table([])
        alignment = LanlAlignment(
            year=2026,
            content="Ref.A1.PT.2024.sample.AB123456,AAAA\n",
        )
        crfs = {
            "205_0107": CrfReference("205_0107", "CRF205_0107", "20HZ2132", "PX262425")
        }

        first = update_reference_dataset(
            base_dir=self.tmp_dir,
            email="test@example.com",
            lanl_alignment_loader=lambda: alignment,
            crf_loader=lambda: crfs,
            genbank_fetcher=lambda accession: ("C" * 9000, None),
            sleep_seconds=0,
        )
        self.assertEqual(first.added_sequences, 1)
        self.assertTrue((self.references_dir / "PX262425-205_0107.fasta").exists())

        def fail_if_called(accession):
            raise AssertionError(f"unexpected fetch for {accession}")

        second = update_reference_dataset(
            base_dir=self.tmp_dir,
            email="test@example.com",
            lanl_alignment_loader=lambda: alignment,
            crf_loader=lambda: crfs,
            genbank_fetcher=fail_if_called,
            sleep_seconds=0,
        )

        table = pd.read_csv(self.sequences_path, sep="\t", dtype=str, keep_default_na=False)
        self.assertTrue((self.references_dir / "PX262425-205_0107.fasta").exists())
        self.assertEqual(set(table["accession"]), {"AB123456", "PX262425"})
        self.assertEqual(second.added_sequences, 0)
        self.assertEqual(second.added_rows, 0)
        self.assertEqual(second.updated_rows, 0)
        self.assertEqual(second.reclassified_references, 0)
        self.assertEqual(second.lanl_fasta_files_removed, 0)
        self.assertTrue(second.up_to_date)

    def test_update_is_idempotent_with_alignment_and_crf_union(self):
        self.write_table([])
        alignment = LanlAlignment(
            year=2026,
            content="Ref.A1.PT.2024.sample.AB123456,AAAA\n",
        )
        crfs = {
            "205_0107": CrfReference("205_0107", "CRF205_0107", "20HZ2132", "PX262425")
        }

        update_reference_dataset(
            base_dir=self.tmp_dir,
            email="test@example.com",
            lanl_alignment_loader=lambda: alignment,
            crf_loader=lambda: crfs,
            genbank_fetcher=lambda accession: ("C" * 9000, None),
            sleep_seconds=0,
        )
        result = update_reference_dataset(
            base_dir=self.tmp_dir,
            email="test@example.com",
            lanl_alignment_loader=lambda: alignment,
            crf_loader=lambda: crfs,
            genbank_fetcher=lambda accession: (_ for _ in ()).throw(
                AssertionError(f"unexpected fetch for {accession}")
            ),
            sleep_seconds=0,
        )

        self.assertEqual(result.added_sequences, 0)
        self.assertEqual(result.added_rows, 0)
        self.assertEqual(result.updated_rows, 0)
        self.assertEqual(result.reclassified_references, 0)
        self.assertEqual(result.lanl_fasta_files_added, 0)
        self.assertEqual(result.lanl_fasta_files_updated, 0)
        self.assertEqual(result.lanl_fasta_files_removed, 0)
        self.assertTrue(result.up_to_date)

    def test_lanl_alignment_refresh_updates_changed_fasta(self):
        self.write_fasta("AB123456-A1.fasta", "AAAA")
        self.write_table([])

        result = update_reference_dataset(
            base_dir=self.tmp_dir,
            lanl_alignment_loader=lambda: LanlAlignment(
                year=2025,
                content="Ref.A1.PT.2024.sample.AB123456,CCCC\n",
            ),
            crf_loader=lambda: {},
        )

        self.assertEqual(result.lanl_fasta_files_added, 0)
        self.assertEqual(result.lanl_fasta_files_updated, 1)
        self.assertIn(
            "CCCC",
            (self.references_dir / "AB123456-A1.fasta").read_text(encoding="utf-8"),
        )

    def test_lanl_alignment_refresh_ignores_fasta_line_wrapping(self):
        self.write_fasta("AB123456-A1.fasta", "ACGTACGT")
        self.write_table([])

        result = update_reference_dataset(
            base_dir=self.tmp_dir,
            lanl_alignment_loader=lambda: LanlAlignment(
                year=2025,
                content="Ref.A1.PT.2024.sample.AB123456,ACGTACGT\n",
            ),
            crf_loader=lambda: {},
        )

        self.assertEqual(result.lanl_fasta_files_added, 0)
        self.assertEqual(result.lanl_fasta_files_updated, 0)

    def test_crf_database_wins_alignment_conflict_with_warning(self):
        with self.assertWarnsRegex(RuntimeWarning, "CRF database label"):
            desired = build_desired_reference_map(
                [
                    LanlReferenceRecord(
                        accession="AB123456",
                        subtype="A1",
                        sequence="AAAA",
                    )
                ],
                (
                    CrfReference("B1", "CRFB1", "sample", "AB123456"),
                ),
            )

        self.assertEqual(desired["AB123456"].subtype, "B1")
        self.assertEqual(desired["AB123456"].source, "crf")
        self.assertEqual(desired["AB123456"].sequence, "AAAA")

    def test_current_crf_database_reclassifies_old_alignment_label(self):
        self.write_table([])
        alignment = LanlAlignment(
            year=2023,
            content=(
                "Ref.151_0107.CN.2022.sample.OR514934,AAAA\n"
                "Ref.151_0107.CN.2022.sample.OR514935,CCCC\n"
                "Ref.151_0107.CN.2022.sample.OR514936,GGGG\n"
            ),
        )
        crfs = {
            "170_0107": (
                CrfReference("170_0107", "CRF170_0107", "sample1", "OR514934"),
                CrfReference("170_0107", "CRF170_0107", "sample2", "OR514935"),
                CrfReference("170_0107", "CRF170_0107", "sample3", "OR514936"),
            )
        }

        result = update_reference_dataset(
            base_dir=self.tmp_dir,
            lanl_alignment_loader=lambda: alignment,
            crf_loader=lambda: crfs,
        )

        table = pd.read_csv(self.sequences_path, sep="\t", dtype=str, keep_default_na=False)
        self.assertEqual(result.reclassified_references, 0)
        self.assertEqual(set(table["subtype"]), {"170_0107"})
        self.assertTrue((self.references_dir / "OR514934-170_0107.fasta").exists())
        self.assertTrue((self.references_dir / "OR514935-170_0107.fasta").exists())
        self.assertTrue((self.references_dir / "OR514936-170_0107.fasta").exists())

    def test_current_crf_database_corrects_requested_historical_alignment_labels(self):
        cases = (
            (
                ["MW110767", "MW110768"],
                "55_01B",
                "149_01B",
            ),
            (
                ["OR514934", "OR514935", "OR514936"],
                "151_0107",
                "170_0107",
            ),
            (
                ["AY371159", "GQ229529"],
                "01_AE",
                "22_01A1",
            ),
            (
                [f"EU1701{suffix:02d}" for suffix in range(35, 56)],
                "42_BF",
                "42_BF1",
            ),
        )

        for accessions, old_label, current_label in cases:
            with self.subTest(current_label=current_label):
                desired = build_desired_reference_map(
                    [
                        LanlReferenceRecord(
                            accession=accession,
                            subtype=old_label,
                            sequence="A" * 9000,
                        )
                        for accession in accessions
                    ],
                    tuple(
                        CrfReference(
                            current_label,
                            f"CRF{current_label}",
                            "sample",
                            accession,
                        )
                        for accession in accessions
                    ),
                )

                self.assertEqual(set(desired), set(accessions))
                self.assertEqual(
                    {reference.subtype for reference in desired.values()},
                    {current_label},
                )
                self.assertEqual(
                    {reference.source for reference in desired.values()},
                    {"crf"},
                )

    def test_same_source_accession_label_conflict_is_error(self):
        with self.assertRaisesRegex(ReconciliationError, "Conflicting crf labels"):
            build_desired_reference_map(
                [],
                (
                    CrfReference("B1", "CRFB1", "sample1", "AB123456"),
                    CrfReference("C1", "CRFC1", "sample2", "AB123456"),
                ),
            )

    def test_same_crf_number_label_alias_warns_and_uses_canonical_label(self):
        with self.assertWarnsRegex(RuntimeWarning, "CRF database aliases"):
            desired = build_desired_reference_map(
                [],
                (
                    CrfReference("42_BF1", "CRF42_BF1", "sample1", "EU170144"),
                    CrfReference("42_BF", "CRF42_BF", "sample2", "EU170144"),
                ),
            )

        self.assertEqual(desired["EU170144"].subtype, "42_BF1")
        self.assertEqual(desired["EU170144"].source, "crf")

    def test_dry_run_reports_same_reclassification_as_real_run(self):
        self.write_fasta("AB123456-OLD_CRF.fasta", "AAAA")
        self.write_table([
            {
                "subtype": "OLD_CRF",
                "country": "PT",
                "year": "2024",
                "isolate": "sample",
                "accession": "AB123456",
                "align sequence": "AAAA",
                "sequence": "AAAA",
                "features": "None",
                "group": "M",
            }
        ])
        alignment = LanlAlignment(
            year=2026,
            content="Ref.NEW_CRF.PT.2024.sample.AB123456,CCCC\n",
        )

        dry_run = update_reference_dataset(
            base_dir=self.tmp_dir,
            dry_run=True,
            lanl_alignment_loader=lambda: alignment,
            crf_loader=lambda: {},
        )

        self.assertTrue((self.references_dir / "AB123456-OLD_CRF.fasta").exists())
        self.assertFalse((self.references_dir / "AB123456-NEW_CRF.fasta").exists())

        real = update_reference_dataset(
            base_dir=self.tmp_dir,
            lanl_alignment_loader=lambda: alignment,
            crf_loader=lambda: {},
        )

        self.assertEqual(dry_run.lanl_fasta_files_added, real.lanl_fasta_files_added)
        self.assertEqual(dry_run.lanl_fasta_files_updated, real.lanl_fasta_files_updated)
        self.assertEqual(dry_run.lanl_fasta_files_removed, real.lanl_fasta_files_removed)
        self.assertEqual(dry_run.reclassified_references, real.reclassified_references)
        self.assertEqual(dry_run.updated_rows, real.updated_rows)
        self.assertEqual(dry_run.total_references, real.total_references)

    def test_divergent_local_sequences_fetch_authoritative_accession(self):
        self.write_fasta("AB123456-A1.fasta", "AAAA")
        self.write_fasta("AB123456-OLD_CRF.fasta", "CCCC")
        self.write_table([])
        fetched = []

        def fetcher(accession):
            fetched.append(accession)
            return "GGGG", {"gag": (1, 4)}

        result = update_reference_dataset(
            base_dir=self.tmp_dir,
            update_lanl_alignment=False,
            crf_loader=lambda: {
                "A1": CrfReference("A1", "CRFA1", "sample", "AB123456")
            },
            genbank_fetcher=fetcher,
            sleep_seconds=0,
        )

        table = pd.read_csv(self.sequences_path, sep="\t", dtype=str, keep_default_na=False)
        self.assertEqual(fetched, ["AB123456"])
        self.assertEqual(result.lanl_fasta_files_updated, 1)
        self.assertEqual(result.lanl_fasta_files_removed, 1)
        self.assertTrue((self.references_dir / "AB123456-A1.fasta").exists())
        self.assertFalse((self.references_dir / "AB123456-OLD_CRF.fasta").exists())
        self.assertEqual(table.loc[0, "sequence"], "GGGG")
        self.assertEqual(table.loc[0, "features"], "{'gag': (1, 4)}")

    def test_crf_only_reference_survives_standalone_alignment_update(self):
        self.write_fasta("PX262425-205_0107.fasta", "CCCC")

        result = update_lanl_reference_alignment(
            self.tmp_dir,
            self.references_dir,
            alignment_loader=lambda: LanlAlignment(
                year=2026,
                content="Ref.A1.PT.2024.sample.AB123456,AAAA\n",
            ),
            crf_loader=lambda: {
                "205_0107": CrfReference("205_0107", "CRF205_0107", "20HZ2132", "PX262425")
            },
        )

        self.assertEqual(result.fasta_files_added, 1)
        self.assertEqual(result.fasta_files_removed, 0)
        self.assertTrue((self.references_dir / "AB123456-A1.fasta").exists())
        self.assertTrue((self.references_dir / "PX262425-205_0107.fasta").exists())

    def test_two_real_runs_leave_no_duplicate_accessions_and_are_idempotent(self):
        self.write_table([])
        alignment = LanlAlignment(
            year=2026,
            content=(
                "Ref.A1.PT.2024.sample.AB123456,AAAA\n"
                "Ref.B1.PT.2024.sample.CD123456,CCCC\n"
            ),
        )
        crfs = {
            "205_0107": CrfReference("205_0107", "CRF205_0107", "20HZ2132", "PX262425")
        }

        update_reference_dataset(
            base_dir=self.tmp_dir,
            email="test@example.com",
            lanl_alignment_loader=lambda: alignment,
            crf_loader=lambda: crfs,
            genbank_fetcher=lambda accession: ("G" * 9000, None),
            sleep_seconds=0,
        )
        result = update_reference_dataset(
            base_dir=self.tmp_dir,
            email="test@example.com",
            lanl_alignment_loader=lambda: alignment,
            crf_loader=lambda: crfs,
            genbank_fetcher=lambda accession: (_ for _ in ()).throw(
                AssertionError(f"unexpected fetch for {accession}")
            ),
            sleep_seconds=0,
        )

        table = pd.read_csv(self.sequences_path, sep="\t", dtype=str, keep_default_na=False)
        fasta_accessions = [
            parse_reference_name(path.stem)["accession"]
            for path in self.references_dir.glob("*.fasta")
        ]
        self.assertEqual(len(fasta_accessions), len(set(fasta_accessions)))
        self.assertEqual(len(table), len(fasta_accessions))
        self.assertEqual(table["accession"].nunique(), len(fasta_accessions))
        self.assertEqual(result.lanl_fasta_files_added, 0)
        self.assertEqual(result.lanl_fasta_files_updated, 0)
        self.assertEqual(result.reclassified_references, 0)
        self.assertEqual(result.lanl_fasta_files_removed, 0)
        self.assertTrue(result.up_to_date)

    def test_extract_lanl_alignment_content_from_pre_block(self):
        html = "<html><body><pre>Ref.A1.PT.2024.sample.AB123456,ACGT</pre></body></html>"

        content = extract_lanl_alignment_content(html)

        self.assertEqual(content, "Ref.A1.PT.2024.sample.AB123456,ACGT\n")

    def test_load_lanl_crf_references_parses_html_rows(self):
        html = """
        <table><tr><td>CRF205_0107</td><td>20HZ2132 (PX262425)</td></tr></table>
        """

        class FakeResponse:
            def __enter__(self):
                return self

            def __exit__(self, *args):
                return None

            def read(self):
                return html.encode("utf-8")

        with self.subTest("parser"):
            from unittest.mock import patch
            with patch("pyhiv.loading.reference_update.urlopen", return_value=FakeResponse()):
                parsed = load_lanl_crf_references()

        self.assertEqual(parsed["205_0107"][0].accession, "PX262425")
        self.assertEqual(parsed["205_0107"][0].strain, "20HZ2132")

    def test_load_lanl_crf_references_uses_parenthesized_accession(self):
        html = """
        <table><tr><td>CRF43_02G</td><td>J11223 (EU697904)</td></tr></table>
        """

        class FakeResponse:
            def __enter__(self):
                return self

            def __exit__(self, *args):
                return None

            def read(self):
                return html.encode("utf-8")

        from unittest.mock import patch
        with patch("pyhiv.loading.reference_update.urlopen", return_value=FakeResponse()):
            parsed = load_lanl_crf_references()

        self.assertEqual(parsed["43_02G"][0].accession, "EU697904")
        self.assertEqual(parsed["43_02G"][0].strain, "J11223")

    def test_load_lanl_crf_references_parses_multiple_reference_strains(self):
        html = """
        <table><tr><td>CRF120_0107</td><td>HES_LS16846 (OK662596), HES_LS17858 (OK662597)</td></tr></table>
        """

        class FakeResponse:
            def __enter__(self):
                return self

            def __exit__(self, *args):
                return None

            def read(self):
                return html.encode("utf-8")

        from unittest.mock import patch
        with patch("pyhiv.loading.reference_update.urlopen", return_value=FakeResponse()):
            parsed = load_lanl_crf_references()

        self.assertEqual(
            [reference.accession for reference in parsed["120_0107"]],
            ["OK662596", "OK662597"],
        )

    def test_load_lanl_crf_references_keeps_first_table_reference(self):
        html = """
        <table>
            <tr><td>CRF151_0107</td><td>220565 (PP982937)</td></tr>
            <tr><td>CRF151_0107 was described based on 151_0107.220415.PP982936</td></tr>
        </table>
        """

        class FakeResponse:
            def __enter__(self):
                return self

            def __exit__(self, *args):
                return None

            def read(self):
                return html.encode("utf-8")

        from unittest.mock import patch
        with patch("pyhiv.loading.reference_update.urlopen", return_value=FakeResponse()):
            parsed = load_lanl_crf_references()

        self.assertEqual(parsed["151_0107"][0].accession, "PP982937")
        self.assertEqual(parsed["151_0107"][0].strain, "220565")
