from pathlib import Path
from unittest import TestCase
import shutil
from unittest.mock import patch

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pyhiv import (
    LOW_SCORE_MARGIN_WARNING,
    NO_SUBTYPING_LABEL,
    PyHIV,
    SEQUENCE_TOO_LONG_WARNING,
    SPLITTING_ALIGNMENT_PREFIX,
    SPLITTING_MODE_HXB2,
    SPLITTING_MODE_NONE,
    SPLITTING_MODE_SUBTYPE,
    SUPPORTED_REFERENCE_GROUPS,
    has_reference_features,
    normalize_splitting_mode,
    normalize_reference_groups,
    reference_features,
    reference_has_features,
    selected_reference_accessions,
    summarize_closest_subtypes,
    subtype_score_warning,
)
from tests import TEST_DIR

DATA_DIR = TEST_DIR / "data" / "fastas"
REFERENCE_BASE = TEST_DIR / "data" / "references"


class TestPyHIV(TestCase):

    def setUp(self):
        """Prepare output directory for results."""
        self.output_dir = TEST_DIR / "output"
        if self.output_dir.exists():
            shutil.rmtree(self.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.output_dir, ignore_errors=True)

    @patch.dict("os.environ", {"REFERENCE_GENOMES_DIR": str(REFERENCE_BASE)})
    def test_pyhiv_with_real_fastas_splitting(self):
        """Run PyHIV on real FASTAs with splitting enabled."""
        PyHIV(
            fastas_dir=str(DATA_DIR),
            subtyping=True,
            splitting=True,
            output_dir=str(self.output_dir),
            n_jobs=12,
            reporting=True,
            alignment_tool="PyFamsa",
        )

        # Check that best alignment files are created
        alignment_files = list(self.output_dir.glob("best_alignment_*.fasta"))
        self.assertTrue(len(alignment_files) > 0)

        # Check that final table exists and has correct columns
        table_file = self.output_dir / "final_table.tsv"
        self.assertTrue(table_file.exists())
        table = pd.read_csv(table_file, sep='\t')
        expected_cols = [
            'File Name', 'Sequence', 'Reference', 'Group', 'Subtype', 'Closest Subtypes',
            'Subtype Score Warning', 'Splitting Reference', 'Most Matching Gene Region', 'Present Gene Regions'
        ]
        self.assertListEqual(list(table.columns), expected_cols)
        self.assertTrue(len(table) > 0)

    @patch.dict("os.environ", {"REFERENCE_GENOMES_DIR": str(REFERENCE_BASE)})
    def test_pyhiv_with_real_fastas_no_splitting(self):
        """Run PyHIV on real FASTAs with splitting disabled."""
        PyHIV(
            fastas_dir=str(DATA_DIR),
            subtyping=True,
            splitting=False,
            output_dir=str(self.output_dir),
            n_jobs=1,
            alignment_tool="PyFamsa",
            show_progress=True,
        )

        table_file = self.output_dir / "final_table.tsv"
        self.assertTrue(table_file.exists())
        table = pd.read_csv(table_file, sep='\t')

        # Columns specific to splitting should be dropped
        self.assertListEqual(
            list(table.columns),
            ['File Name', 'Sequence', 'Reference', 'Group', 'Subtype', 'Closest Subtypes', 'Subtype Score Warning'],
        )

    @patch("pyhiv.align_with_references")
    @patch("pyhiv.read_input_fastas")
    @patch("pyhiv.validate_reference_paths")
    @patch("pyhiv.get_reference_paths")
    def test_subtyping_without_splitting_uses_all_reference_fastas(
        self,
        mock_paths,
        mock_validate,
        mock_read_fastas,
        mock_align,
    ):
        """When no splitting is requested, subtype against every FASTA in the reference folder."""
        mock_paths.return_value = {
            "SEQUENCES_WITH_LOCATION": REFERENCE_BASE / "sequences_with_locations.tsv",
            "REFERENCE_GENOMES_FASTAS_DIR": REFERENCE_BASE / "reference_fastas",
            "HXB2_GENOME_FASTA_DIR": REFERENCE_BASE / "HXB2_fasta",
        }
        query = SeqRecord(Seq("AAA"), id="query")
        query.annotations["source_file"] = "sample.fasta"
        mock_read_fastas.return_value = [query]
        mock_align.return_value = ("AAA", "AAA", "not_in_table-X", [(3, "not_in_table-X")])

        PyHIV(
            fastas_dir=str(DATA_DIR),
            subtyping=True,
            splitting=False,
            output_dir=str(self.output_dir),
            reporting=False,
        )

        self.assertIsNone(mock_align.call_args.kwargs["allowed_reference_accessions"])

    @patch("pyhiv.align_with_references")
    @patch("pyhiv.read_input_fastas")
    @patch("pyhiv.validate_reference_paths")
    @patch("pyhiv.get_reference_paths")
    def test_subtyping_without_splitting_can_filter_reference_groups(
        self,
        mock_paths,
        mock_validate,
        mock_read_fastas,
        mock_align,
    ):
        """Explicit reference groups should still filter subtyping without splitting."""
        mock_paths.return_value = {
            "SEQUENCES_WITH_LOCATION": REFERENCE_BASE / "sequences_with_locations.tsv",
            "REFERENCE_GENOMES_FASTAS_DIR": REFERENCE_BASE / "reference_fastas",
            "HXB2_GENOME_FASTA_DIR": REFERENCE_BASE / "HXB2_fasta",
        }
        query = SeqRecord(Seq("AAA"), id="query")
        query.annotations["source_file"] = "sample.fasta"
        mock_read_fastas.return_value = [query]
        mock_align.return_value = ("AAA", "AAA", "AB253421-A1", [(3, "AB253421-A1")])

        PyHIV(
            fastas_dir=str(DATA_DIR),
            subtyping=True,
            splitting=False,
            output_dir=str(self.output_dir),
            reporting=False,
            reference_groups=("M",),
        )

        self.assertEqual(
            mock_align.call_args.kwargs["allowed_reference_accessions"],
            selected_reference_accessions(
                pd.read_csv(REFERENCE_BASE / "sequences_with_locations.tsv", sep="\t"),
                ("M",),
            ),
        )

    @patch("pyhiv.align_with_references", return_value=None)
    @patch("pyhiv.read_input_fastas")
    @patch("pyhiv.validate_reference_paths")
    @patch("pyhiv.get_reference_paths")
    def test_subtyping_with_splitting_filters_to_listed_reference_fastas(
        self,
        mock_paths,
        mock_validate,
        mock_read_fastas,
        mock_align,
    ):
        """When splitting is requested, ignore reference FASTAs absent from sequences_with_locations."""
        mock_paths.return_value = {
            "SEQUENCES_WITH_LOCATION": REFERENCE_BASE / "sequences_with_locations.tsv",
            "REFERENCE_GENOMES_FASTAS_DIR": REFERENCE_BASE / "reference_fastas",
            "HXB2_GENOME_FASTA_DIR": REFERENCE_BASE / "HXB2_fasta",
        }
        mock_read_fastas.return_value = [SeqRecord(Seq("AAA"), id="query")]

        PyHIV(
            fastas_dir=str(DATA_DIR),
            subtyping=True,
            splitting=True,
            output_dir=str(self.output_dir),
            reporting=False,
            reference_groups=("M", "N", "O", "P"),
        )

        allowed = mock_align.call_args.kwargs["allowed_reference_accessions"]
        listed_accessions = set(
            pd.read_csv(REFERENCE_BASE / "sequences_with_locations.tsv", sep="\t")["accession"].astype(str)
        )
        self.assertEqual(allowed, listed_accessions)

    @patch("pyhiv.align_with_references")
    @patch("pyhiv.read_input_fastas")
    @patch("pyhiv.validate_reference_paths")
    @patch("pyhiv.get_reference_paths")
    def test_subtyping_with_splitting_falls_back_to_hxb2_when_reference_has_no_features(
        self,
        mock_paths,
        mock_validate,
        mock_read_fastas,
        mock_align,
    ):
        """Subtype against the best reference, but split against HXB2 if its features are missing."""
        sequences_with_locations = self.output_dir / "sequences_with_locations.tsv"
        pd.DataFrame([
            {
                "accession": "ACCNO",
                "group": "M",
                "subtype": "B",
                "features": "None",
            },
            {
                "accession": "K03455",
                "group": "M",
                "subtype": "B",
                "features": "{'gag': (1, 3)}",
            },
        ]).to_csv(sequences_with_locations, sep="\t", index=False)

        mock_paths.return_value = {
            "SEQUENCES_WITH_LOCATION": sequences_with_locations,
            "REFERENCE_GENOMES_FASTAS_DIR": REFERENCE_BASE / "reference_fastas",
            "HXB2_GENOME_FASTA_DIR": REFERENCE_BASE / "HXB2_fasta",
        }
        query = SeqRecord(Seq("AAA"), id="query")
        query.annotations["source_file"] = "sample.fasta"
        mock_read_fastas.return_value = [query]
        mock_align.side_effect = [
            ("AAA", "AAA", "ACCNO-B", [(5, "ACCNO-B"), (4, "K03455-B")]),
            ("AAA", "AAA", "K03455-B"),
        ]

        PyHIV(
            fastas_dir=str(DATA_DIR),
            subtyping=True,
            splitting=True,
            output_dir=str(self.output_dir),
            reporting=False,
        )

        self.assertEqual(mock_align.call_count, 2)
        self.assertEqual(
            mock_align.call_args_list[0].kwargs["references_dir"],
            REFERENCE_BASE / "reference_fastas",
        )
        self.assertEqual(
            mock_align.call_args_list[0].kwargs["allowed_reference_accessions"],
            {"ACCNO", "K03455"},
        )
        self.assertEqual(
            mock_align.call_args_list[1].kwargs["references_dir"],
            REFERENCE_BASE / "HXB2_fasta",
        )

        table = pd.read_csv(self.output_dir / "final_table.tsv", sep='\t')
        self.assertEqual(table.loc[0, "File Name"], "sample.fasta")
        self.assertEqual(table.loc[0, "Reference"], "ACCNO")
        self.assertEqual(table.loc[0, "Splitting Reference"], "K03455")
        self.assertEqual(table.loc[0, "Group"], "M")
        self.assertEqual(table.loc[0, "Subtype"], "B")
        self.assertTrue((self.output_dir / f"best_alignment_sample_query.fasta").exists())
        self.assertTrue((self.output_dir / f"{SPLITTING_ALIGNMENT_PREFIX}_sample_query.fasta").exists())

    @patch("pyhiv.align_with_references")
    @patch("pyhiv.read_input_fastas")
    @patch("pyhiv.validate_reference_paths")
    @patch("pyhiv.get_reference_paths")
    def test_splitting_disambiguates_features_that_share_an_output_file(
        self,
        mock_paths,
        mock_validate,
        mock_read_fastas,
        mock_align,
    ):
        """Two annotated features that map to the same output file must not overwrite each other."""
        sequences_with_locations = self.output_dir / "sequences_with_locations.tsv"
        pd.DataFrame([
            {
                "accession": "ACCNO",
                "group": "M",
                "subtype": "B",
                # "gag-pol" and "gag-pol fusion" both map to the same
                # ("gag-pol",) output path/suffix in FEATURE_OUTPUT_PATHS.
                "features": "{'gag-pol': (1, 3), 'gag-pol fusion': (4, 6)}",
            },
        ]).to_csv(sequences_with_locations, sep="\t", index=False)

        mock_paths.return_value = {
            "SEQUENCES_WITH_LOCATION": sequences_with_locations,
            "REFERENCE_GENOMES_FASTAS_DIR": REFERENCE_BASE / "reference_fastas",
            "HXB2_GENOME_FASTA_DIR": REFERENCE_BASE / "HXB2_fasta",
        }
        query = SeqRecord(Seq("AAA"), id="query")
        query.annotations["source_file"] = "sample.fasta"
        mock_read_fastas.return_value = [query]
        mock_align.return_value = ("ACGTAC", "ACGTAC", "ACCNO-B", [(6, "ACCNO-B")])

        PyHIV(
            fastas_dir=str(DATA_DIR),
            subtyping=True,
            splitting=True,
            output_dir=str(self.output_dir),
            reporting=False,
        )

        gene_dir = self.output_dir / "gag-pol"
        first_file = gene_dir / "sample_query_gag-pol.fasta"
        second_file = gene_dir / "sample_query_gag-pol-fusion.fasta"

        self.assertTrue(first_file.exists())
        self.assertTrue(second_file.exists())
        self.assertIn("ACG", first_file.read_text())
        self.assertIn("TAC", second_file.read_text())

    @patch("pyhiv.align_with_references")
    @patch("pyhiv.read_input_fastas")
    @patch("pyhiv.validate_reference_paths")
    @patch("pyhiv.get_reference_paths")
    def test_hxb2_fallback_alignment_failure_records_placeholder_row(
        self,
        mock_paths,
        mock_validate,
        mock_read_fastas,
        mock_align,
    ):
        """When the HXB2 fallback alignment fails, the row keeps subtyping data but placeholders for splitting."""
        sequences_with_locations = self.output_dir / "sequences_with_locations.tsv"
        pd.DataFrame([
            {
                "accession": "ACCNO",
                "group": "M",
                "subtype": "B",
                "features": "None",
            },
            {
                "accession": "K03455",
                "group": "M",
                "subtype": "B",
                "features": "{'gag': (1, 3)}",
            },
        ]).to_csv(sequences_with_locations, sep="\t", index=False)

        mock_paths.return_value = {
            "SEQUENCES_WITH_LOCATION": sequences_with_locations,
            "REFERENCE_GENOMES_FASTAS_DIR": REFERENCE_BASE / "reference_fastas",
            "HXB2_GENOME_FASTA_DIR": REFERENCE_BASE / "HXB2_fasta",
        }
        query = SeqRecord(Seq("AAA"), id="query")
        query.annotations["source_file"] = "sample.fasta"
        mock_read_fastas.return_value = [query]
        mock_align.side_effect = [
            ("AAA", "AAA", "ACCNO-B", [(5, "ACCNO-B"), (4, "K03455-B")]),
            None,
        ]

        PyHIV(
            fastas_dir=str(DATA_DIR),
            subtyping=True,
            splitting=True,
            output_dir=str(self.output_dir),
            reporting=False,
        )

        self.assertEqual(mock_align.call_count, 2)

        table = pd.read_csv(self.output_dir / "final_table.tsv", sep='\t')
        self.assertEqual(len(table), 1)
        self.assertEqual(table.loc[0, "Reference"], "ACCNO")
        self.assertEqual(table.loc[0, "Group"], "M")
        self.assertEqual(table.loc[0, "Subtype"], "B")
        self.assertEqual(table.loc[0, "Splitting Reference"], "-")
        self.assertEqual(table.loc[0, "Most Matching Gene Region"], "-")
        self.assertEqual(table.loc[0, "Present Gene Regions"], "-")
        self.assertTrue((self.output_dir / "best_alignment_sample_query.fasta").exists())
        self.assertFalse((self.output_dir / f"{SPLITTING_ALIGNMENT_PREFIX}_sample_query.fasta").exists())

    @patch("pyhiv.reference_features", return_value={"gag": (1, 3)})
    @patch("pyhiv.align_with_references")
    @patch("pyhiv.read_input_fastas")
    @patch("pyhiv.validate_reference_paths")
    @patch("pyhiv.get_reference_paths")
    def test_subtyping_can_split_against_hxb2_reference(
        self,
        mock_paths,
        mock_validate,
        mock_read_fastas,
        mock_align,
        mock_reference_features,
    ):
        """Subtyping can keep the subtype call while using HXB2 for gene splitting."""
        mock_paths.return_value = {
            "SEQUENCES_WITH_LOCATION": REFERENCE_BASE / "sequences_with_locations.tsv",
            "REFERENCE_GENOMES_FASTAS_DIR": REFERENCE_BASE / "reference_fastas",
            "HXB2_GENOME_FASTA_DIR": REFERENCE_BASE / "HXB2_fasta",
        }
        mock_read_fastas.return_value = [SeqRecord(Seq("AAA"), id="query")]
        mock_align.side_effect = [
            ("AAA", "AAA", "AB253421-A1", [(3, "AB253421-A1")]),
            ("AAA", "AAA", "K03455-B",),
        ]

        PyHIV(
            fastas_dir=str(DATA_DIR),
            subtyping=True,
            splitting="hxb2",
            output_dir=str(self.output_dir),
            reporting=False,
        )

        self.assertEqual(mock_align.call_count, 2)
        self.assertEqual(
            mock_align.call_args_list[0].kwargs["references_dir"],
            REFERENCE_BASE / "reference_fastas",
        )
        self.assertEqual(
            mock_align.call_args_list[1].kwargs["references_dir"],
            REFERENCE_BASE / "HXB2_fasta",
        )
        mock_reference_features.assert_called_once()
        self.assertEqual(mock_reference_features.call_args.args[1], "K03455")

        table = pd.read_csv(self.output_dir / "final_table.tsv", sep='\t')
        self.assertEqual(table.loc[0, "Reference"], "AB253421")
        self.assertEqual(table.loc[0, "Splitting Reference"], "K03455")
        self.assertTrue((self.output_dir / f"best_alignment_query.fasta").exists())
        self.assertTrue((self.output_dir / f"{SPLITTING_ALIGNMENT_PREFIX}_query.fasta").exists())

    @patch("pyhiv.reference_features", return_value={"gag": (1, 3)})
    @patch("pyhiv.align_with_references")
    @patch("pyhiv.read_input_fastas")
    @patch("pyhiv.validate_reference_paths")
    @patch("pyhiv.get_reference_paths")
    def test_no_subtyping_forces_subtype_splitting_to_hxb2(
        self,
        mock_paths,
        mock_validate,
        mock_read_fastas,
        mock_align,
        mock_reference_features,
    ):
        """Without subtyping, active splitting always uses HXB2 coordinates."""
        mock_paths.return_value = {
            "SEQUENCES_WITH_LOCATION": REFERENCE_BASE / "sequences_with_locations.tsv",
            "REFERENCE_GENOMES_FASTAS_DIR": REFERENCE_BASE / "reference_fastas",
            "HXB2_GENOME_FASTA_DIR": REFERENCE_BASE / "HXB2_fasta",
        }
        mock_read_fastas.return_value = [SeqRecord(Seq("AAA"), id="query")]
        mock_align.return_value = ("AAA", "AAA", "K03455-B", [(3, "K03455-B")])

        PyHIV(
            fastas_dir=str(DATA_DIR),
            subtyping=False,
            splitting="subtype",
            output_dir=str(self.output_dir),
            reporting=False,
        )

        mock_align.assert_called_once()
        self.assertEqual(
            mock_align.call_args.kwargs["references_dir"],
            REFERENCE_BASE / "HXB2_fasta",
        )
        mock_reference_features.assert_called_once()
        self.assertEqual(mock_reference_features.call_args.args[1], "K03455")

        table = pd.read_csv(self.output_dir / "final_table.tsv", sep='\t')
        self.assertEqual(table.loc[0, "Reference"], "K03455")
        self.assertEqual(table.loc[0, "Splitting Reference"], "K03455")
        self.assertEqual(table.loc[0, "Group"], NO_SUBTYPING_LABEL)
        self.assertEqual(table.loc[0, "Subtype"], NO_SUBTYPING_LABEL)

    @patch.dict("os.environ", {"REFERENCE_GENOMES_DIR": str(REFERENCE_BASE)})
    def test_pyhiv_no_subtyping_uses_uniform_table_labels(self):
        """HXB2 alignment should not report reference group/subtype as subtype results."""
        PyHIV(
            fastas_dir=str(DATA_DIR),
            subtyping=False,
            splitting=True,
            output_dir=str(self.output_dir),
            n_jobs=1,
            reporting=False,
            alignment_tool="PyFamsa",
        )

        table = pd.read_csv(self.output_dir / "final_table.tsv", sep='\t')
        self.assertTrue(len(table) > 0)
        self.assertTrue((table["Group"] == NO_SUBTYPING_LABEL).all())
        self.assertTrue((table["Subtype"] == NO_SUBTYPING_LABEL).all())
        self.assertTrue((table["Closest Subtypes"] == NO_SUBTYPING_LABEL).all())

    @patch.dict("os.environ", {"REFERENCE_GENOMES_DIR": str(REFERENCE_BASE)})
    @patch("pyhiv.align_with_references", return_value=None)
    def test_best_alignment_is_none(self, mock_align):
        """Test that PyHIV handles cases when align_with_references returns None."""
        # Run PyHIV on a real FASTA directory, alignment will always return None
        PyHIV(
            fastas_dir=str(DATA_DIR),
            subtyping=True,
            splitting=True,
            output_dir=str(self.output_dir),
            n_jobs=1
        )

        # The output directory should exist but no best_alignment files should be created
        alignment_files = list(Path(self.output_dir).glob("best_alignment_*.fasta"))
        self.assertEqual(len(alignment_files), 0)

        # The final_table.tsv should still be created but empty (or just headers)
        table_file = Path(self.output_dir) / "final_table.tsv"
        self.assertTrue(table_file.exists())
        table = pd.read_csv(table_file, sep='\t')
        self.assertEqual(len(table), 0)  # No rows because no alignments succeeded

    @patch("pyhiv.align_with_references")
    @patch("pyhiv.read_input_fastas")
    @patch("pyhiv.validate_reference_paths")
    @patch("pyhiv.get_reference_paths")
    @patch("logging.warning")
    def test_sequence_longer_than_hiv_genome_is_skipped(
        self,
        mock_warning,
        mock_paths,
        mock_validate,
        mock_read_fastas,
        mock_align,
    ):
        """Sequences with more than 12000 nt should be skipped before alignment."""
        mock_paths.return_value = {
            "SEQUENCES_WITH_LOCATION": REFERENCE_BASE / "sequences_with_locations.tsv",
            "REFERENCE_GENOMES_FASTAS_DIR": REFERENCE_BASE / "reference_fastas",
            "HXB2_GENOME_FASTA_DIR": REFERENCE_BASE / "HXB2_fasta",
        }
        mock_read_fastas.return_value = [SeqRecord(Seq("A" * 12001), id="too_long")]

        PyHIV(
            fastas_dir=str(DATA_DIR),
            output_dir=str(self.output_dir),
            reporting=False,
        )

        mock_validate.assert_called_once()
        mock_align.assert_not_called()
        warning_calls = [call_args[0] for call_args in mock_warning.call_args_list]
        self.assertIn(
            ("%s Skipping sequence '%s'.", SEQUENCE_TOO_LONG_WARNING, "too_long"),
            warning_calls,
        )

        table_file = Path(self.output_dir) / "final_table.tsv"
        self.assertTrue(table_file.exists())
        table = pd.read_csv(table_file, sep='\t')
        self.assertEqual(len(table), 0)

    def test_normalize_reference_groups_accepts_comma_separated_groups(self):
        self.assertEqual(normalize_reference_groups("M,N,O,P"), ("M", "N", "O", "P"))

    def test_selected_reference_accessions_filters_by_group(self):
        reference_sequences = pd.DataFrame([
            {"accession": "acc_m", "group": "M"},
            {"accession": "acc_o", "group": "O"},
        ])

        self.assertEqual(
            selected_reference_accessions(reference_sequences, ("O",)),
            {"acc_o"},
        )

    def test_selected_reference_accessions_without_group_returns_all_listed(self):
        reference_sequences = pd.DataFrame([
            {"accession": "acc_1"},
            {"accession": "acc_2"},
        ])

        self.assertEqual(
            selected_reference_accessions(reference_sequences, ("M",)),
            {"acc_1", "acc_2"},
        )

    def test_selected_reference_accessions_can_require_features(self):
        reference_sequences = pd.DataFrame([
            {"accession": "acc_with_features", "group": "M", "features": "{'gag': (1, 10)}"},
            {"accession": "acc_without_features", "group": "M", "features": "None"},
            {"accession": "acc_other_group", "group": "N", "features": "{'env': (20, 30)}"},
        ])

        self.assertEqual(
            selected_reference_accessions(reference_sequences, ("M",), require_features=True),
            {"acc_with_features"},
        )

    def test_reference_has_features_handles_none_values(self):
        self.assertFalse(reference_has_features(None))
        self.assertFalse(reference_has_features("None"))
        self.assertFalse(reference_has_features("{}"))
        self.assertTrue(reference_has_features("{'gag': (1, 10)}"))

    def test_reference_has_features_accepts_dict_input(self):
        self.assertFalse(reference_has_features({}))
        self.assertTrue(reference_has_features({"gag": (1, 10)}))

    def test_reference_has_features_treats_unparsable_text_as_present(self):
        self.assertTrue(reference_has_features("not a python literal("))

    def test_reference_has_features_handles_none_literal_result(self):
        self.assertFalse(reference_has_features("(None)"))

    def test_reference_has_features_treats_non_dict_literal_as_present(self):
        self.assertTrue(reference_has_features("['gag', 'pol']"))

    def test_reference_features_raises_for_unknown_accession(self):
        reference_sequences = pd.DataFrame([
            {"accession": "acc1", "features": "{'gag': (1, 10)}"},
        ])
        with self.assertRaises(ValueError):
            reference_features(reference_sequences, "missing_accession")

    def test_has_reference_features_returns_false_without_features_column(self):
        reference_sequences = pd.DataFrame([{"accession": "acc1"}])
        self.assertFalse(has_reference_features(reference_sequences, "acc1"))

    def test_normalize_reference_groups_rejects_empty_group_list(self):
        with self.assertRaises(ValueError):
            normalize_reference_groups(",, ;")

    def test_normalize_reference_groups_all_alias_returns_supported_groups(self):
        self.assertEqual(normalize_reference_groups("all"), SUPPORTED_REFERENCE_GROUPS)

    def test_selected_reference_accessions_requires_features_column(self):
        reference_sequences = pd.DataFrame([{"accession": "acc1", "group": "M"}])
        with self.assertRaises(ValueError):
            selected_reference_accessions(reference_sequences, ("M",), require_features=True)

    def test_selected_reference_accessions_raises_when_none_match(self):
        reference_sequences = pd.DataFrame([{"accession": "acc1", "group": "N"}])
        with self.assertRaises(ValueError):
            selected_reference_accessions(reference_sequences, ("M",))

    def test_has_reference_features_looks_up_accession_features(self):
        reference_sequences = pd.DataFrame([
            {"accession": "acc_with_features", "features": "{'gag': (1, 10)}"},
            {"accession": "acc_without_features", "features": "None"},
        ])

        self.assertTrue(has_reference_features(reference_sequences, "acc_with_features"))
        self.assertFalse(has_reference_features(reference_sequences, "acc_without_features"))
        self.assertFalse(has_reference_features(reference_sequences, "missing"))

    def test_normalize_splitting_mode_rejects_unknown_mode(self):
        with self.assertRaises(ValueError):
            normalize_splitting_mode("not-a-real-mode")

    def test_normalize_reference_groups_rejects_unsupported_group(self):
        with self.assertRaises(ValueError):
            normalize_reference_groups("Z")

    def test_normalize_splitting_mode_accepts_new_modes(self):
        self.assertEqual(normalize_splitting_mode(True, subtyping=True), SPLITTING_MODE_SUBTYPE)
        self.assertEqual(normalize_splitting_mode(False, subtyping=True), SPLITTING_MODE_NONE)
        self.assertEqual(normalize_splitting_mode("reference", subtyping=True), SPLITTING_MODE_HXB2)
        self.assertEqual(normalize_splitting_mode("hxb2", subtyping=True), SPLITTING_MODE_HXB2)
        self.assertEqual(normalize_splitting_mode("subtype", subtyping=True), SPLITTING_MODE_SUBTYPE)
        self.assertEqual(normalize_splitting_mode("subtype", subtyping=False), SPLITTING_MODE_HXB2)

    def test_summarize_closest_subtypes_returns_top_unique_subtypes(self):
        metadata_by_accession = {
            "acc_b1": {"group": "M", "subtype": "B"},
            "acc_b2": {"group": "M", "subtype": "B"},
            "acc_c": {"group": "M", "subtype": "C"},
            "acc_o": {"group": "O", "subtype": "O"},
        }

        result = summarize_closest_subtypes(
            [
                (100, "acc_b1-B"),
                (98, "acc_b2-B"),
                (95, "acc_c-C"),
                (90, "acc_o-O"),
            ],
            metadata_by_accession,
            top_n=3,
        )

        self.assertEqual(result, "M:B (score=100); M:C (score=95); O:O (score=90)")

    def test_subtype_score_warning_flags_low_relative_margin(self):
        metadata_by_accession = {
            "acc_b": {"group": "M", "subtype": "B"},
            "acc_c": {"group": "M", "subtype": "C"},
        }

        self.assertEqual(
            subtype_score_warning(
                [(100, "acc_b-B"), (99, "acc_c-C")],
                metadata_by_accession,
            ),
            LOW_SCORE_MARGIN_WARNING,
        )
        self.assertEqual(
            subtype_score_warning(
                [(100, "acc_b-B"), (98, "acc_c-C")],
                metadata_by_accession,
            ),
            "",
        )

    def test_subtype_score_warning_returns_empty_when_top_score_not_positive(self):
        metadata_by_accession = {
            "acc_b": {"group": "M", "subtype": "B"},
            "acc_c": {"group": "M", "subtype": "C"},
        }

        self.assertEqual(
            subtype_score_warning(
                [(0, "acc_b-B"), (0, "acc_c-C")],
                metadata_by_accession,
            ),
            "",
        )
