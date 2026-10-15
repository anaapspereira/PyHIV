import csv
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import patch

from pyhiv import write_other_region_mutation_files
from pyhiv.mutations import (
    call_other_region_mutations,
    hxb2_reference_nt_for_other_region,
    load_other_region_resistance_catalog,
    write_other_region_mutation_matrices_tsv,
    write_other_region_mutations_tsv,
    write_other_region_resistance_tsv,
    write_other_region_warnings_tsv,
)


class TestOtherRegionMutationCaller(TestCase):
    def setUp(self):
        load_other_region_resistance_catalog.cache_clear()

    def tearDown(self):
        load_other_region_resistance_catalog.cache_clear()

    def vif_m1v_nt(self):
        sequence = list(hxb2_reference_nt_for_other_region("VIF"))
        sequence[0] = "G"
        return "".join(sequence)

    def read_tsv(self, path):
        with open(path, newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_calls_other_region_mutation_against_hxb2_without_resistance_marker(self):
        calls = call_other_region_mutations(self.vif_m1v_nt(), "sample1", "vif")

        self.assertEqual(len(calls), 1)
        call = calls[0]
        self.assertEqual(call.region, "VIF")
        self.assertEqual(call.position, 1)
        self.assertEqual(call.mutation, "VIF:M1V")
        self.assertEqual(call.reference_system, "HXB2")
        self.assertFalse(call.resistance_marker)
        self.assertEqual(call.resistance_evaluation_status, "NOT_EVALUATED_NO_CATALOG_MARKERS")
        self.assertIn("mutations are reported against HXB2 only", call.resistance_evaluation_note)
        self.assertEqual(call.drug_class, "")
        self.assertEqual(call.drug, "")

    def test_external_catalog_marks_other_region_resistance(self):
        with TemporaryDirectory() as directory:
            catalog = Path(directory) / "other_region_resistance.tsv"
            catalog.write_text(
                "region\tposition\tref_aa\talt_aa\tdrug_class\tdrug\tresistance_scope\tevidence_level\tsource\tsource_url\tnote\n"
                "VIF\t1\tM\tV\tAPOBEC antagonist\tExampleDrug\tcandidate\tliterature\tTest source\thttps://example.test\tTest note\n"
            )
            with patch.dict(os.environ, {"PYHIV_OTHER_REGION_RESISTANCE_CATALOG": str(catalog)}):
                load_other_region_resistance_catalog.cache_clear()
                calls = call_other_region_mutations(self.vif_m1v_nt(), "sample1", "VIF")

        self.assertEqual(len(calls), 1)
        call = calls[0]
        self.assertTrue(call.resistance_marker)
        self.assertEqual(call.resistance_evaluation_status, "CATALOG_MARKER_MATCH")
        self.assertEqual(call.drug_class, "APOBEC antagonist")
        self.assertEqual(call.drug, "ExampleDrug")
        self.assertEqual(call.resistance_source, "Test source")

    def test_writes_separate_mutation_resistance_and_matrix_outputs(self):
        calls = call_other_region_mutations(self.vif_m1v_nt(), "sample1", "VIF")
        with TemporaryDirectory() as directory:
            output_dir = Path(directory)
            write_other_region_mutations_tsv(calls, output_dir / "other_region_mutations.tsv")
            write_other_region_resistance_tsv(calls, output_dir / "other_region_resistance.tsv")
            write_other_region_warnings_tsv(
                [("sample1", "VIF")],
                output_dir / "other_region_warnings.tsv",
            )
            matrix_paths = write_other_region_mutation_matrices_tsv(
                calls,
                output_dir,
                sequence_region_pairs=[("sample1", "VIF")],
            )

            mutation_rows = self.read_tsv(output_dir / "other_region_mutations.tsv")
            resistance_rows = self.read_tsv(output_dir / "other_region_resistance.tsv")
            matrix_rows = self.read_tsv(output_dir / "other_region_mutation_matrix_VIF.tsv")
            warning_rows = self.read_tsv(output_dir / "other_region_warnings.tsv")

        self.assertEqual(mutation_rows[0]["mutation"], "VIF:M1V")
        self.assertEqual(mutation_rows[0]["resistance_marker"], "False")
        self.assertEqual(mutation_rows[0]["resistance_evaluation_status"], "NOT_EVALUATED_NO_CATALOG_MARKERS")
        self.assertEqual(resistance_rows, [])
        self.assertEqual(warning_rows[0]["warning_type"], "NOT_EVALUATED_NO_CATALOG_MARKERS")
        self.assertEqual([path.name for path in matrix_paths], ["other_region_mutation_matrix_VIF.tsv"])
        self.assertEqual(matrix_rows[0]["Sequence"], "REFERENCE")
        self.assertEqual(matrix_rows[1]["Sequence"], "sample1")
        self.assertEqual(matrix_rows[1]["VIF:M1"], "M1V")

    def test_pipeline_helper_writes_other_region_files_without_final_summary_inputs(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            fasta_path = root / "VIF.fasta"
            fasta_path.write_text(f">sample1\n{self.vif_m1v_nt()}\n")
            calls = write_other_region_mutation_files([fasta_path], root / "out")

            self.assertEqual([call.mutation for call in calls], ["VIF:M1V"])
            self.assertTrue((root / "out" / "other_region_mutations.tsv").exists())
            self.assertTrue((root / "out" / "other_region_resistance.tsv").exists())
            self.assertTrue((root / "out" / "other_region_warnings.tsv").exists())
            self.assertTrue((root / "out" / "other_region_mutation_matrix_VIF.tsv").exists())

from pyhiv.mutations import infer_gene_targets_from_path, infer_other_region_targets_from_path
from pyhiv.mutations.caller import call_mutations_for_fasta_files_with_qc, hxb2_genome_sequence, sequence_gene_pairs_from_fasta_files


class TestBroadSplitRegionRouting(TestCase):
    def test_broad_pol_file_routes_to_pr_rt_and_in(self):
        self.assertEqual(
            infer_gene_targets_from_path(Path("results/pol/sample_pol.fasta")),
            ("PR", "RT", "IN"),
        )

    def test_broad_gag_pol_file_routes_to_gag_and_pol_targets(self):
        cases = [
            Path("results/gag-pol/sample_gag-pol.fasta"),
            Path("results/gag-pol.fasta"),
            Path("results/gag-pol-fusion/sample_gag-pol-fusion.fasta"),
            Path("results/gag_pol precursor/sample_gag_pol_precursor.fasta"),
        ]
        for path in cases:
            with self.subTest(path=path):
                self.assertEqual(
                    infer_gene_targets_from_path(path),
                    ("CA", "PR", "RT", "IN"),
                )
                self.assertEqual(
                    infer_other_region_targets_from_path(path),
                    ("MA", "NC", "P6"),
                )

    def test_specific_pol_child_is_not_treated_as_full_pol(self):
        path = Path("results/pol/p15-rnase-h/sample_p15-rnase-h.fasta")
        self.assertEqual(infer_gene_targets_from_path(path), ())
        self.assertEqual(infer_other_region_targets_from_path(path), ("RH",))

    def test_broad_gag_routes_ca_and_other_gag_regions_but_capsid_stays_specific(self):
        self.assertEqual(
            infer_gene_targets_from_path(Path("results/gag/sample_gag.fasta")),
            ("CA",),
        )
        self.assertEqual(
            infer_other_region_targets_from_path(Path("results/gag/sample_gag.fasta")),
            ("MA", "NC", "P6"),
        )
        self.assertEqual(
            infer_other_region_targets_from_path(Path("results/gag/p24-capsid/sample_p24-capsid.fasta")),
            (),
        )

    def test_broad_pol_fasta_is_called_for_all_consensus_pol_targets(self):
        pol_sequence = hxb2_genome_sequence()[2253 - 1:5093]
        with TemporaryDirectory() as directory:
            path = Path(directory) / "pol" / "sample_pol.fasta"
            path.parent.mkdir()
            path.write_text(f">sample1\n{pol_sequence}\n")

            calls, summaries, position_qc = call_mutations_for_fasta_files_with_qc([path], sequence_type="nt")
            pairs = sequence_gene_pairs_from_fasta_files([path])

        self.assertEqual(pairs, [("sample1", "PR"), ("sample1", "RT"), ("sample1", "IN")])
        self.assertEqual({summary.gene for summary in summaries}, {"PR", "RT", "IN"})
        self.assertEqual({row.gene for row in position_qc}, {"PR", "RT", "IN"})
        self.assertTrue({call.gene for call in calls}.issubset({"PR", "RT", "IN"}))
