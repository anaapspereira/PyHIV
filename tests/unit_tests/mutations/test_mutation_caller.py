import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import patch

from pyhiv.mutations import (
    MUTATION_INPUT_QC_TSV_COLUMNS,
    MUTATION_POSITION_QC_TSV_COLUMNS,
    call_mutations,
    call_mutations_with_qc,
    mutation_input_qc_for_fasta_files,
    validate_mutation_input_sequence,
    write_drm_screening_tsv,
    write_mutation_input_qc_tsv,
    write_mutation_matrices_tsv,
    write_mutation_position_qc_tsv,
    write_mutations_tsv,
)
from pyhiv.mutations.caller import (
    DEFAULT_NT_ALIGNMENT_BACKEND,
    DRM_SCREENING_NOT_COVERED,
    DRM_SCREENING_RESOLVED_MUTANT,
    DRM_SCREENING_RESOLVED_WT,
    DRM_SCREENING_UNRESOLVED,
    align_dna_to_hxb2_columns,
    align_nucleotides_semiglobal_codonaware,
    build_mutation_call,
    codon_calls_from_raw_nucleotide_alignment,
    coverage_confidence_not_covered_codons,
    drm_catalog_positions,
    drm_screening_status_by_sequence,
    hxb2_reference_nt_for_gene,
    infer_gene_from_path,
    internal_gap_score,
    mutations_from_alignment,
    mutation_position_qc_rows,
    mutations_from_hxb2_aligned_dna,
    mutations_from_reference_aligned_dna,
    protein_sequence_for_calling,
    summarize_drm_screening,
)
from pyhiv.mutations.references import reference_for_gene


class TestMutationCaller(TestCase):
    def semiglobal_primary_score(self, query_aligned, ref_aligned):
        score = 0
        index = 0
        alignment_length = len(query_aligned)
        while index < alignment_length:
            query_nt = query_aligned[index]
            ref_nt = ref_aligned[index]
            if query_nt != "-" and ref_nt != "-":
                score += 2 if query_nt == ref_nt else -2
                index += 1
                continue
            gap_start = index
            gap_side = "query" if query_nt == "-" else "reference"
            while index < alignment_length:
                same_query_gap = gap_side == "query" and query_aligned[index] == "-"
                same_ref_gap = gap_side == "reference" and ref_aligned[index] == "-"
                if not (same_query_gap or same_ref_gap):
                    break
                index += 1
            if gap_start == 0 or index == alignment_length:
                continue
            gap_length = index - gap_start
            score += -20 + -1 * (gap_length - 1)
        return score

    def gap_lengths(self, query_aligned, ref_aligned, gap_side):
        lengths = []
        index = 0
        alignment_length = len(query_aligned)
        while index < alignment_length:
            is_gap = query_aligned[index] == "-" if gap_side == "query" else ref_aligned[index] == "-"
            if not is_gap:
                index += 1
                continue
            start = index
            while index < alignment_length:
                is_gap = query_aligned[index] == "-" if gap_side == "query" else ref_aligned[index] == "-"
                if not is_gap:
                    break
                index += 1
            lengths.append(index - start)
        return lengths

    def drm_screening_codons(self, gene='RT'):
        return [
            {
                'alt_aa': aa,
                'codon_status': 'OBSERVED',
                'hxb2_ref_codon': '',
                'query_codon': '',
                'possible_alt_aas': aa,
                'inserted_nts': '',
                'deleted_nt_count': 0,
                'phase_at_start': 0,
                'phase_at_end': 0,
            }
            for aa in reference_for_gene(gene)
        ]

    def drm_screening_position_statuses(self, summary):
        return {position.position: position.status for position in summary.positions}

    def test_drm_screening_complete_without_drms(self):
        summary = summarize_drm_screening(
            sequence_id='complete_wt',
            gene='RT',
            codon_calls=self.drm_screening_codons('RT'),
            mutation_calls=[],
        )

        self.assertEqual(summary.status, 'COMPLETE')
        self.assertEqual(summary.gene, 'RT')
        self.assertTrue(summary.positions)
        self.assertTrue(all(
            position.status == DRM_SCREENING_RESOLVED_WT
            for position in summary.positions
        ))

    def test_drm_screening_complete_with_drm(self):
        codons = self.drm_screening_codons('RT')
        codons[102]['alt_aa'] = 'N'
        call = build_mutation_call(
            sequence_id='complete_drm', gene='RT', position=103, ref_aa='K', alt_aa='N',
            mutation_type='substitution', qc_status='PASS',
            reference_system='HIVDB_CONSENSUS_B',
        )

        summary = summarize_drm_screening(
            sequence_id='complete_drm',
            gene='RT',
            codon_calls=codons,
            mutation_calls=[call],
        )

        self.assertEqual(summary.status, 'COMPLETE')
        self.assertEqual(
            self.drm_screening_position_statuses(summary)[103],
            DRM_SCREENING_RESOLVED_MUTANT,
        )

    def test_drm_screening_marks_drm_position_not_covered(self):
        codons = self.drm_screening_codons('RT')
        codons[102]['alt_aa'] = 'X'
        codons[102]['codon_status'] = 'NOT_COVERED'

        summary = summarize_drm_screening(
            sequence_id='not_covered',
            gene='RT',
            codon_calls=codons,
            mutation_calls=[],
        )

        self.assertEqual(summary.status, 'PARTIAL')
        self.assertEqual(
            self.drm_screening_position_statuses(summary)[103],
            DRM_SCREENING_NOT_COVERED,
        )

    def test_drm_screening_marks_partial_codon_unresolved(self):
        codons = self.drm_screening_codons('RT')
        codons[102]['alt_aa'] = 'X'
        codons[102]['codon_status'] = 'PARTIAL_CODON'
        call = build_mutation_call(
            sequence_id='partial_codon', gene='RT', position=103, ref_aa='K', alt_aa='X',
            mutation_type='partial_codon', qc_status='WARN',
            codon_status='PARTIAL_CODON',
            reference_system='HIVDB_CONSENSUS_B',
        )

        summary = summarize_drm_screening(
            sequence_id='partial_codon',
            gene='RT',
            codon_calls=codons,
            mutation_calls=[call],
        )

        self.assertEqual(summary.status, 'PARTIAL')
        self.assertEqual(
            self.drm_screening_position_statuses(summary)[103],
            DRM_SCREENING_UNRESOLVED,
        )

    def test_drm_screening_resolved_mixture_is_resolved_mutant(self):
        codons = self.drm_screening_codons('RT')
        codons[102]['alt_aa'] = 'KN'
        codons[102]['possible_alt_aas'] = 'KN'
        call = build_mutation_call(
            sequence_id='mixture', gene='RT', position=103, ref_aa='K', alt_aa='KN',
            mutation_type='substitution', qc_status='PASS',
            reference_system='HIVDB_CONSENSUS_B',
        )

        summary = summarize_drm_screening(
            sequence_id='mixture',
            gene='RT',
            codon_calls=codons,
            mutation_calls=[call],
        )

        self.assertEqual(summary.status, 'COMPLETE')
        self.assertEqual(
            self.drm_screening_position_statuses(summary)[103],
            DRM_SCREENING_RESOLVED_MUTANT,
        )

    def test_drm_screening_combines_resolved_and_unresolved_positions(self):
        codons = self.drm_screening_codons('RT')
        unresolved_position = next(position for position in drm_catalog_positions('RT') if position != 103)
        codons[102]['alt_aa'] = 'N'
        codons[unresolved_position - 1]['alt_aa'] = 'X'
        codons[unresolved_position - 1]['codon_status'] = 'PARTIAL_CODON'
        call = build_mutation_call(
            sequence_id='mixed', gene='RT', position=103, ref_aa='K', alt_aa='N',
            mutation_type='substitution', qc_status='WARN',
            reference_system='HIVDB_CONSENSUS_B',
        )

        summary = summarize_drm_screening(
            sequence_id='mixed',
            gene='RT',
            codon_calls=codons,
            mutation_calls=[call],
        )
        statuses = self.drm_screening_position_statuses(summary)

        self.assertEqual(summary.status, 'PARTIAL')
        self.assertEqual(statuses[103], DRM_SCREENING_RESOLVED_MUTANT)
        self.assertEqual(statuses[unresolved_position], DRM_SCREENING_UNRESOLVED)

    def test_mixture_drm_tri_state_calls(self):
        cases = [
            ('WT+DRM', 'KN', True, 'COMPLETE'),
            ('non-DRM+DRM', 'QN', True, 'COMPLETE'),
            ('DRM+DRM', 'NS', True, 'COMPLETE'),
            ('DRM+X', 'NX', True, 'PARTIAL'),
            ('non-DRM+X', 'QX', None, 'PARTIAL'),
        ]
        for label, alt_aa, expected_is_drm, expected_evaluation in cases:
            with self.subTest(label=label):
                call = build_mutation_call(
                    sequence_id=label, gene='RT', position=103, ref_aa='K', alt_aa=alt_aa,
                    mutation_type='substitution', qc_status='PASS',
                    reference_system='HIVDB_CONSENSUS_B',
                )

                self.assertIs(call.is_drm, expected_is_drm)
                self.assertEqual(call.drm_evaluation_status, expected_evaluation)
                self.assertTrue(call.mixture)

    def test_ambiguous_resolved_iupac_mixture_keeps_drm_evidence(self):
        codons = self.drm_screening_codons('RT')
        codons[102]['alt_aa'] = 'KN'
        codons[102]['query_codon'] = 'AAN'
        codons[102]['possible_alt_aas'] = 'KN'
        codons[102]['codon_status'] = 'AMBIGUOUS'
        call = build_mutation_call(
            sequence_id='ambiguous', gene='RT', position=103, ref_aa='K', alt_aa='KN',
            mutation_type='substitution', qc_status='WARN',
            codon_status='AMBIGUOUS', query_codon='AAN', possible_alt_aas='KN',
            reference_system='HIVDB_CONSENSUS_B',
        )

        summary = summarize_drm_screening(
            sequence_id='ambiguous',
            gene='RT',
            codon_calls=codons,
            mutation_calls=[call],
        )

        self.assertIs(call.is_drm, True)
        self.assertEqual(call.drm_evaluation_status, 'COMPLETE')
        self.assertTrue(call.mixture)
        self.assertEqual(summary.status, 'COMPLETE')
        self.assertEqual(
            self.drm_screening_position_statuses(summary)[103],
            DRM_SCREENING_RESOLVED_MUTANT,
        )

    def test_x_call_keeps_possible_drm_evidence_separate(self):
        call = build_mutation_call(
            sequence_id='possible', gene='RT', position=215, ref_aa='T', alt_aa='X',
            mutation_type='substitution', qc_status='WARN',
            codon_status='AMBIGUOUS', query_codon='WYM', possible_alt_aas='FILST',
            reference_system='HIVDB_CONSENSUS_B',
        )

        self.assertEqual(call.alt_aa, 'X')
        self.assertIsNone(call.is_drm)
        self.assertEqual(call.drm_status, 'UNRESOLVED')
        self.assertTrue(call.possible_is_drm)
        self.assertEqual(call.possible_drm_class, 'NRTI')
        self.assertEqual(call.possible_drug_class, 'NRTI')
        self.assertTrue(any(component.aa == 'F' and component.is_drm for component in call.possible_drm_components))

    def test_partial_codon_possibilities_do_not_become_possible_drm(self):
        call = build_mutation_call(
            sequence_id='partial', gene='RT', position=215, ref_aa='T', alt_aa='X',
            mutation_type='partial_codon', qc_status='WARN',
            codon_status='PARTIAL_CODON', query_codon='W--', possible_alt_aas='FILST',
            reference_system='HIVDB_CONSENSUS_B',
        )

        self.assertIsNone(call.is_drm)
        self.assertFalse(call.possible_is_drm)
        self.assertEqual(call.possible_drm_components, ())

    def test_deletion_and_insertion_drm_events_are_annotated(self):
        deletion = build_mutation_call(
            sequence_id='deletion', gene='RT', position=69, ref_aa='T', alt_aa='-',
            mutation_type='deletion', qc_status='PASS',
            codon_status='DELETION', deleted_nt_count=3,
            reference_system='HIVDB_CONSENSUS_B',
        )
        insertion = build_mutation_call(
            sequence_id='insertion', gene='RT', position=69, ref_aa='-', alt_aa='A',
            mutation_type='insertion', qc_status='PASS',
            codon_status='INSERTION', inserted_nts='GCT',
            reference_system='HIVDB_CONSENSUS_B',
        )

        self.assertTrue(deletion.deletion)
        self.assertTrue(deletion.is_drm)
        self.assertEqual(deletion.mutation, 'RT:T69del')
        self.assertTrue(insertion.insertion)
        self.assertTrue(insertion.is_drm)
        self.assertEqual(insertion.mutation, 'RT:69insA')

    def test_partial_codon_deleted_nt_count_is_qc_not_amino_acid_deletion(self):
        call = build_mutation_call(
            sequence_id='partial', gene='PR', position=36, ref_aa='M', alt_aa='X',
            mutation_type='partial_codon', qc_status='WARN',
            codon_status='PARTIAL_CODON', deleted_nt_count=1,
            reference_system='HIVDB_CONSENSUS_B',
        )

        self.assertFalse(call.deletion)
        self.assertEqual(call.deleted_nt_count, 1)
        self.assertIn('PARTIAL_CODON', call.codon_status)
        self.assertEqual(call.drm_evaluation_status, 'UNRESOLVED')

    def test_drm_screening_summarizes_by_sequence_from_position_screening(self):
        complete = summarize_drm_screening(
            sequence_id='complete',
            gene='RT',
            codon_calls=self.drm_screening_codons('RT'),
            mutation_calls=[],
        )
        partial_codons = self.drm_screening_codons('RT')
        partial_codons[102]['alt_aa'] = 'X'
        partial_codons[102]['codon_status'] = 'PARTIAL_CODON'
        partial = summarize_drm_screening(
            sequence_id='partial',
            gene='RT',
            codon_calls=partial_codons,
            mutation_calls=[],
        )
        unresolved_codons = self.drm_screening_codons('RT')
        for position in drm_catalog_positions('RT'):
            unresolved_codons[position - 1]['alt_aa'] = 'X'
            unresolved_codons[position - 1]['codon_status'] = 'NOT_COVERED'
        unresolved = summarize_drm_screening(
            sequence_id='unresolved',
            gene='RT',
            codon_calls=unresolved_codons,
            mutation_calls=[],
        )

        self.assertEqual(
            drm_screening_status_by_sequence([complete, partial, unresolved]),
            {'complete': 'COMPLETE', 'partial': 'PARTIAL', 'unresolved': 'UNRESOLVED'},
        )

    def test_writes_drm_screening_position_and_summary_tsvs(self):
        codons = self.drm_screening_codons('RT')
        codons[102]['alt_aa'] = 'N'
        call = build_mutation_call(
            sequence_id='SEQ001', gene='RT', position=103, ref_aa='K', alt_aa='N',
            mutation_type='substitution', qc_status='PASS',
            reference_system='HIVDB_CONSENSUS_B',
        )
        summary = summarize_drm_screening(
            sequence_id='SEQ001',
            gene='RT',
            codon_calls=codons,
            mutation_calls=[call],
        )

        with TemporaryDirectory() as directory:
            positions_path = Path(directory) / 'drm_screening.tsv'
            summary_path = Path(directory) / 'drm_screening_summary.tsv'
            write_drm_screening_tsv([summary], positions_path, summary_path)
            with positions_path.open() as stream:
                positions = list(csv.DictReader(stream, delimiter='\t'))
            with summary_path.open() as stream:
                summaries = list(csv.DictReader(stream, delimiter='\t'))

        self.assertTrue(positions)
        by_position = {int(row['position']): row for row in positions}
        self.assertEqual(by_position[103]['sequence_id'], 'SEQ001')
        self.assertEqual(by_position[103]['gene'], 'RT')
        self.assertEqual(by_position[103]['status'], DRM_SCREENING_RESOLVED_MUTANT)
        self.assertEqual(summaries, [{
            'sequence_id': 'SEQ001',
            'gene': 'RT',
            'drm_evaluation_status': 'COMPLETE',
            'resolved_wt': str(len(summary.positions) - 1),
            'resolved_mutant': '1',
            'unresolved': '0',
            'not_covered': '0',
            'positions': str(len(summary.positions)),
        }])

    def test_append_drm_summary_reports_sequence_evaluation_status(self):
        import pandas as pd
        from pyhiv import append_drm_summary

        partial_codons = self.drm_screening_codons('RT')
        partial_codons[102]['alt_aa'] = 'X'
        partial_codons[102]['codon_status'] = 'NOT_COVERED'
        partial_summary = summarize_drm_screening(
            sequence_id='partial',
            gene='RT',
            codon_calls=partial_codons,
            mutation_calls=[],
        )
        complete = build_mutation_call(
            sequence_id='complete', gene='RT', position=103, ref_aa='K', alt_aa='N',
            mutation_type='substitution', qc_status='PASS',
            reference_system='HIVDB_CONSENSUS_B',
        )
        accessory = build_mutation_call(
            sequence_id='complete', gene='IN', position=128, ref_aa='A', alt_aa='T',
            mutation_type='substitution', qc_status='PASS',
            reference_system='HIVDB_CONSENSUS_B',
        )
        complete_codons = self.drm_screening_codons('RT')
        complete_codons[102]['alt_aa'] = 'N'
        complete_summary = summarize_drm_screening(
            sequence_id='complete',
            gene='RT',
            codon_calls=complete_codons,
            mutation_calls=[complete],
        )
        final_table = pd.DataFrame({'Sequence': ['partial', 'complete', 'no_calls']})

        result = append_drm_summary(
            final_table, [complete, accessory], [partial_summary, complete_summary]
        )

        by_sequence = result.set_index('Sequence')
        self.assertEqual(by_sequence.loc['partial', 'DRM Evaluation Status'], 'PARTIAL')
        self.assertEqual(by_sequence.loc['complete', 'DRM Evaluation Status'], 'COMPLETE')
        self.assertEqual(by_sequence.loc['no_calls', 'DRM Evaluation Status'], 'COMPLETE')
        self.assertEqual(by_sequence.loc['complete', 'DRM Count'], 1)
        self.assertEqual(by_sequence.loc['complete', 'Accessory Mutation Count'], 1)
        self.assertEqual(by_sequence.loc['complete', 'Accessory Mutations'], 'IN:A128T')
        self.assertEqual(by_sequence.loc['no_calls', 'Accessory Mutations'], '-')

    def test_position_qc_rows_cover_reference_positions_for_identity(self):
        calls, screening, rows = call_mutations_with_qc(
            reference_for_gene('PR'),
            sequence_id='SEQ001',
            gene='PR',
            sequence_type='aa',
        )

        self.assertEqual(calls, [])
        self.assertEqual(screening.status, 'COMPLETE')
        self.assertEqual(len(rows), len(reference_for_gene('PR')))
        self.assertEqual(rows[0].sequence_id, 'SEQ001')
        self.assertEqual(rows[0].gene, 'PR')
        self.assertEqual(rows[0].HXB2_position, 1)
        self.assertEqual(rows[0].Consensus_B_ref_aa, 'P')
        self.assertEqual(rows[0].observed_aa, 'P')
        self.assertEqual(rows[0].mutation, '-')
        self.assertEqual(rows[0].codon_status, 'OBSERVED')
        self.assertEqual(rows[0].coverage_status, 'COVERED')
        self.assertEqual(rows[0].ambiguity, 'NONE')
        self.assertFalse(rows[0].insertion)
        self.assertFalse(rows[0].deletion)
        self.assertFalse(rows[0].stop)
        self.assertFalse(rows[0].is_drm)
        self.assertEqual(rows[0].qc_status, 'PASS')

    def test_position_qc_rows_preserve_uncertainty_states(self):
        codons = self.drm_screening_codons('PR')
        codons[0] = {**codons[0], 'alt_aa': 'X', 'query_codon': '', 'codon_status': 'NOT_COVERED'}
        codons[1] = {**codons[1], 'alt_aa': 'X', 'query_codon': 'A--', 'codon_status': 'PARTIAL_CODON'}
        codons[2] = {**codons[2], 'alt_aa': 'X', 'query_codon': 'NNN', 'codon_status': 'AMBIGUOUS'}
        codons[3] = {**codons[3], 'alt_aa': '-', 'query_codon': '---', 'codon_status': 'DELETION'}
        insertion = build_mutation_call(
            sequence_id='SEQ001', gene='PR', position=5, ref_aa='-', alt_aa='A',
            mutation_type='insertion', qc_status='WARN',
            reference_system='HIVDB_CONSENSUS_B', codon_status='INSERTION',
            query_codon='GCT', inserted_nts='GCT',
        )

        rows = mutation_position_qc_rows(
            sequence_id='SEQ001',
            gene='PR',
            codon_calls=codons,
            mutation_calls=[insertion],
            reference_protein=reference_for_gene('PR'),
            qc_status='WARN',
        )

        self.assertEqual(rows[0].coverage_status, 'NOT_COVERED')
        self.assertEqual(rows[0].ambiguity, 'NOT_COVERED')
        self.assertEqual(rows[0].mutation, 'NOT_COVERED')
        self.assertEqual(rows[1].coverage_status, 'PARTIAL')
        self.assertEqual(rows[1].ambiguity, 'PARTIAL_CODON')
        self.assertEqual(rows[2].coverage_status, 'COVERED')
        self.assertEqual(rows[2].ambiguity, 'AMBIGUOUS')
        self.assertIsNone(rows[2].is_drm)
        self.assertTrue(rows[3].deletion)
        self.assertEqual(rows[3].codon_status, 'DELETION')
        self.assertTrue(rows[4].insertion)
        self.assertIn('INSERTION', rows[4].codon_status)
        self.assertIn('ins', rows[4].mutation)

    def test_mutation_input_qc_passes_iupac_nt_symbols(self):
        row = validate_mutation_input_sequence(
            "ACGTRYSWKMBDHVN-",
            sequence_id="SEQ001",
            gene="PR",
            sequence_type="nt",
        )

        self.assertEqual(row.sequence_id, "SEQ001")
        self.assertEqual(row.gene, "PR")
        self.assertEqual(row.sequence_type, "nt")
        self.assertEqual(row.input_qc_status, "PASS")
        self.assertEqual(row.invalid_symbols, "")
        self.assertEqual(row.invalid_symbol_count, 0)
        self.assertEqual(row.invalid_positions, "")
        self.assertEqual(row.invalid_runs, "")
        self.assertEqual(row.invalid_run_count, 0)
        self.assertEqual(row.invalid_run_length_class, "none")

    def test_mutation_input_qc_warns_on_non_iupac_nt_symbols(self):
        row = validate_mutation_input_sequence(
            "AC~T~~A~~~",
            sequence_id="SEQ001",
            gene="RT",
            sequence_type="nt",
        )

        self.assertEqual(row.input_qc_status, "WARN")
        self.assertEqual(row.invalid_symbols, "~")
        self.assertEqual(row.invalid_symbol_count, 6)
        self.assertEqual(row.invalid_positions, "3,5,6,8,9,10")
        self.assertEqual(row.invalid_runs, "3-3:~;5-6:~~;8-10:~~~")
        self.assertEqual(row.invalid_run_count, 3)
        self.assertEqual(row.invalid_run_length_class, "mixed")

    def test_mutation_input_qc_auto_flags_nt_like_sequence_with_invalid_symbol(self):
        row = validate_mutation_input_sequence(
            "ACGT~ACGT",
            sequence_id="SEQ001",
            gene="IN",
            sequence_type="auto",
        )

        self.assertEqual(row.sequence_type, "nt")
        self.assertEqual(row.input_qc_status, "WARN")
        self.assertEqual(row.invalid_positions, "5")

    def test_writes_mutation_input_qc_tsv(self):
        with TemporaryDirectory() as directory:
            fasta = Path(directory) / "rt.fasta"
            fasta.write_text(">SEQ001\nAC~T~~A~~~\n")
            rows = mutation_input_qc_for_fasta_files(
                [fasta],
                gene="RT",
                sequence_type="nt",
            )
            path = Path(directory) / "mutation_input_qc.tsv"
            write_mutation_input_qc_tsv(rows, path)
            with path.open() as stream:
                output_rows = list(csv.DictReader(stream, delimiter="	"))

        self.assertEqual(output_rows[0]["sequence_id"], "SEQ001")
        self.assertEqual(output_rows[0]["gene"], "RT")
        self.assertEqual(output_rows[0]["sequence_type"], "nt")
        self.assertEqual(output_rows[0]["input_qc_status"], "WARN")
        self.assertEqual(output_rows[0]["invalid_symbols"], "~")
        self.assertEqual(output_rows[0]["invalid_symbol_count"], "6")
        self.assertEqual(output_rows[0]["invalid_positions"], "3,5,6,8,9,10")
        self.assertEqual(output_rows[0]["invalid_runs"], "3-3:~;5-6:~~;8-10:~~~")
        self.assertEqual(output_rows[0]["invalid_run_count"], "3")
        self.assertEqual(output_rows[0]["invalid_run_length_class"], "mixed")
        self.assertEqual(list(output_rows[0]), MUTATION_INPUT_QC_TSV_COLUMNS)

    def test_writes_mutation_position_qc_tsv(self):
        _calls, _screening, rows = call_mutations_with_qc(
            reference_for_gene('PR'),
            sequence_id='SEQ001',
            gene='PR',
            sequence_type='aa',
        )
        with TemporaryDirectory() as directory:
            path = Path(directory) / 'mutation_position_qc.tsv'
            write_mutation_position_qc_tsv(rows, path)
            with path.open() as stream:
                output_rows = list(csv.DictReader(stream, delimiter='	'))

        self.assertEqual(output_rows[0]['sequence_id'], 'SEQ001')
        self.assertEqual(output_rows[0]['gene'], 'PR')
        self.assertEqual(output_rows[0]['HXB2_position'], '1')
        self.assertEqual(output_rows[0]['Consensus_B_ref_aa'], 'P')
        self.assertEqual(output_rows[0]['mutation'], '-')
        self.assertEqual(output_rows[0]['coverage_status'], 'COVERED')
        self.assertEqual(output_rows[0]['ambiguity'], 'NONE')
        self.assertEqual(output_rows[0]['is_drm'], 'false')
        self.assertEqual(list(output_rows[0]), MUTATION_POSITION_QC_TSV_COLUMNS)

    def test_writes_gene_mutation_matrix_with_reference_row(self):
        calls = [
            build_mutation_call(
                sequence_id='SEQ001', gene='RT', position=184, ref_aa='M', alt_aa='V',
                mutation_type='substitution', qc_status='PASS',
                reference_system='HIVDB_CONSENSUS_B',
            ),
            build_mutation_call(
                sequence_id='SEQ001', gene='RT', position=184, ref_aa='M', alt_aa='I',
                mutation_type='substitution', qc_status='PASS',
                reference_system='HIVDB_CONSENSUS_B',
            ),
        ]
        with TemporaryDirectory() as directory:
            output_dir = Path(directory)
            paths = write_mutation_matrices_tsv(
                calls,
                output_dir,
                sequence_gene_pairs=[('SEQ001', 'RT'), ('SEQ002', 'RT')],
            )
            matrix_path = output_dir / 'mutation_matrix_RT.tsv'
            with matrix_path.open() as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))

        self.assertEqual(paths, [matrix_path])
        self.assertEqual(rows[0]['Sequence'], 'REFERENCE')
        self.assertEqual(rows[0]['RT:M184'], 'M')
        self.assertEqual(rows[1]['Sequence'], 'SEQ001')
        self.assertEqual(rows[1]['RT:M184'], 'M184I, M184V')
        self.assertEqual(rows[2]['Sequence'], 'SEQ002')
        self.assertEqual(rows[2]['RT:M184'], '-')

    def test_annotation_uncertainty_and_provenance_survive_tsv_export(self):
        from pyhiv.mutations.drm import UNRESOLVED_CODON_STATUSES, LOCAL_DRM_VERSION

        calls = [build_mutation_call(
            sequence_id=status, gene='PR', position=7, ref_aa='Q', alt_aa='V',
            mutation_type='substitution', qc_status='WARN',
            codon_status=status, possible_alt_aas='V',
            reference_system='HIVDB_CONSENSUS_B',
        ) for status in sorted(UNRESOLVED_CODON_STATUSES)]
        for call in calls:
            self.assertIsNone(call.is_drm)
            self.assertEqual(call.drm_status, 'UNRESOLVED')
        with TemporaryDirectory() as directory:
            path = Path(directory) / 'calls.tsv'
            write_mutations_tsv(calls, path)
            with path.open() as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
        for row in rows:
            self.assertEqual(row['is_drm'], 'UNRESOLVED')
            self.assertEqual(row['drm_status'], 'UNRESOLVED')
            self.assertEqual(row['drm_evaluation_status'], 'UNRESOLVED')
            self.assertEqual(row['drm_version'], LOCAL_DRM_VERSION)
            self.assertEqual(len(row['drm_catalog_sha256']), 64)
            component = json.loads(row['drm_components'])[0]
            self.assertIsNone(component['is_drm'])
            self.assertEqual(component['drm_status'], 'UNRESOLVED')
            self.assertEqual(row['reference_system'], 'HIVDB_CONSENSUS_B')
            self.assertEqual(row['coordinate_system'], 'HXB2')

    def test_drm_catalogue_cannot_change_alignment_or_mutation_calls(self):
        from dataclasses import asdict
        from pyhiv.mutations.drm import DRMAnnotation

        # Same rule across several genes and positions, without knowing DRMs.
        for gene in ('PR', 'RT', 'IN'):
            sequence = list(hxb2_reference_nt_for_gene(gene))
            for index in range(12, len(sequence) - 12, 37):
                sequence[index] = next(base for base in 'ACGT' if base != sequence[index])
            results = []
            for annotation in (DRMAnnotation(False, 'test', '', ''),
                               DRMAnnotation(True, 'test', 'test', 'test')):
                with patch('pyhiv.mutations.caller.annotate_drm', return_value=annotation):
                    calls = call_mutations(''.join(sequence), 'test', gene, sequence_type='nt')
                results.append([
                    {key: value for key, value in asdict(call).items()
                     if key not in {'is_drm', 'drm_source', 'drm_class', 'drug_class'}}
                    for call in calls
                ])
            with self.subTest(gene=gene):
                self.assertTrue(results[0])
                self.assertEqual(results[0], results[1])

    def test_calls_rt_substitutions_and_deletion(self):
        sequence = list(reference_for_gene("RT"))
        sequence[102] = "N"
        sequence[183] = "V"
        del sequence[68]

        calls = call_mutations("".join(sequence), sequence_id="SEQ001", gene="RT", sequence_type="aa")
        by_mutation = {call.mutation: call for call in calls}

        self.assertIn("RT:K103N", by_mutation)
        self.assertIn("RT:M184V", by_mutation)
        self.assertIn("RT:T69del", by_mutation)
        self.assertTrue(by_mutation["RT:K103N"].is_drm)
        self.assertEqual(by_mutation["RT:K103N"].drm_class, "NNRTI")
        self.assertTrue(by_mutation["RT:M184V"].is_drm)
        self.assertEqual(by_mutation["RT:M184V"].drm_class, "NRTI")
        self.assertTrue(by_mutation["RT:T69del"].deletion)

    def test_reports_insertion_from_alignment(self):
        calls = mutations_from_alignment(
            sequence_id="SEQ001",
            gene="PR",
            query_aligned="PQXIT",
            ref_aligned="PQ-IT",
            qc_status="PASS",
        )

        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0].mutation, "PR:2insX")
        self.assertTrue(calls[0].insertion)
        self.assertEqual(calls[0].position, 2)

    def test_marks_stop(self):
        sequence = list(reference_for_gene("RT"))
        sequence[183] = "*"

        calls = call_mutations("".join(sequence), sequence_id="SEQ001", gene="RT", sequence_type="aa")

        self.assertEqual(calls[0].mutation, "RT:M184*")
        self.assertTrue(calls[0].stop)

    def test_translates_gapped_nucleotide_alignment_without_shifting_frame(self):
        protein, qc_status = protein_sequence_for_calling(
            "ATG-GTT",
            reference_protein="MV",
            sequence_type="nt",
        )

        self.assertEqual(protein, "MV")
        self.assertEqual(qc_status, "PASS")

    def test_calls_reference_aligned_codon_deletion_without_frameshift(self):
        sequence = list(hxb2_reference_nt_for_gene("PR"))
        sequence[3:6] = "---"

        calls = mutations_from_reference_aligned_dna(
            sequence_id="SEQ001",
            gene="PR",
            sequence="".join(sequence),
            reference_protein=reference_for_gene("PR"),
        )

        by_mutation = {call.mutation: call for call in calls}
        self.assertIn("PR:Q2del", by_mutation)
        self.assertEqual(by_mutation["PR:Q2del"].qc_status, "PASS")
        self.assertEqual(by_mutation["PR:Q2del"].hxb2_ref_codon, "CAG")
        self.assertEqual(by_mutation["PR:Q2del"].query_codon, "---")
        self.assertEqual(by_mutation["PR:Q2del"].codon_status, "DELETION")

    def test_one_nt_deletion_marks_frameshift_without_position_specific_rules(self):
        for position in (36, 63):
            sequence = list(hxb2_reference_nt_for_gene("PR"))
            sequence[(position - 1) * 3 + 1] = "-"

            calls = mutations_from_reference_aligned_dna(
                sequence_id=f"SEQ{position}",
                gene="PR",
                sequence="".join(sequence),
                reference_protein=reference_for_gene("PR"),
            )
            call = next(call for call in calls if call.position == position)

            self.assertEqual(call.alt_aa, "X")
            self.assertEqual(call.mutation_type, "frameshift")
            self.assertIn("PARTIAL_CODON", call.codon_status)
            self.assertIn("FRAMESHIFT", call.codon_status)

    def test_frameshift_is_marked_until_compensating_indel_restores_phase(self):
        hxb2_ref = hxb2_reference_nt_for_gene("PR")
        query = list(hxb2_ref)
        query[(10 - 1) * 3 + 1] = "-"
        insertion_column = 11 * 3
        hxb2_ref_aligned = hxb2_ref[:insertion_column] + "-" + hxb2_ref[insertion_column:]
        query_aligned = "".join(query[:insertion_column]) + "A" + "".join(query[insertion_column:])

        calls = mutations_from_hxb2_aligned_dna(
            sequence_id="SEQ001",
            gene="PR",
            query_aligned=query_aligned,
            hxb2_ref_aligned=hxb2_ref_aligned,
            reference_protein=reference_for_gene("PR"),
        )
        by_position = {call.position: call for call in calls}

        self.assertIn("FRAMESHIFT", by_position[10].codon_status)
        self.assertIn("FRAMESHIFT_RESTORED", by_position[11].codon_status)
        self.assertNotIn(12, by_position)

    def test_codon_deletion_is_not_treated_as_frameshift_at_multiple_positions(self):
        for position in (2, 63):
            sequence = list(hxb2_reference_nt_for_gene("PR"))
            start = (position - 1) * 3
            sequence[start:start + 3] = "---"

            calls = mutations_from_reference_aligned_dna(
                sequence_id=f"SEQ{position}",
                gene="PR",
                sequence="".join(sequence),
                reference_protein=reference_for_gene("PR"),
            )
            call = next(call for call in calls if call.position == position)

            self.assertTrue(call.deletion)
            self.assertEqual(call.codon_status, "DELETION")

    def test_preserves_iupac_amino_acid_mixture_in_coding_region(self):
        sequence = list(hxb2_reference_nt_for_gene("PR"))
        sequence[42:45] = "ATN"

        calls = mutations_from_reference_aligned_dna(
            sequence_id="SEQ001",
            gene="PR",
            sequence="".join(sequence),
            reference_protein=reference_for_gene("PR"),
        )
        by_mutation = {call.mutation: call for call in calls}

        self.assertIn("PR:I15IM", by_mutation)
        self.assertEqual(by_mutation["PR:I15IM"].possible_alt_aas, "IM")

    def test_preserves_full_iupac_possible_amino_acids_when_final_call_is_x(self):
        for position in (15, 72):
            sequence = list(hxb2_reference_nt_for_gene("PR"))
            start = (position - 1) * 3
            sequence[start:start + 3] = "NNK"

            calls = mutations_from_reference_aligned_dna(
                sequence_id=f"SEQ{position}",
                gene="PR",
                sequence="".join(sequence),
                reference_protein=reference_for_gene("PR"),
            )
            call = next(call for call in calls if call.position == position)

            self.assertEqual(call.alt_aa, "X")
            self.assertGreater(len(call.possible_alt_aas), 4)
            self.assertNotEqual(call.possible_alt_aas, "X")

    def test_reports_codon_insertions_from_existing_hxb2_alignment(self):
        hxb2_ref = hxb2_reference_nt_for_gene("PR")
        for position in (2, 63):
            insertion_column = position * 3
            hxb2_ref_aligned = hxb2_ref[:insertion_column] + "---" + hxb2_ref[insertion_column:]
            query_aligned = hxb2_ref[:insertion_column] + "GCT" + hxb2_ref[insertion_column:]

            calls = mutations_from_hxb2_aligned_dna(
                sequence_id=f"SEQ{position}",
                gene="PR",
                query_aligned=query_aligned,
                hxb2_ref_aligned=hxb2_ref_aligned,
                reference_protein=reference_for_gene("PR"),
            )

            self.assertIn(f"PR:{position}insA", {call.mutation for call in calls})

    def test_ignores_terminal_unsequenced_codons(self):
        sequence = list(hxb2_reference_nt_for_gene("PR"))
        sequence[:3] = "NNN"
        sequence[-3:] = "NNN"

        calls = mutations_from_reference_aligned_dna(
            sequence_id="SEQ001",
            gene="PR",
            sequence="".join(sequence),
            reference_protein=reference_for_gene("PR"),
        )

        mutations = {call.mutation for call in calls}
        self.assertNotIn("PR:P1X", mutations)
        self.assertNotIn("PR:F99X", mutations)

    def test_uses_hxb2_coordinates_but_consensus_b_reference_amino_acids(self):
        calls = call_mutations(
            hxb2_reference_nt_for_gene("PR"),
            sequence_id="HXB2",
            gene="PR",
            sequence_type="nt",
        )

        self.assertIn("PR:I3V", {call.mutation for call in calls})

    def test_raw_polymorphic_pr_does_not_invent_frameshifts(self):
        sequence = list(hxb2_reference_nt_for_gene("PR"))
        # Introduce many substitutions without changing sequence length. The
        # raw aligner must keep these as mismatches rather than inventing 1-nt
        # indels that would create long frameshifts.
        for index in range(5, len(sequence) - 5, 11):
            original = sequence[index]
            sequence[index] = next(base for base in "ACGT" if base != original)

        calls = call_mutations(
            "".join(sequence),
            sequence_id="POLYMORPHIC",
            gene="PR",
            sequence_type="nt",
        )

        self.assertFalse(any(call.mutation_type == "frameshift" for call in calls))

    def test_raw_one_nt_deletion_is_local_uncertainty_not_propagated(self):
        sequence = hxb2_reference_nt_for_gene("PR")
        deleted = sequence[:90] + sequence[91:]

        calls = call_mutations(
            deleted,
            sequence_id="RAW_FS",
            gene="PR",
            sequence_type="nt",
        )

        self.assertTrue(any("PARTIAL_CODON" in call.codon_status for call in calls))
        self.assertFalse(any(call.mutation_type == "frameshift" for call in calls))
        self.assertLess(sum(call.alt_aa == "X" for call in calls), 10)

    def test_raw_three_nt_deletion_stays_in_frame(self):
        sequence = hxb2_reference_nt_for_gene("PR")
        deleted = sequence[:90] + sequence[93:]

        calls = call_mutations(
            deleted,
            sequence_id="RAW_INFRAME_DEL",
            gene="PR",
            sequence_type="nt",
        )

        self.assertTrue(any(call.deletion for call in calls))
        self.assertFalse(any(call.mutation_type == "frameshift" for call in calls))

    def test_codonaware_gap_score_keeps_terminal_gaps_free(self):
        self.assertEqual(internal_gap_score(start=0, length=10, sequence_length=100), 0)
        self.assertEqual(internal_gap_score(start=100, length=10, sequence_length=100), 0)

    def test_codonaware_gap_score_adds_frameshift_penalty_to_nontriplet_internal_gaps(self):
        in_frame = internal_gap_score(start=10, length=3, sequence_length=100)
        one_nt = internal_gap_score(start=10, length=1, sequence_length=100)
        two_nt = internal_gap_score(start=10, length=2, sequence_length=100)

        self.assertEqual(in_frame, -22)
        self.assertEqual(one_nt, -32)
        self.assertEqual(two_nt, -33)

    def test_semiglobal_codonaware_removes_terminal_query_overhangs(self):
        query_aligned, ref_aligned = align_nucleotides_semiglobal_codonaware(
            query="TTTAAACCCGGGTTT",
            reference="AAACCCGGG",
        )

        self.assertEqual(query_aligned.replace("-", ""), "AAACCCGGG")
        self.assertEqual(ref_aligned.replace("-", ""), "AAACCCGGG")
        self.assertEqual(query_aligned, ref_aligned)

    def test_semiglobal_codonaware_backend_is_experimental_and_not_default(self):
        reference = hxb2_reference_nt_for_gene("PR")
        query_aligned, ref_aligned = align_dna_to_hxb2_columns(
            reference,
            reference,
            alignment_backend="semiglobal_codonaware",
        )

        self.assertEqual(query_aligned, reference)
        self.assertEqual(ref_aligned, reference)

    def test_semiglobal_codonaware_aligns_identity_exactly(self):
        reference = hxb2_reference_nt_for_gene("PR")
        query_aligned, ref_aligned = align_nucleotides_semiglobal_codonaware(
            query=reference,
            reference=reference,
        )

        self.assertEqual(query_aligned, reference)
        self.assertEqual(ref_aligned, reference)

    def test_semiglobal_codonaware_handles_internal_one_two_and_three_nt_insertions(self):
        reference = hxb2_reference_nt_for_gene("PR")[:180]
        for inserted in ["A", "AG", "AGT"]:
            with self.subTest(inserted=inserted):
                query = reference[:90] + inserted + reference[90:]
                query_aligned, ref_aligned = align_nucleotides_semiglobal_codonaware(
                    query=query,
                    reference=reference,
                )

                self.assertEqual(query_aligned.replace("-", ""), query)
                self.assertEqual(ref_aligned.replace("-", ""), reference)
                self.assertIn(len(inserted), self.gap_lengths(query_aligned, ref_aligned, "reference"))

    def test_semiglobal_codonaware_handles_internal_one_two_and_three_nt_deletions(self):
        reference = hxb2_reference_nt_for_gene("PR")[:180]
        for deletion_length in [1, 2, 3]:
            with self.subTest(deletion_length=deletion_length):
                query = reference[:90] + reference[90 + deletion_length:]
                query_aligned, ref_aligned = align_nucleotides_semiglobal_codonaware(
                    query=query,
                    reference=reference,
                )

                self.assertEqual(query_aligned.replace("-", ""), query)
                self.assertEqual(ref_aligned.replace("-", ""), reference)
                self.assertIn(deletion_length, self.gap_lengths(query_aligned, ref_aligned, "query"))

    def test_semiglobal_codonaware_allows_large_internal_deletions(self):
        reference = hxb2_reference_nt_for_gene("RT")[:900]
        query = reference[:300] + reference[450:]

        query_aligned, ref_aligned = align_nucleotides_semiglobal_codonaware(
            query=query,
            reference=reference,
        )

        self.assertEqual(query_aligned.replace("-", ""), query)
        self.assertEqual(ref_aligned.replace("-", ""), reference)
        self.assertIn(150, self.gap_lengths(query_aligned, ref_aligned, "query"))

    def test_semiglobal_codonaware_allows_real_frameshift(self):
        reference = hxb2_reference_nt_for_gene("PR")
        query = reference[:90] + reference[91:]

        calls = call_mutations(
            query,
            sequence_id="CODONAWARE_REAL_FRAMESHIFT",
            gene="PR",
            sequence_type="nt",
            alignment_backend="semiglobal_codonaware",
        )

        self.assertTrue(any(call.deleted_nt_count == 1 for call in calls))
        self.assertTrue(any(call.codon_status == "PARTIAL_CODON" for call in calls))

    def test_semiglobal_codonaware_never_lowers_primary_score(self):
        reference = hxb2_reference_nt_for_gene("PR")
        queries = [
            reference,
            "TTT" + reference + "GGG",
            reference[:90] + "A" + reference[90:],
            reference[:90] + "AG" + reference[90:],
            reference[:90] + "AGT" + reference[90:],
            reference[:90] + reference[91:],
            reference[:90] + reference[92:],
            reference[:90] + reference[93:],
            reference[:60] + reference[150:],
        ]
        for query in queries:
            with self.subTest(query_length=len(query)):
                baseline_query_aligned, baseline_ref_aligned = align_dna_to_hxb2_columns(
                    query,
                    reference,
                    alignment_backend="semiglobal_mismatch2_open20",
                )
                codonaware_query_aligned, codonaware_ref_aligned = align_dna_to_hxb2_columns(
                    query,
                    reference,
                    alignment_backend="semiglobal_codonaware",
                )

                self.assertGreaterEqual(
                    self.semiglobal_primary_score(codonaware_query_aligned, codonaware_ref_aligned),
                    self.semiglobal_primary_score(baseline_query_aligned, baseline_ref_aligned),
                )

    def test_semiglobal_codonaware_local_is_experimental_and_preserves_primary_score(self):
        reference = hxb2_reference_nt_for_gene("PR")
        queries = [
            reference,
            "TTT" + reference + "GGG",
            reference[:90] + "A" + reference[90:],
            reference[:90] + "AG" + reference[90:],
            reference[:90] + reference[91:],
            reference[:60] + reference[150:],
        ]
        for query in queries:
            with self.subTest(query_length=len(query)):
                baseline_query_aligned, baseline_ref_aligned = align_dna_to_hxb2_columns(
                    query,
                    reference,
                    alignment_backend="semiglobal_mismatch2_open20",
                )
                local_query_aligned, local_ref_aligned = align_dna_to_hxb2_columns(
                    query,
                    reference,
                    alignment_backend="semiglobal_codonaware_local",
                )

                self.assertEqual(local_ref_aligned.replace("-", ""), reference)
                self.assertIn(local_query_aligned.replace("-", ""), query.strip().upper().replace("-", ""))
                self.assertGreaterEqual(
                    self.semiglobal_primary_score(local_query_aligned, local_ref_aligned),
                    self.semiglobal_primary_score(baseline_query_aligned, baseline_ref_aligned),
                )

    def test_semiglobal_segmented_is_experimental_and_handles_identity(self):
        reference = hxb2_reference_nt_for_gene("PR")
        query_aligned, ref_aligned = align_dna_to_hxb2_columns(
            reference,
            reference,
            alignment_backend="semiglobal_segmented",
        )

        self.assertEqual(query_aligned, reference)
        self.assertEqual(ref_aligned, reference)

    def test_semiglobal_segmented_trims_terminal_overhangs(self):
        reference = hxb2_reference_nt_for_gene("PR")
        query = "TTT" + reference + "GGG"
        query_aligned, ref_aligned = align_dna_to_hxb2_columns(
            query,
            reference,
            alignment_backend="semiglobal_segmented",
        )

        self.assertEqual(query_aligned.replace("-", ""), reference)
        self.assertEqual(ref_aligned.replace("-", ""), reference)

    def test_semiglobal_segmented_preserves_query_reference_and_allows_large_deletion(self):
        reference = hxb2_reference_nt_for_gene("RT")[:900]
        query = reference[:300] + reference[450:]
        query_aligned, ref_aligned = align_dna_to_hxb2_columns(
            query,
            reference,
            alignment_backend="semiglobal_segmented",
        )

        self.assertEqual(query_aligned.replace("-", ""), query)
        self.assertEqual(ref_aligned.replace("-", ""), reference)
        self.assertIn(150, self.gap_lengths(query_aligned, ref_aligned, "query"))

    def test_semiglobal_anchorblocks_is_experimental_and_handles_identity(self):
        reference = hxb2_reference_nt_for_gene("PR")
        query_aligned, ref_aligned = align_dna_to_hxb2_columns(
            reference,
            reference,
            alignment_backend="semiglobal_anchorblocks",
        )

        self.assertEqual(query_aligned, reference)
        self.assertEqual(ref_aligned, reference)

    def test_semiglobal_anchorblocks_marks_unsupported_internal_interval_with_n(self):
        reference = hxb2_reference_nt_for_gene("RT")[:900]
        query = reference[:300] + reference[450:]
        query_aligned, ref_aligned = align_dna_to_hxb2_columns(
            query,
            reference,
            alignment_backend="semiglobal_anchorblocks",
        )

        self.assertEqual(query_aligned.replace("-", "").replace("N", ""), query)
        self.assertEqual(ref_aligned.replace("-", ""), reference)
        self.assertIn("N" * 45, query_aligned)

    def test_semiglobal_mismatch2_open20_trimmed_handles_identity(self):
        reference = hxb2_reference_nt_for_gene("PR")
        query_aligned, ref_aligned = align_dna_to_hxb2_columns(
            reference,
            reference,
            alignment_backend="semiglobal_mismatch2_open20_trimmed",
        )

        self.assertEqual(query_aligned, reference)
        self.assertEqual(ref_aligned, reference)

    def test_semiglobal_mismatch2_open20_trimmed_removes_terminal_overhangs(self):
        reference = hxb2_reference_nt_for_gene("PR")
        query = "TTT" + reference + "GGG"
        query_aligned, ref_aligned = align_dna_to_hxb2_columns(
            query,
            reference,
            alignment_backend="semiglobal_mismatch2_open20_trimmed",
        )

        self.assertEqual(query_aligned, reference)
        self.assertEqual(ref_aligned, reference)

    def test_semiglobal_mismatch2_open20_trimmed_preserves_internal_deletion(self):
        reference = hxb2_reference_nt_for_gene("PR")
        query = reference[:60] + reference[90:]
        query_aligned, ref_aligned = align_dna_to_hxb2_columns(
            query,
            reference,
            alignment_backend="semiglobal_mismatch2_open20_trimmed",
        )

        self.assertEqual(query_aligned.replace("-", ""), query)
        self.assertEqual(ref_aligned.replace("-", ""), reference)
        self.assertIn("-" * 30, query_aligned)

    def test_semiglobal_mismatch2_open20_trimmed_does_not_report_terminal_overhang_insertions(self):
        reference = hxb2_reference_nt_for_gene("PR")
        calls = call_mutations(
            "TTT" + reference + "GGG",
            sequence_id="SEQ001",
            gene="PR",
            sequence_type="nt",
            alignment_backend="semiglobal_mismatch2_open20_trimmed",
        )

        self.assertFalse(any(call.insertion for call in calls))

    def test_nt_mutation_calling_default_uses_trimmed_semiglobal_backend(self):
        reference = hxb2_reference_nt_for_gene("PR")
        calls = call_mutations(
            "TTT" + reference + "GGG",
            sequence_id="SEQ001",
            gene="PR",
            sequence_type="nt",
        )

        self.assertFalse(any(call.insertion for call in calls))
        self.assertEqual(
            DEFAULT_NT_ALIGNMENT_BACKEND,
            "semiglobal_mismatch2_open20_trimmed",
        )

    def test_current_alignment_benchmark_result_is_regression_fixture(self):
        benchmark = (
            Path(__file__).resolve().parents[3]
            / "analysis_outputs"
            / "pol2024_final_trimmed_2026-09-25"
            / "comparison_mismatch2_vs_trimmed.tsv"
        )
        with benchmark.open() as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))

        trimmed = {
            row["gene"]: row
            for row in rows
            if row["variant"] == "semiglobal_mismatch2_open20_trimmed"
        }

        self.assertEqual(set(trimmed), {"PR", "RT", "IN"})
        self.assertEqual(trimmed["PR"]["insertion_calls"], "0")
        self.assertEqual(trimmed["RT"]["insertion_calls"], "4")
        self.assertEqual(trimmed["IN"]["insertion_calls"], "58")
        self.assertEqual(trimmed["PR"]["drm_only_pyhiv"], "0")
        self.assertEqual(trimmed["PR"]["drm_only_sierra"], "0")

    def test_builtin_alignment_backend_remains_available_explicitly(self):
        reference = hxb2_reference_nt_for_gene("PR")
        calls = call_mutations(
            "TTT" + reference + "GGG",
            sequence_id="SEQ001",
            gene="PR",
            sequence_type="nt",
            alignment_backend="builtin",
        )

        self.assertTrue(any(call.insertion for call in calls))

    def test_semiglobal_coverage_confidence_keeps_trimmed_baseline_alignment(self):
        reference = hxb2_reference_nt_for_gene("PR")
        query = reference[:60] + reference[180:]
        baseline_query, baseline_ref = align_dna_to_hxb2_columns(
            query,
            reference,
            alignment_backend="semiglobal_mismatch2_open20",
        )
        coverage_query, coverage_ref = align_dna_to_hxb2_columns(
            query,
            reference,
            alignment_backend="semiglobal_coverage_confidence",
        )
        ref_columns = [index for index, base in enumerate(baseline_ref) if base != "-"]
        baseline_query = baseline_query[ref_columns[0]:ref_columns[-1] + 1]
        baseline_ref = baseline_ref[ref_columns[0]:ref_columns[-1] + 1]

        self.assertEqual(coverage_query, baseline_query)
        self.assertEqual(coverage_ref, baseline_ref)

    def test_semiglobal_coverage_confidence_does_not_mask_identity(self):
        reference = hxb2_reference_nt_for_gene("PR")
        query_aligned, ref_aligned = align_dna_to_hxb2_columns(
            reference,
            reference,
            alignment_backend="semiglobal_coverage_confidence",
        )
        mask = coverage_confidence_not_covered_codons(
            query_aligned,
            ref_aligned,
            reference_for_gene("PR"),
        )

        self.assertEqual(mask, set())

    def test_semiglobal_coverage_confidence_marks_long_anchorless_region_not_covered(self):
        reference = hxb2_reference_nt_for_gene("PR")
        reference_protein = reference_for_gene("PR")
        query = reference[:60] + reference[180:]
        query_aligned, ref_aligned = align_dna_to_hxb2_columns(
            query,
            reference,
            alignment_backend="semiglobal_coverage_confidence",
        )
        mask = coverage_confidence_not_covered_codons(
            query_aligned,
            ref_aligned,
            reference_protein,
        )
        codons = codon_calls_from_raw_nucleotide_alignment(
            query_aligned,
            ref_aligned,
            reference_protein,
            force_not_covered_codons=mask,
        )

        self.assertTrue(mask)
        self.assertTrue(any(call["codon_status"] == "NOT_COVERED" for call in codons))
        self.assertFalse(any(
            call["codon_status"] == "DELETION"
            for index, call in enumerate(codons)
            if index in mask
        ))

    def test_infers_gene_from_pyhiv_split_path(self):
        self.assertEqual(infer_gene_from_path(Path("results/pol/p51-rt/sample.fasta")), "RT")
        self.assertEqual(infer_gene_from_path(Path("results/pol/p51 (rt)/sample.fasta")), "RT")
        self.assertEqual(infer_gene_from_path(Path("results/pol/protease/sample.fasta")), "PR")
        self.assertEqual(infer_gene_from_path(Path("results/pol/p31 (integrase)/sample.fasta")), "IN")
        self.assertEqual(infer_gene_from_path(Path("results/gag/p24-capsid/sample.fasta")), "CA")
        self.assertEqual(infer_gene_from_path(Path("results/gag/p24 (capsid)/sample.fasta")), "CA")

    def test_writes_mutations_tsv(self):
        sequence = list(reference_for_gene("RT"))
        sequence[183] = "V"
        calls = call_mutations("".join(sequence), sequence_id="SEQ001", gene="RT", sequence_type="aa")

        output = Path("tests") / "output_mutations.tsv"
        try:
            write_mutations_tsv(calls, output)
            with open(output, newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        finally:
            output.unlink(missing_ok=True)

        self.assertEqual(rows[0]["sequence_id"], "SEQ001")
        self.assertEqual(rows[0]["mutation"], "RT:M184V")
        self.assertEqual(rows[0]["is_drm"], "true")
        self.assertEqual(rows[0]["reference_system"], "HIVDB_CONSENSUS_B")
