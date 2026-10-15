"""Exhaustive annotation checks against raw official data, without alignments."""
import hashlib
import json
from itertools import combinations
from pathlib import Path
from shutil import copytree
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import patch

from pyhiv.mutations.drm import (
    LOCAL_DRM_SOURCE, LOCAL_DRM_VERSION, UNRESOLVED_CODON_STATUSES,
    annotate_drm, load_catalog,
)
from pyhiv.mutations.references import CONSENSUS_B_AMINO_ACIDS

DATA = Path(__file__).resolve().parents[3] / 'src/pyhiv/mutations/data'
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


class TestDRMAnnotation(TestCase):
    @classmethod
    def setUpClass(cls):
        # Expected results come straight from BOTH official JSON files, never
        # from the production catalogue loader (which could omit entries).
        cls.raw_drms = json.loads((DATA / 'drms_hiv1.json').read_text())
        cls.raw_types = json.loads((DATA / 'mutation-type-pairs_hiv1.json').read_text())
        cls.membership = {}
        cls.types = {}
        cls.types_by_drug = {}
        for drug, rows in cls.raw_drms.items():
            for row in rows:
                key = row['gene'], row['position'], row['aa']
                cls.membership.setdefault(key, set()).add(drug)
        for row in cls.raw_types:
            for aa in row['aas']:
                key = row['gene'], row['position'], aa
                cls.types.setdefault(key, set()).add(row['mutationType'])
                cls.types_by_drug.setdefault((*key, row['drugClass']), set()).add(row['mutationType'])

    def assert_official_annotation(self, gene, position, aa):
        key = gene, position, aa
        expected_drugs = self.membership.get(key, set())
        expected_types = self.types.get(key, {'Other'})
        expected_drm_classes = set()
        for drug in expected_drugs:
            expected_drm_classes.update(self.types_by_drug.get((*key, drug), {'Other'}))
        actual = annotate_drm(gene, position, aa)
        self.assertIs(actual.is_drm, bool(expected_drugs), key)
        self.assertEqual(actual.drm_status, 'DRM' if expected_drugs else 'NON_DRM', key)
        self.assertEqual(actual.drug_class, ';'.join(sorted(expected_drugs)), key)
        self.assertEqual(actual.drm_class, ';'.join(sorted(expected_drm_classes)), key)
        self.assertEqual(actual.mutation_types, ';'.join(sorted(expected_types)), key)
        self.assertIs(actual.is_accessory, 'Accessory' in expected_types, key)
        self.assertEqual(len(actual.components), 1)
        self.assertIs(actual.components[0].is_drm, actual.is_drm)
        self.assertEqual(actual.components[0].mutation_types, actual.mutation_types)
        return actual

    def test_every_official_drm_and_mutation_type_entry(self):
        for gene, position, aa in self.membership.keys() | self.types.keys():
            self.assert_official_annotation(gene, position, aa)

    def test_all_residues_and_indel_events_at_every_supported_position(self):
        # Both positive and negative membership over the complete protein
        # domain: 27,094 gene/position/residue-or-event combinations.
        for gene, reference in CONSENSUS_B_AMINO_ACIDS.items():
            for position in range(1, len(reference) + 1):
                for aa in AMINO_ACIDS + '*-_':
                    self.assert_official_annotation(gene, position, aa)

    def test_all_two_residue_mixtures_at_officially_annotated_positions(self):
        positions = {(gene, pos) for gene, pos, _ in self.membership.keys() | self.types.keys()}
        for gene, position in sorted(positions):
            for a, b in combinations(AMINO_ACIDS, 2):
                actual = annotate_drm(gene, position, a + b)
                reverse = annotate_drm(gene, position, b + a + a)
                self.assertEqual(actual, reverse, (gene, position, a, b))
                expected = [bool(self.membership.get((gene, position, aa))) for aa in (a, b)]
                self.assertIs(actual.is_drm, any(expected))
                self.assertEqual([c.aa for c in actual.components], sorted((a, b)))
                self.assertEqual([c.is_drm for c in actual.components], expected)
                for component in actual.components:
                    single = annotate_drm(gene, position, component.aa)
                    self.assertEqual(component, single.components[0])
                drugs = self.membership.get((gene, position, a), set()) | self.membership.get((gene, position, b), set())
                self.assertEqual(actual.drug_class, ';'.join(sorted(drugs)))

    def test_accessory_annotation_does_not_require_drm_membership(self):
        cases = [key for key, labels in self.types.items()
                 if 'Accessory' in labels and key not in self.membership]
        self.assertTrue(cases, 'The official snapshot includes non-DRM accessory variants')
        for key in cases:
            result = self.assert_official_annotation(*key)
            self.assertIs(result.is_drm, False)
            self.assertIs(result.is_accessory, True)

    def test_other_type_does_not_imply_non_drm(self):
        cases = [key for key in self.membership if self.types.get(key, {'Other'}) == {'Other'}]
        self.assertTrue(cases)
        for key in cases:
            self.assertIs(self.assert_official_annotation(*key).is_drm, True)

    def test_unresolved_status_overrides_every_known_drm(self):
        for gene, position, aa in self.membership:
            for status in UNRESOLVED_CODON_STATUSES:
                result = annotate_drm(gene, position, aa, codon_status=status + '|AMBIGUOUS')
                self.assertIsNone(result.is_drm)
                self.assertEqual(result.drm_status, 'UNRESOLVED')
                self.assertIsNone(result.is_accessory)
                self.assertEqual(result.drug_class, '')
                self.assertEqual(result.mutation_types, '')
                self.assertTrue(all(c.is_drm is None for c in result.components))
                self.assertEqual(result.drm_version, LOCAL_DRM_VERSION)
                self.assertEqual(len(result.drm_catalog_sha256), 64)

    def test_unknown_symbols_and_empty_input_are_unresolved(self):
        for aa in ('', 'X', 'B', 'Z', 'J', '?'):
            result = annotate_drm('PR', 1, aa)
            self.assertIsNone(result.is_drm)
            self.assertEqual(result.drm_status, 'UNRESOLVED')

    def test_unknown_mixture_component_preserves_resolved_component_evidence(self):
        gene, position, aa = next(key for key in self.membership if key[2] in AMINO_ACIDS)
        result = annotate_drm(gene, position, aa + 'X')
        self.assertIs(result.is_drm, True)
        self.assertEqual(result.drm_evaluation_status, 'PARTIAL')
        components = {c.aa: c for c in result.components}
        self.assertTrue(components[aa].is_drm)
        self.assertIsNone(components['X'].is_drm)

    def test_unresolved_mixture_without_confirmed_drm_remains_unresolved(self):
        result = annotate_drm('RT', 103, 'QX')
        self.assertIsNone(result.is_drm)
        self.assertEqual(result.drm_status, 'UNRESOLVED')
        self.assertEqual(result.drm_evaluation_status, 'PARTIAL')

    def test_accessory_mixture_uses_same_tri_state_logic(self):
        result = annotate_drm('PR', 10, 'FX')
        self.assertIs(result.is_drm, True)
        self.assertIs(result.is_accessory, True)
        self.assertEqual(result.drm_evaluation_status, 'PARTIAL')

    def test_ambiguous_resolved_mixture_is_evaluated_normally(self):
        for gene, position, aa in self.membership:
            if aa in AMINO_ACIDS:
                self.assertEqual(
                    annotate_drm(gene, position, aa, codon_status='AMBIGUOUS'),
                    annotate_drm(gene, position, aa),
                )

    def test_insertion_annotation_uses_event_not_inserted_amino_acids(self):
        for gene, position, aa in self.membership:
            if aa == '_':
                for inserted in ('A', 'VV', 'MK'):
                    actual = annotate_drm(gene, position, inserted, mutation_type='insertion')
                    self.assertEqual(actual, annotate_drm(gene, position, '_'))
                self.assertIsNone(annotate_drm(gene, position, 'X', mutation_type='insertion').is_drm)

    def test_provenance_and_hash_cover_both_official_datasets(self):
        metadata = json.loads((DATA / 'metadata.json').read_text())
        hashes = {name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
                  for name in ('drms_hiv1.json', 'mutation-type-pairs_hiv1.json')}
        self.assertEqual(hashes, metadata['sha256'])
        expected_hash = hashlib.sha256(json.dumps(hashes, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        for aa in ('A', 'X'):
            annotation = annotate_drm('PR', 1, aa)
            self.assertEqual(annotation.drm_source, LOCAL_DRM_SOURCE)
            self.assertEqual(annotation.drm_version, metadata['revision'])
            self.assertEqual(annotation.drm_catalog_sha256, expected_hash)

    def test_corruption_of_either_dataset_is_rejected(self):
        for name in ('drms_hiv1.json', 'mutation-type-pairs_hiv1.json'):
            with self.subTest(name=name), TemporaryDirectory() as directory:
                root = Path(directory)
                copytree(DATA, root / 'data')
                (root / 'data' / name).write_text('{}')
                load_catalog.cache_clear()
                try:
                    with patch('pyhiv.mutations.drm.files', return_value=root):
                        with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
                            annotate_drm('PR', 1, 'X')
                finally:
                    load_catalog.cache_clear()

    def test_incomplete_manifest_is_rejected(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            copytree(DATA, root / 'data')
            metadata = json.loads((root / 'data/metadata.json').read_text())
            metadata['sha256'].pop('mutation-type-pairs_hiv1.json')
            (root / 'data/metadata.json').write_text(json.dumps(metadata))
            load_catalog.cache_clear()
            try:
                with patch('pyhiv.mutations.drm.files', return_value=root):
                    with self.assertRaisesRegex(ValueError, 'both official datasets'):
                        load_catalog()
            finally:
                load_catalog.cache_clear()
