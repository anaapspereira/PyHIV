from pathlib import Path
from unittest import TestCase, mock

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pyhiv.align.align_with_reference as alignment


class TestResolveWorkerCount(TestCase):
    def test_resolve_worker_count_uses_explicit_n_jobs(self):
        self.assertEqual(alignment.resolve_worker_count(4), 4)

    def test_resolve_worker_count_falls_back_to_cpu_count(self):
        with mock.patch("pyhiv.align.align_with_reference.os.cpu_count", return_value=8):
            self.assertEqual(alignment.resolve_worker_count(None), 8)

    def test_resolve_worker_count_falls_back_to_one_when_cpu_count_unknown(self):
        with mock.patch("pyhiv.align.align_with_reference.os.cpu_count", return_value=None):
            self.assertEqual(alignment.resolve_worker_count(0), 1)


class TestCalculateAlignmentScore(TestCase):
    def test_valid_score(self):
        """Should count matching bases ignoring case and gaps."""
        seq1 = "ATGC"
        seq2 = "Atg-"
        result = alignment.calculate_alignment_score(seq1, seq2)
        self.assertEqual(result, 3)

    def test_handles_mismatched_lengths(self):
        """Should log error and return 0 when the two sequences have different lengths."""
        with mock.patch("logging.error") as mock_log:
            result = alignment.calculate_alignment_score("AAAA", "AAA")
            mock_log.assert_called_once()
            self.assertEqual(result, 0)


class TestProcessAlignment(TestCase):
    @mock.patch("pyhiv.align.align_with_reference.calculate_alignment_score", return_value=42)
    @mock.patch("pyhiv.align.align_with_reference.align_sequences", return_value=("AAA", "AAA"))
    def test_process_alignment_success(self, mock_align, mock_score):
        """Should return tuple with score, test, ref, and name."""
        test_seq = SeqRecord(Seq("AAA"), id="test")
        ref_seq = SeqRecord(Seq("AAA"), id="ref")
        ref_seq.name = "REF1"

        result = alignment.process_alignment(test_seq, ref_seq)
        self.assertEqual(result, (42, "AAA", "AAA", "REF1"))
        mock_align.assert_called_once_with(test_seq, ref_seq, alignment.DEFAULT_ALIGNMENT_TOOL)
        mock_score.assert_called_once()

    @mock.patch("logging.error")
    @mock.patch("pyhiv.align.align_with_reference.align_sequences", side_effect=Exception("boom"))
    def test_process_alignment_failure(self, mock_align, mock_log):
        """Should log error and return None when alignment fails."""
        test_seq = SeqRecord(Seq("AAA"), id="test")
        ref_seq = SeqRecord(Seq("AAA"), id="ref")
        ref_seq.name = "REF1"

        result = alignment.process_alignment(test_seq, ref_seq)
        self.assertIsNone(result)
        mock_log.assert_called_once()
        self.assertIn("Failed to process REF1", mock_log.call_args[0][0])

    @mock.patch("pyhiv.align.align_with_reference.calculate_alignment_score", return_value=42)
    @mock.patch("pyhiv.align.align_with_reference.align_sequences", return_value=("AAA", "AAA"))
    def test_process_alignment_uses_selected_tool(self, mock_align, mock_score):
        """Should pass the selected alignment tool to the dispatcher."""
        test_seq = SeqRecord(Seq("AAA"), id="test")
        ref_seq = SeqRecord(Seq("AAA"), id="ref")
        ref_seq.name = "REF1"

        result = alignment.process_alignment(test_seq, ref_seq, alignment_tool="PyFamsa")

        self.assertEqual(result, (42, "AAA", "AAA", "REF1"))
        mock_align.assert_called_once_with(test_seq, ref_seq, "PyFamsa")


class TestKmerShortlist(TestCase):
    def test_kmers_filters_ambiguous_bases_and_gaps(self):
        self.assertEqual(alignment.kmers("AAC-GTN", 3), {"AAC"})

    def test_kmers_rejects_non_positive_size(self):
        with self.assertRaises(ValueError):
            alignment.kmers("ACGT", 0)

    def test_kmers_returns_empty_set_for_short_sequence(self):
        self.assertEqual(alignment.kmers("AC", 5), set())

    def test_kmer_shortlist_falls_back_to_prefix_when_query_has_no_kmers(self):
        query = SeqRecord(Seq("NNNN"), id="query")
        references = [
            SeqRecord(Seq("AAAA"), id="ref1", name="ref1"),
            SeqRecord(Seq("CCCC"), id="ref2", name="ref2"),
            SeqRecord(Seq("GGGG"), id="ref3", name="ref3"),
        ]

        result = alignment.kmer_shortlist(query, references, size=4, top_k=2)

        self.assertEqual(result, references[:2])

    def test_kmer_shortlist_ranks_by_query_containment(self):
        query = SeqRecord(Seq("AAAACCCCGGGG"), id="query")
        ref_low = SeqRecord(Seq("TTTTAAAATTTT"), id="ref_low", name="ref_low")
        ref_high = SeqRecord(Seq("GGGGAAAACCCC"), id="ref_high", name="ref_high")
        ref_none = SeqRecord(Seq("TTTTTTTTTTTT"), id="ref_none", name="ref_none")

        result = alignment.kmer_shortlist(
            query,
            [ref_low, ref_none, ref_high],
            size=4,
            top_k=2,
        )

        self.assertEqual([record.name for record in result], ["ref_high", "ref_low"])

    def test_kmer_shortlist_can_be_disabled_with_non_positive_top_k(self):
        references = [
            SeqRecord(Seq("AAAA"), id="ref1", name="ref1"),
            SeqRecord(Seq("CCCC"), id="ref2", name="ref2"),
        ]

        self.assertEqual(
            alignment.kmer_shortlist(
                SeqRecord(Seq("AAAA"), id="query"),
                references,
                size=2,
                top_k=0,
            ),
            references,
        )

    def test_kmer_shortlist_uses_precomputed_reference_kmers(self):
        query = SeqRecord(Seq("AAAA"), id="query")
        ref1 = SeqRecord(Seq("CCCC"), id="ref1", name="ref1")
        ref2 = SeqRecord(Seq("GGGG"), id="ref2", name="ref2")

        result = alignment.kmer_shortlist(
            query,
            [ref1, ref2],
            size=4,
            top_k=1,
            reference_kmers=[{"AAAA"}, set()],
        )

        self.assertEqual(result, [ref1])

    def test_kmer_shortlist_rejects_mismatched_reference_kmers(self):
        with self.assertRaises(ValueError):
            alignment.kmer_shortlist(
                SeqRecord(Seq("AAAA"), id="query"),
                [
                    SeqRecord(Seq("AAAA"), id="ref1"),
                    SeqRecord(Seq("CCCC"), id="ref2"),
                ],
                size=4,
                top_k=1,
                reference_kmers=[],
            )


class TestAlignWithReferences(TestCase):

    def setUp(self):
        self.ref_dir = Path("tmp_refs")
        self.ref_dir.mkdir(exist_ok=True)
        self.test_seq = SeqRecord(Seq("AAA"), id="test")

    def tearDown(self):
        for f in self.ref_dir.glob("*"):
            f.unlink()
        self.ref_dir.rmdir()

    @mock.patch("logging.error")
    def test_invalid_reference_dir(self, mock_log):
        """Should return None and log error if reference dir is invalid."""
        result = alignment.align_with_references(self.test_seq, references_dir="not_a_path")
        self.assertIsNone(result)
        mock_log.assert_called_once_with("Invalid reference directory provided.")

    @mock.patch("logging.error")
    def test_reference_dir_does_not_exist(self, mock_log):
        """Should return None if directory path does not exist."""
        fake_path = self.ref_dir / "missing"
        result = alignment.align_with_references(self.test_seq, references_dir=fake_path)
        self.assertIsNone(result)
        mock_log.assert_called_once_with("Invalid reference directory provided.")

    @mock.patch("logging.error")
    def test_no_valid_references(self, mock_log):
        """Should log error if no FASTA files found."""
        result = alignment.align_with_references(self.test_seq, references_dir=self.ref_dir)
        self.assertIsNone(result)
        mock_log.assert_called_with("No valid reference sequences found.")

    @mock.patch("pyhiv.align.align_with_reference.as_completed")
    @mock.patch("pyhiv.align.align_with_reference.ProcessPoolExecutor")
    @mock.patch("pyhiv.align.align_with_reference.SeqIO.parse")
    def test_successful_alignment_with_best_score(self, mock_parse, mock_executor, mock_as_completed):
        """Test that the best alignment is selected correctly."""
        # Create fake FASTA files
        fasta1 = self.ref_dir / "ref1.fasta"
        fasta1.write_text(">ref1\nAAA\n")
        fasta2 = self.ref_dir / "ref2.fasta"
        fasta2.write_text(">ref2\nAAA\n")

        ref1 = SeqRecord(Seq("AAA"), id="ref1")
        ref2 = SeqRecord(Seq("AAA"), id="ref2")
        mock_parse.side_effect = [[ref1], [ref2]]

        # Dummy future objects
        class DummyFuture:
            def __init__(self, result_value):
                self._result_value = result_value

            def result(self):
                return self._result_value

        dummy_futures = [DummyFuture((5, "TEST", "REF", "refA")),
                         DummyFuture((10, "TEST2", "REF2", "refB"))]

        # Mock as_completed to just yield our dummy futures
        mock_as_completed.return_value = dummy_futures

        # Mock ProcessPoolExecutor to just return dummy futures mapping
        dummy_executor = mock.MagicMock()
        dummy_executor.__enter__.return_value = dummy_executor
        dummy_executor.__exit__.return_value = False
        dummy_executor.submit.side_effect = lambda *a, **kw: DummyFuture((10, "TEST2", "REF2", "refB"))
        mock_executor.return_value = dummy_executor

        result = alignment.align_with_references(self.test_seq, references_dir=self.ref_dir, n_jobs=1)

        self.assertEqual(result, ("TEST2", "REF2", "refB"))

    @mock.patch("pyhiv.align.align_with_reference.as_completed")
    @mock.patch("pyhiv.align.align_with_reference.ProcessPoolExecutor")
    @mock.patch("pyhiv.align.align_with_reference.SeqIO.parse")
    def test_alignment_scores_can_be_returned(self, mock_parse, mock_executor, mock_as_completed):
        fasta1 = self.ref_dir / "ref1.fasta"
        fasta1.write_text(">ref1\nAAA\n")
        fasta2 = self.ref_dir / "ref2.fasta"
        fasta2.write_text(">ref2\nAAA\n")

        ref1 = SeqRecord(Seq("AAA"), id="ref1")
        ref2 = SeqRecord(Seq("AAA"), id="ref2")
        mock_parse.side_effect = [[ref1], [ref2]]

        class DummyFuture:
            def __init__(self, result_value):
                self._result_value = result_value

            def result(self):
                return self._result_value

        dummy_futures = [
            DummyFuture((5, "TEST", "REF", "refA")),
            DummyFuture((10, "TEST2", "REF2", "refB")),
        ]

        dummy_executor = mock.MagicMock()
        dummy_executor.__enter__.return_value = dummy_executor
        dummy_executor.__exit__.return_value = False
        mock_executor.return_value = dummy_executor
        mock_as_completed.return_value = dummy_futures

        result = alignment.align_with_references(
            self.test_seq,
            references_dir=self.ref_dir,
            n_jobs=1,
            include_alignment_scores=True,
        )

        self.assertEqual(result, ("TEST2", "REF2", "refB", [(10, "refB"), (5, "refA")]))

    @mock.patch("pyhiv.align.align_with_reference.as_completed")
    @mock.patch("pyhiv.align.align_with_reference.ProcessPoolExecutor")
    def test_align_with_references_only_submits_kmer_shortlist(self, mock_executor, mock_as_completed):
        fasta1 = self.ref_dir / "ref1.fasta"
        fasta1.write_text(">ref1\nTTTTTTTTTT\n")
        fasta2 = self.ref_dir / "ref2.fasta"
        fasta2.write_text(">ref2\nAAAACCCCGG\n")
        fasta3 = self.ref_dir / "ref3.fasta"
        fasta3.write_text(">ref3\nCCCCGGGGTT\n")

        class DummyFuture:
            def __init__(self, result_value):
                self._result_value = result_value

            def result(self):
                return self._result_value

        dummy_executor = mock.MagicMock()
        dummy_executor.__enter__.return_value = dummy_executor
        dummy_executor.__exit__.return_value = False
        dummy_executor.submit.return_value = DummyFuture((10, "TEST", "REF", "ref2"))
        mock_executor.return_value = dummy_executor
        mock_as_completed.return_value = [dummy_executor.submit.return_value]

        result = alignment.align_with_references(
            SeqRecord(Seq("AAAACCCC"), id="test"),
            references_dir=self.ref_dir,
            n_jobs=1,
            kmer_size=4,
            reference_top_k=1,
        )

        self.assertEqual(result, ("TEST", "REF", "ref2"))
        self.assertEqual(dummy_executor.submit.call_count, 1)
        submitted_ref = dummy_executor.submit.call_args.args[2]
        self.assertEqual(submitted_ref.id, "ref2")

    @mock.patch("pyhiv.align.align_with_reference.SeqIO.parse", side_effect=Exception("bad parse"))
    @mock.patch("logging.error")
    def test_seqio_parse_error_logged(self, mock_log, mock_parse):
        """Should log error and skip bad FASTA file."""
        fasta_file = self.ref_dir / "bad.fasta"
        fasta_file.write_text(">bad\nAAA\n")

        result = alignment.align_with_references(self.test_seq, references_dir=self.ref_dir, n_jobs=1)
        self.assertIsNone(result)
        self.assertIn("No valid reference sequences found.", mock_log.call_args[0][0])

    @mock.patch("pyhiv.align.align_with_reference.as_completed")
    @mock.patch("pyhiv.align.align_with_reference.ProcessPoolExecutor")
    def test_align_with_references_filters_reference_accessions(self, mock_executor, mock_as_completed):
        (self.ref_dir / "keep-A1.fasta").write_text(">keep-A1\nAAA\n")
        (self.ref_dir / "skip-B.fasta").write_text(">skip-B\nAAA\n")

        class DummyFuture:
            def result(self):
                return (10, "TEST", "REF", "keep-A1")

        dummy_future = DummyFuture()
        dummy_executor = mock.MagicMock()
        dummy_executor.__enter__.return_value = dummy_executor
        dummy_executor.__exit__.return_value = False
        dummy_executor.submit.return_value = dummy_future
        mock_executor.return_value = dummy_executor
        mock_as_completed.return_value = [dummy_future]

        result = alignment.align_with_references(
            self.test_seq,
            references_dir=self.ref_dir,
            n_jobs=1,
            allowed_reference_accessions={"keep"},
        )

        self.assertEqual(result, ("TEST", "REF", "keep-A1"))
        self.assertEqual(dummy_executor.submit.call_count, 1)
        submitted_ref = dummy_executor.submit.call_args.args[2]
        self.assertEqual(submitted_ref.id, "keep-A1")

    @mock.patch("pyhiv.align.align_with_reference.kmer_shortlist", return_value=[])
    @mock.patch("logging.error")
    def test_no_candidates_selected_by_kmer_prefilter(self, mock_log, mock_shortlist):
        """Should log error and return None when the k-mer prefilter finds no candidates."""
        (self.ref_dir / "ref1.fasta").write_text(">ref1\nAAA\n")

        result = alignment.align_with_references(self.test_seq, references_dir=self.ref_dir)

        self.assertIsNone(result)
        mock_log.assert_called_with("No reference sequences selected by k-mer prefilter.")

    @mock.patch("pyhiv.align.align_with_reference.as_completed")
    @mock.patch("pyhiv.align.align_with_reference.ProcessPoolExecutor")
    def test_returns_none_when_no_alignment_results(self, mock_executor, mock_as_completed):
        """Should return None when every submitted alignment task fails."""
        (self.ref_dir / "ref1.fasta").write_text(">ref1\nAAA\n")

        class DummyFuture:
            def result(self):
                return None

        dummy_future = DummyFuture()
        dummy_executor = mock.MagicMock()
        dummy_executor.__enter__.return_value = dummy_executor
        dummy_executor.__exit__.return_value = False
        dummy_executor.submit.return_value = dummy_future
        mock_executor.return_value = dummy_executor
        mock_as_completed.return_value = [dummy_future]

        result = alignment.align_with_references(self.test_seq, references_dir=self.ref_dir, n_jobs=1)

        self.assertIsNone(result)

    @mock.patch("pyhiv.align.align_with_reference.as_completed")
    @mock.patch("pyhiv.align.align_with_reference.ProcessPoolExecutor")
    def test_reuses_caller_supplied_executor_instead_of_owning_one(self, mock_executor_cls, mock_as_completed):
        """Should submit tasks to a caller-supplied executor and not create its own pool."""
        (self.ref_dir / "ref1.fasta").write_text(">ref1\nAAA\n")

        class DummyFuture:
            def result(self):
                return (10, "TEST", "REF", "ref1")

        dummy_future = DummyFuture()
        shared_executor = mock.MagicMock()
        shared_executor.submit.return_value = dummy_future
        mock_as_completed.return_value = [dummy_future]

        result = alignment.align_with_references(
            self.test_seq,
            references_dir=self.ref_dir,
            executor=shared_executor,
        )

        self.assertEqual(result, ("TEST", "REF", "ref1"))
        shared_executor.submit.assert_called_once()
        mock_executor_cls.assert_not_called()
