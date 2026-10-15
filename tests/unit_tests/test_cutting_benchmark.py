from scripts.compare_pyhiv_cuts_to_genbank import (
    CutMapping,
    ExpectedInterval,
    OriginalSequence,
    build_cut_mapping,
    build_expected_genbank_intervals,
    compare_predictions_to_expected,
    infer_gene_from_genbank_row,
    location_end_partial_or_uncertain,
    location_start_partial_or_uncertain,
    normalize_accession,
    summarize_boundary_aware,
    summarize_group,
    to_1_based_end,
    to_1_based_start,
)


def make_genbank_row(
    *,
    accession="AB123456.1",
    feature_index="1",
    feature_type="CDS",
    gene="gag",
    product="gag protein",
    note="",
    location=None,
    start=100,
    end=200,
    part_count="1",
    is_compound="False",
    partial_or_uncertain="False",
):
    return {
        "accession": accession,
        "feature_index": feature_index,
        "feature_type": feature_type,
        "gene": gene,
        "product": product,
        "note": note,
        "location": location if location is not None else f"[{start}:{end}](+)",
        "start_0_based": str(start),
        "end_0_based_exclusive": str(end),
        "part_count": part_count,
        "is_compound": is_compound,
        "partial_or_uncertain": partial_or_uncertain,
    }


def test_normalize_accession_removes_version_suffix():
    assert normalize_accession("AB123456.1") == "AB123456"
    assert normalize_accession("AB123456") == "AB123456"
    assert normalize_accession("K03455.2") == "K03455"


def test_coordinate_conversion_reports_genbank_style_coordinates():
    assert to_1_based_start(100) == 101
    assert to_1_based_end(200) == 200


def test_location_boundary_uncertainty_is_boundary_specific():
    assert location_start_partial_or_uncertain("[<2097:5121](+)") is True
    assert location_end_partial_or_uncertain("[<2097:5121](+)") is False
    assert location_start_partial_or_uncertain("[8392:>8936](+)") is False
    assert location_end_partial_or_uncertain("[8392:>8936](+)") is True


def test_gene_inference_prefers_primary_gene_qualifier():
    assert infer_gene_from_genbank_row(make_genbank_row(gene="gag")) == "gag"
    assert infer_gene_from_genbank_row(make_genbank_row(gene="pol")) == "pol"
    assert infer_gene_from_genbank_row(make_genbank_row(gene="env")) == "env"


def test_gene_inference_uses_explicit_product_rules_when_gene_missing():
    assert infer_gene_from_genbank_row(make_genbank_row(gene="", product="Gag polyprotein")) == "gag"
    assert infer_gene_from_genbank_row(make_genbank_row(gene="", product="envelope glycoprotein")) == "env"


def test_gene_inference_ignores_note_location_and_gag_pol_for_primary_regions():
    assert (
        infer_gene_from_genbank_row(
            make_genbank_row(gene="", product="", note="overlaps pol region")
        )
        == ""
    )
    assert infer_gene_from_genbank_row(make_genbank_row(gene="gag-pol", product="")) == ""


def test_genbank_feature_selection_single_valid_cds():
    expected, qc = build_expected_genbank_intervals(
        [make_genbank_row()],
        dataset_accessions=["AB123456"],
        regions=["gag"],
    )

    selected = expected[("AB123456", "gag")]
    assert selected.expected_start_0_based == 100
    assert selected.expected_end_0_based_exclusive == 200
    assert selected.feature_source == "cds_fallback"
    assert qc[0].selection_status == "cds_fallback"


def test_genbank_feature_selection_keeps_fuzzy_start_cds_but_counts_it():
    expected, qc = build_expected_genbank_intervals(
        [make_genbank_row(location="<100..500")],
        dataset_accessions=["AB123456"],
        regions=["gag"],
    )

    selected = expected[("AB123456", "gag")]
    assert selected.feature_source == "cds_fallback"
    assert qc[0].n_fuzzy_location == 1
    assert qc[0].fuzzy_location is True
    assert selected.start_partial_or_uncertain is True
    assert selected.end_partial_or_uncertain is False
    assert qc[0].selection_status == "cds_fallback"


def test_genbank_feature_selection_keeps_fuzzy_end_cds_but_counts_it():
    expected, qc = build_expected_genbank_intervals(
        [make_genbank_row(location="100..>500")],
        dataset_accessions=["AB123456"],
        regions=["gag"],
    )

    selected = expected[("AB123456", "gag")]
    assert selected.feature_source == "cds_fallback"
    assert qc[0].n_fuzzy_location == 1
    assert qc[0].fuzzy_location is True
    assert selected.start_partial_or_uncertain is False
    assert selected.end_partial_or_uncertain is True
    assert qc[0].selection_status == "cds_fallback"


def test_genbank_feature_selection_keeps_partial_cds_but_counts_it():
    expected, qc = build_expected_genbank_intervals(
        [make_genbank_row(partial_or_uncertain="True")],
        dataset_accessions=["AB123456"],
        regions=["gag"],
    )

    selected = expected[("AB123456", "gag")]
    assert selected.feature_source == "cds_fallback"
    assert qc[0].n_partial_or_uncertain == 1
    assert qc[0].partial_or_uncertain is True
    assert selected.start_partial_or_uncertain is True
    assert selected.end_partial_or_uncertain is True
    assert qc[0].selection_status == "cds_fallback"


def test_genbank_feature_selection_excludes_compound_cds():
    expected, qc = build_expected_genbank_intervals(
        [make_genbank_row(is_compound="True")],
        dataset_accessions=["AB123456"],
        regions=["gag"],
    )

    assert expected == {}
    assert qc[0].n_compound == 1
    assert qc[0].selection_status == "no_valid_gene"


def test_genbank_feature_selection_excludes_multipart_cds():
    expected, qc = build_expected_genbank_intervals(
        [make_genbank_row(part_count="2")],
        dataset_accessions=["AB123456"],
        regions=["gag"],
    )

    assert expected == {}
    assert qc[0].n_compound == 1
    assert qc[0].selection_status == "no_valid_gene"


def test_genbank_feature_selection_keeps_uncertain_cds_but_counts_it():
    expected, qc = build_expected_genbank_intervals(
        [make_genbank_row(partial_or_uncertain="True", start=99, end=200)],
        dataset_accessions=["AB123456"],
        regions=["gag"],
    )

    selected = expected[("AB123456", "gag")]
    assert selected.expected_start_0_based == 99
    assert qc[0].n_partial_or_uncertain == 1
    assert qc[0].selection_status == "cds_fallback"


def test_genbank_feature_selection_marks_two_valid_cds_as_ambiguous():
    expected, qc = build_expected_genbank_intervals(
        [
            make_genbank_row(feature_index="1"),
            make_genbank_row(feature_index="2", start=300, end=400),
        ],
        dataset_accessions=["AB123456"],
        regions=["gag"],
    )

    assert expected == {}
    assert qc[0].selection_status == "ambiguous_multiple_cds"


def test_genbank_feature_selection_prefers_gene_over_cds():
    expected, qc = build_expected_genbank_intervals(
        [
            make_genbank_row(feature_index="1", feature_type="gene", product="", start=90, end=210),
            make_genbank_row(feature_index="2", feature_type="CDS", start=100, end=200),
        ],
        dataset_accessions=["AB123456"],
        regions=["gag"],
    )

    selected = expected[("AB123456", "gag")]
    assert selected.feature_type == "gene"
    assert selected.expected_start_0_based == 90
    assert selected.feature_source == "gene"
    assert qc[0].selection_status == "valid_gene"


def test_genbank_feature_selection_uses_gene_by_default():
    expected, qc = build_expected_genbank_intervals(
        [make_genbank_row(feature_type="gene", product="")],
        dataset_accessions=["AB123456"],
        regions=["gag"],
    )

    selected = expected[("AB123456", "gag")]
    assert selected.feature_type == "gene"
    assert selected.feature_source == "gene"
    assert qc[0].selection_status == "valid_gene"


def test_genbank_feature_selection_cds_fallback_can_be_disabled():
    expected, qc = build_expected_genbank_intervals(
        [make_genbank_row()],
        dataset_accessions=["AB123456"],
        regions=["gag"],
        allow_cds_fallback=False,
    )

    assert expected == {}
    assert qc[0].selection_status == "no_valid_gene"


def test_exact_cut_mapping_recovers_start_and_end():
    originals = {
        "AB123456": OriginalSequence(
            fasta_file="AB123456.fasta",
            accession="AB123456",
            original_accession="AB123456.1",
            sequence_id="AB123456.1",
            sequence="NNNATGCNN",
        )
    }

    mapping = build_cut_mapping(
        tool="tool",
        result_dir="results",
        region_dir="gag",
        result_file="AB123456.1_gag.fasta",
        sequence_id="AB123456.1",
        raw_sequence="ATGC",
        originals=originals,
        max_map_error_rate=0.15,
    )

    assert mapping.match_status == "exact_sequence_match"
    assert mapping.start_0_based == 3
    assert mapping.end_0_based_exclusive == 7


def test_approximate_cut_mapping_recovers_position_with_small_difference():
    originals = {
        "AB123456": OriginalSequence(
            fasta_file="AB123456.fasta",
            accession="AB123456",
            original_accession="AB123456.1",
            sequence_id="AB123456.1",
            sequence="TTTTACGTACGTGGGG",
        )
    }

    mapping = build_cut_mapping(
        tool="tool",
        result_dir="results",
        region_dir="gag",
        result_file="AB123456.1_gag.fasta",
        sequence_id="AB123456.1",
        raw_sequence="ACGTTCGT",
        originals=originals,
        max_map_error_rate=0.25,
    )

    assert mapping.match_status == "approximate_match"
    assert mapping.start_0_based == 4
    assert mapping.end_0_based_exclusive == 12
    assert mapping.map_edit_distance == 1


def make_expected_interval(gene="pol"):
    return ExpectedInterval(
        accession="AB123456",
        original_accession="AB123456.1",
        gene=gene,
        expected_start_0_based=100,
        expected_end_0_based_exclusive=200,
        expected_length_nt=100,
        feature_type="CDS",
        feature_source="CDS",
        feature_index="1",
        qualifier_gene=gene,
        product="",
        note="",
        location="[100:200](+)",
        start_partial_or_uncertain=False,
        end_partial_or_uncertain=False,
    )


def make_cut_mapping(gene="pol"):
    return CutMapping(
        tool="tool1",
        result_dir="results",
        region_dir=gene,
        canonical_region=gene,
        result_file=f"AB123456.1_{gene}.fasta",
        accession="AB123456",
        original_accession="AB123456.1",
        sequence_id="AB123456.1",
        cut_length_raw=100,
        cut_length_ungapped=100,
        start_0_based=100,
        end_0_based_exclusive=200,
        match_count=1,
        match_status="exact_sequence_match",
        map_edit_distance=0,
        map_error_rate=0.0,
        map_normalized_edit_similarity=1.0,
    )


def test_multiple_outputs_do_not_count_as_successful_mapping_or_accuracy():
    expected = {("AB123456", "pol"): make_expected_interval()}
    details = compare_predictions_to_expected(
        expected_intervals=expected,
        mappings=[make_cut_mapping(), make_cut_mapping()],
        tools=["tool1"],
        tolerance=3,
    )

    assert len(details) == 1
    assert details[0]["match_status"] == "multiple_outputs"
    assert details[0]["mapping_successful"] is False
    assert details[0]["exact_boundary_match"] is False
    assert details[0]["within_tolerance"] is False

    summary = summarize_group("tool1", "pol", details)
    assert summary["n_mapped"] == 0
    assert summary["exact_boundary_n"] == 0
    assert summary["within_tolerance_n"] == 0
    assert summary["median_abs_start_error_nt"] == ""
    assert summary["median_abs_end_error_nt"] == ""
    assert summary["median_iou"] == ""


def test_summary_denominators_separate_coverage_from_accuracy():
    rows = []
    for index in range(10):
        mapped = index < 7
        exact = index < 5
        rows.append(
            {
                "output_present": index < 8,
                "mapping_successful": mapped,
                "exact_boundary_match": exact if mapped else False,
                "within_tolerance": exact if mapped else False,
                "abs_start_error_nt": 0 if mapped else "",
                "abs_end_error_nt": 0 if mapped else "",
                "abs_length_error_nt": 0 if mapped else "",
                "iou": 1.0 if mapped else "",
            }
        )

    summary = summarize_group("tool", "gag", rows)

    assert summary["n_expected"] == 10
    assert summary["n_output"] == 8
    assert summary["n_mapped"] == 7
    assert summary["coverage_pct"] == 70.0
    assert summary["exact_boundary_n"] == 5
    assert summary["exact_pct_of_mapped"] == 71.429
    assert summary["overall_exact_recovery_pct"] == 50.0


def test_boundary_aware_summary_excludes_only_fuzzy_boundary():
    expected = make_expected_interval()
    expected = ExpectedInterval(
        accession=expected.accession,
        original_accession=expected.original_accession,
        gene=expected.gene,
        expected_start_0_based=expected.expected_start_0_based,
        expected_end_0_based_exclusive=expected.expected_end_0_based_exclusive,
        expected_length_nt=expected.expected_length_nt,
        feature_type=expected.feature_type,
        feature_source=expected.feature_source,
        feature_index=expected.feature_index,
        qualifier_gene=expected.qualifier_gene,
        product=expected.product,
        note=expected.note,
        location="[<100:200](+)",
        start_partial_or_uncertain=True,
        end_partial_or_uncertain=False,
    )
    mapping = make_cut_mapping()
    details = compare_predictions_to_expected(
        expected_intervals={("AB123456", "pol"): expected},
        mappings=[mapping],
        tools=["tool1"],
        tolerance=3,
    )

    summary = summarize_boundary_aware(details, ["pol"])

    assert summary[0]["n_valid"] == 1
    assert summary[0]["n_fuzzy"] == 1
    assert summary[0]["start_comparable_n"] == 0
    assert summary[0]["end_comparable_n"] == 1
    assert summary[0]["end_exact_n"] == 1
    assert summary[0]["end_within_10nt_n"] == 1
