import copy

import pytest

from workflow_presentation import diagnosis_dashboard


def saved_result():
    return {"config": {"gap_max": 12345, "dry_run": False},
            "directory": "/example/run", "diagnosis": "/example/run/filter_diagnosis.tsv"}


def test_saved_evidence_is_not_reclassified_and_html_is_escaped():
    rows = [{"path": '/example/<script>alert(1)</script>.gbk', "pass": "True",
             "pf00501": "2", "amp_binding": "0", "condensation": "2",
             "min_core_gap": "25000", "failed_rules": ""},
            {"path": "/example/b.gbk", "pass": "False", "min_core_gap": "",
             "failed_rules": "domain_gate;core_gap_gate"}]
    before = copy.deepcopy(rows)
    report = diagnosis_dashboard(saved_result(), rows)
    assert '<script>' not in report and '&lt;script&gt;' in report
    assert '>12,345<' in report and '>25000<' in report
    assert 'dxp-status-pass">Pass' in report  # Preserve saved decision even when the displayed cutoff differs.
    assert '>Not found<' in report and 'width:50.00%' in report
    assert 'Showing 2 of 2 rows' in report and rows == before


def test_zero_matches_and_preview_do_not_claim_an_export():
    result = saved_result()
    result['config']['dry_run'] = True
    report = diagnosis_dashboard(result, [])
    assert 'No matched regions' in report
    assert 'files were not exported' in report
    assert 'width:0.00%' in report
    assert 'nan' not in report.lower()


def test_table_limit_preserves_complete_counts():
    rows = [{"path": f"/{i}.gbk", "pass": "True", "failed_rules": ""} for i in range(25)]
    report = diagnosis_dashboard(saved_result(), rows)
    assert 'Showing 20 of 25 rows' in report
    assert '>19.gbk<' in report and '>20.gbk<' not in report
    assert '>25<' in report
    with pytest.raises(ValueError, match='positive'):
        diagnosis_dashboard(saved_result(), rows, limit=0)
