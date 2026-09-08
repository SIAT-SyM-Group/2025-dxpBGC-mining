"""Display saved workflow evidence without rerunning or changing selection."""

from collections import Counter as _UICounter
from html import escape as _ui_escape
from pathlib import Path as _UIPath


_UI_GATES = {
    "has_pks_at": "PKS_AT annotation",
    "has_ketoacyl_synthase": "Ketoacyl synthase annotation",
    "domain_gate": "A / AMP-binding and condensation domains",
    "core_gap_gate": "NRPS–hglE-KS core distance",
    "has_exact_nrps": "Exact NRPS product",
    "parse_or_read_error": "File / annotation error",
}

_UI_STYLE = """
.dxp-report{font:15px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Arial,sans-serif;
color:#193247;background:#fff;border:1px solid #dce5ec;border-radius:16px;
max-width:1100px;box-sizing:border-box;overflow:hidden;text-align:left;margin:16px 0}
.dxp-report *{box-sizing:border-box}.dxp-report h2,.dxp-report h3,.dxp-report p{margin:0}
.dxp-report header{padding:26px 28px;background:#10283f;color:#fff}
.dxp-report header p{color:#bcd5df;font-size:13px;letter-spacing:1.4px;text-transform:uppercase}
.dxp-report h2{font-size:28px;line-height:1.3;font-weight:700;margin:7px 0}
.dxp-report header small{color:#dce7ed;font-size:14px}
.dxp-report .dxp-body{padding:24px 28px}
.dxp-report .dxp-metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(130px,1fr));
gap:12px;margin-bottom:26px}.dxp-report .dxp-metric{padding:16px;background:#f2f6f9;border-radius:10px}
.dxp-report .dxp-metric strong{display:block;font-size:32px;line-height:1.3;color:#10283f}
.dxp-report .dxp-metric span{display:block;font-size:13px;color:#526777}
.dxp-report .dxp-pass{background:#e6f4f0}.dxp-report .dxp-pass strong{color:#146b60}
.dxp-report h3{font-size:18px;margin:22px 0 6px}.dxp-report .dxp-note{font-size:13px;color:#526777}
.dxp-report .dxp-gates{display:grid;gap:12px;margin:17px 0 25px}
.dxp-report .dxp-gate-label{display:flex;justify-content:space-between;gap:16px;font-size:14px}
.dxp-report .dxp-track{height:7px;background:#edf1f5;border-radius:8px;margin-top:5px;overflow:hidden}
.dxp-report .dxp-bar{height:100%;background:#748d9e;border-radius:8px}
.dxp-report .dxp-scroll{overflow-x:auto;margin:15px 0;border:1px solid #dce5ec;border-radius:9px}
.dxp-report table{border-collapse:collapse;min-width:670px;width:100%;font-size:13px}
.dxp-report th{background:#f2f6f9;color:#405b6d;font-weight:600;text-align:left}
.dxp-report th,.dxp-report td{padding:12px 13px;border-bottom:1px solid #e7edf1;vertical-align:top}
.dxp-report tr:last-child td{border-bottom:0}.dxp-report td:first-child{max-width:240px;overflow-wrap:anywhere}
.dxp-report .dxp-status{display:inline-block;border-radius:5px;padding:2px 8px;white-space:nowrap;
color:#405b6d;background:#edf1f5;font-weight:600}.dxp-report .dxp-status-pass{background:#e6f4f0;color:#146b60}
.dxp-report details{margin-top:12px}.dxp-report summary{cursor:pointer;font-weight:600}
.dxp-report code{font-size:12px;overflow-wrap:anywhere;white-space:normal;color:#405b6d}
.dxp-report footer{padding:18px 28px;background:#f5f8fa;border-top:1px solid #dce5ec;font-size:13px;color:#526777}
.dxp-report .dxp-empty{padding:18px;background:#f2f6f9;border-radius:10px;margin:16px 0}
@media(max-width:540px){.dxp-report .dxp-body,.dxp-report header,.dxp-report footer{padding:18px}
.dxp-report h2{font-size:24px}.dxp-report .dxp-metrics{grid-template-columns:repeat(2,minmax(0,1fr))}}
"""


def diagnosis_dashboard(result, rows, limit=20):
    """Return standalone HTML derived only from saved diagnosis rows and config.

    No files are read or written here. Region names, paths and annotation text
    are escaped; the report needs no JavaScript, CDN, plotting library or GPU.
    """
    if limit < 1:
        raise ValueError("The table preview limit must be positive.")
    rows = list(rows)
    passed = sum(str(row.get("pass", "")).lower() == "true" for row in rows)
    matched = len(rows)
    gap = int(result["config"]["gap_max"])
    failures = _UICounter(gate for row in rows
                          for gate in row.get("failed_rules", "").split(";") if gate)
    e = lambda value: _ui_escape(str(value), quote=True)
    metrics = "".join(
        f'<div class="dxp-metric {accent}"><strong>{e(value)}</strong><span>{label}</span></div>'
        for value, label, accent in [
            (f"{matched:,}", "Matched regions", ""),
            (f"{passed:,}", "Passing regions", "dxp-pass"),
            (f"{matched - passed:,}", "Excluded regions", ""),
            (f"{gap:,}", "Core-gap cutoff · bp", ""),
        ]
    )
    gate_names = [*list(_UI_GATES)[:5], *sorted(set(failures) - set(list(_UI_GATES)[:5]))]
    gates = "".join(
        '<div><div class="dxp-gate-label">'
        f'<span>{e(_UI_GATES.get(gate, gate))}</span><strong>{failures[gate]:,} / {matched:,}</strong></div>'
        f'<div class="dxp-track" aria-hidden="true"><div class="dxp-bar" '
        f'style="width:{100 * failures[gate] / matched if matched else 0:.2f}%"></div></div></div>'
        for gate in gate_names
    )
    preview = []
    for row in rows[:limit]:
        is_pass = str(row.get("pass", "")).lower() == "true"
        status = "Pass" if is_pass else "Excluded"
        reasons = "; ".join(_UI_GATES.get(gate, gate)
                            for gate in row.get("failed_rules", "").split(";") if gate)
        domains = " / ".join(e(row.get(key, "—")) for key in ("pf00501", "amp_binding", "condensation"))
        core_gap = row.get("min_core_gap", "")
        preview.append(
            f'<tr><td title="{e(row.get("path", ""))}">{e(_UIPath(row.get("path", "")).name)}</td>'
            f'<td><span class="dxp-status {"dxp-status-pass" if is_pass else ""}">{status}</span></td>'
            f'<td>{domains}</td><td>{e(core_gap if core_gap != "" else "Not found")}</td>'
            f'<td>{e(reasons or "All gates passed")}</td></tr>'
        )
    table = (
        '<div class="dxp-scroll" tabindex="0" role="region" aria-label="Region evidence; scroll horizontally on small screens">'
        '<table><thead><tr><th scope="col">Region</th><th scope="col">Decision</th>'
        '<th scope="col">A / AMP / C</th><th scope="col">Gap · bp</th>'
        '<th scope="col">Failed conditions</th></tr></thead><tbody>' + "".join(preview) + '</tbody></table></div>'
        if rows else '<p class="dxp-empty">No matched regions in this run. Inspect the cblaster result and region-matching inputs.</p>'
    )
    dry = bool(result.get("config", {}).get("dry_run", False))
    export_note = "Preview only · files were not exported." if dry else "Use the passing list to locate the exported candidates."
    paths = "".join(f'<p>{label}<br><code>{e(path)}</code></p>' for label, path in [
        ("Complete evidence table", result["diagnosis"]),
        ("Authoritative passing list", _UIPath(result["directory"]) / "kept_region_files.abs.txt"),
    ])
    return (
        '<!doctype html><html lang="en"><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        '<title>DxpBGC · candidate evidence</title><style>' + _UI_STYLE + '</style>'
        '<section class="dxp-report" aria-label="DxpBGC candidate evidence report">'
        '<header><p>DxpBGC / Evidence review</p><h2>Candidate evidence</h2>'
        '<small>Read from the saved diagnosis · selection is unchanged</small></header>'
        '<div class="dxp-body"><div class="dxp-metrics">' + metrics + '</div>'
        '<h3>Where regions failed</h3><p class="dxp-note">Counts are independent: a region can fail more than one condition.</p>'
        '<div class="dxp-gates">' + gates + '</div><h3>Region evidence</h3>'
        f'<p class="dxp-note">Showing {min(limit, matched):,} of {matched:,} rows. '
        'A / AMP / C = PF00501 / AMP-binding / condensation feature counts.</p>' + table +
        '<details><summary>Result files</summary>' + paths + '</details></div>'
        '<footer>' + export_note + ' Computational candidates do not establish compound identity, bioactivity or experimental biosynthesis.</footer>'
        '</section></html>'
    )
