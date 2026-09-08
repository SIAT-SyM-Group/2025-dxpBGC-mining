"""Offline regression tests. No NCBI requests or bioinformatics tool installations."""
import ast
import copy
import csv
import io
import json
import sqlite3
import subprocess
import sys
import urllib.request
from pathlib import Path

import nbformat
import pytest
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import workflow_support as workflow


BLAST_ROW = "q\tref|WP_000000001.1|\t90\t100\t10\t0\t1\t100\t1\t100\t1e-30\t200\n"


class FakeClock:
    def __init__(self):
        self.now = 1000.0

    def time(self):
        return self.now

    def sleep(self, seconds):
        self.now += seconds


class Response:
    def __init__(self, text, status=200, headers=None):
        self.text, self.status_code, self.headers = text, status, headers or {}

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}")


class FakeBlast:
    def __init__(self, clock, events):
        self.clock, self.events, self.calls = clock, iter(events), []

    def request(self, method, url, **kwargs):
        self.calls.append((method, self.clock.time(), kwargs))
        event = next(self.events)
        if isinstance(event, Exception):
            raise event
        return event if isinstance(event, Response) else Response(event)


def query(tmp_path):
    path = tmp_path / "query.faa"
    path.write_text(">q\nMACDEFGHIK\n")
    return path


def run_blast(tmp_path, clock, session, **kwargs):
    return workflow.run_blast(query(tmp_path), tmp_path / "blast", "test@example.org",
        session=session, clock=clock.time, sleep=clock.sleep, report=lambda _: None, **kwargs)


def region(path, gap=25000, domains=2, duplicate_domains=False, pks_at=True, versioned_pfams=True):
    """An explicitly synthetic antiSMASH-shaped GenBank fixture, not research data."""
    record = SeqRecord(Seq("A" * 60000), id="SYNTHETIC", name="SYNTHETIC",
                       description="Synthetic regression fixture")
    record.annotations.update(molecule_type="DNA", organism="Synthetic fixture")
    record.features = [SeqFeature(FeatureLocation(0, 60000), type="source", qualifiers={
        "note": (["PKS_AT; ketoacyl synthase"] if pks_at else ["ketoacyl synthase"])})]
    record.features += [SeqFeature(FeatureLocation(100, 1000), type="protocluster", qualifiers={
        "product": ["NRPS"], "core_location": ["[100:1000](+)"]}),
        SeqFeature(FeatureLocation(1000 + gap, 2000 + gap), type="protocluster", qualifiers={
        "product": ["hglE-KS"], "core_location": [f"[{1000 + gap}:{2000 + gap}](+)"]})]
    for number in range(domains):
        for pfam, name, offset in (("PF00501", "AMP-binding", 0), ("PF00668", "Condensation_LCL", 30)):
            start = 5000 + 100 * number + offset
            feature = SeqFeature(FeatureLocation(start, start + 24), type="PFAM_domain", qualifiers={
                "db_xref": [pfam + (".1" if versioned_pfams else "")], "aSDomain": [name],
                "protein_start": [str(number * 50 + offset)], "protein_end": [str(number * 50 + offset + 8)]})
            record.features.append(feature)
            if duplicate_domains:
                duplicate = copy.deepcopy(feature)
                duplicate.qualifiers["tool"] = ["clusterhmmer"]
                record.features.append(duplicate)
    # A translated CDS lets a local DB fixture model both scaffold and protein records.
    record.features.append(SeqFeature(FeatureLocation(6000, 6030), type="CDS",
                                      qualifiers={"locus_tag": ["synthetic_gene"], "translation": ["MABCDEFGHI"]}))
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(record, path, "genbank")
    return path


def test_credentials_precedence_and_optional_secrets():
    env = {"NCBI_EMAIL": "env@example.org", "NCBI_API_KEY": "env-secret"}
    assert workflow.resolve_credentials("form@example.org", environ=env) == ("form@example.org", "env-secret")
    assert workflow.resolve_credentials(environ=env) == ("env@example.org", "env-secret")
    assert workflow.resolve_credentials(environ={}, secret_getter=lambda key: env[key]) == ("env@example.org", "env-secret")
    with pytest.raises(ValueError, match="contact email"):
        workflow.resolve_credentials(environ={})


@pytest.mark.parametrize("bad", ["", ">a\n", ">\nMACD", ">a\nMACD\n>a other\nMKK", "prefix\n>a\nMACD", ">a\nMA-AC", ">a\nMA*AC", ">a\n1234"])
def test_invalid_fasta_does_not_modify_existing_file(tmp_path, bad):
    path = tmp_path / "protein.faa"
    path.write_text("original")
    with pytest.raises(ValueError):
        workflow.write_protein_fasta(bad, path)
    assert path.read_text() == "original"


def test_fasta_indented_headers_ambiguities_and_parent_creation(tmp_path):
    path = tmp_path / "nested" / "query.faa"
    records = workflow.write_protein_fasta("  >alpha note\nmacd\n >beta\nBXZJUO", path, width=3)
    assert records == [("alpha note", "MACD"), ("beta", "BXZJUO")]
    assert list(SeqIO.parse(path, "fasta"))[1].id == "beta"


def test_blast_validates_every_row():
    assert workflow.parse_blast_tabular("<pre>" + BLAST_ROW + "</pre>")[0][1] == "ref|WP_000000001.1|"
    with pytest.raises(ValueError):
        workflow.parse_blast_tabular(BLAST_ROW + BLAST_ROW.replace("1e-30", "nan"))
    with pytest.raises(ValueError):
        workflow.parse_blast_tabular("<html>temporarily unavailable</html>")


def test_blast_resume_without_second_submission_and_obey_intervals(tmp_path):
    clock = FakeClock()
    session = FakeBlast(clock, ["RID = TEST123\nRTOE = 10", "Status=WAITING",
                                "Status=READY\nThereAreHits=yes", BLAST_ROW])
    with pytest.raises(TimeoutError):
        run_blast(tmp_path, clock, session, max_wait_hours=0.02, poll_seconds=1)
    checkpoint = next((tmp_path / "blast").glob("*/job.json"))
    assert json.loads(checkpoint.read_text())["rid"] == "TEST123"
    output = run_blast(tmp_path, clock, session, poll_seconds=1)
    assert output.read_text() == BLAST_ROW
    assert [method for method, _, _ in session.calls].count("POST") == 1
    times = [at for _, at, _ in session.calls]
    assert all(b - a >= 10 for a, b in zip(times, times[1:]))
    polls = [at for _, at, kw in session.calls if kw.get("params", {}).get("FORMAT_OBJECT") == "SearchInfo"]
    assert polls[1] - polls[0] >= 60
    # An unchanged complete result does not contact NCBI at all.
    assert run_blast(tmp_path, clock, session) == output
    assert len(session.calls) == 4


def test_blast_ambiguous_post_is_not_automatically_retried(tmp_path):
    clock = FakeClock()
    session = FakeBlast(clock, [requests.Timeout("uncertain submission")])
    with pytest.raises(RuntimeError, match="checkpoint"):
        run_blast(tmp_path, clock, session)
    with pytest.raises(RuntimeError, match="no confirmed RID"):
        run_blast(tmp_path, clock, session)
    assert len(session.calls) == 1


def test_blast_get_retry_no_hits_and_resume_cache(tmp_path):
    clock = FakeClock()
    session = FakeBlast(clock, ["RID = EMPTY123", Response("busy", 429, {"Retry-After": "65"}),
                                "Status=READY\nThereAreHits=no"])
    path = run_blast(tmp_path, clock, session)
    assert path.read_text() == ""
    assert session.calls[-1][1] - session.calls[-2][1] >= 65
    assert json.loads((path.parent / "job.json").read_text())["no_hits"] is True
    assert run_blast(tmp_path, clock, session) == path
    assert len(session.calls) == 3


def test_blast_corrupt_output_refetches_existing_rid(tmp_path):
    clock = FakeClock()
    session = FakeBlast(clock, ["RID = RID123", "Status=READY\nThereAreHits=yes", BLAST_ROW,
                                "Status=READY\nThereAreHits=yes", BLAST_ROW])
    output = run_blast(tmp_path, clock, session)
    output.write_text("corrupt")
    assert run_blast(tmp_path, clock, session).read_text() == BLAST_ROW
    assert sum(call[0] == "POST" for call in session.calls) == 1


def test_blast_force_new_archives_previous_result(tmp_path):
    clock = FakeClock()
    session = FakeBlast(clock, ["RID = RID123", "Status=READY\nThereAreHits=yes", BLAST_ROW,
                                "RID = RID456", "Status=READY\nThereAreHits=no"])
    first = run_blast(tmp_path, clock, session)
    run_blast(tmp_path, clock, session, force_new=True)
    assert json.loads((first.parent / "job.json").read_text())["rid"] == "RID456"
    assert next(first.parent.glob("previous-*/hits.tsv")).read_text() == BLAST_ROW


def test_entrez_retries_are_limited_and_key_is_not_in_url(monkeypatch):
    clock, calls, events = FakeClock(), [], iter([Response("busy", 503), Response("ACGT")])
    client = workflow.EntrezClient("test@example.org", "secret")
    client.limiter = workflow.RateLimiter(0.11, clock.time, clock.sleep)

    class Session:
        def post(self, url, **kwargs):
            calls.append((url, clock.time(), kwargs))
            return next(events)

    client.local.session = Session()
    monkeypatch.setattr(workflow.time, "sleep", clock.sleep)
    assert client.fetch(db="nuccore", id="NC_123.1") == "ACGT"
    assert calls[1][1] - calls[0][1] >= 0.11
    assert all("secret" not in url for url, _, _ in calls)
    assert calls[0][2]["data"]["api_key"] == "secret"
    assert calls[0][2]["timeout"] == (15, 90)


def test_ipg_selection_keeps_order_refseq_and_missing_taxid_policy():
    xml = '<IPGReportSet><CDS accver="AB123.1" start="9" stop="3" strand="-" taxid="4" org="A"/><CDS accver="NZ_123.1" start="3" stop="9" strand="+" taxid="5" org="B"/></IPGReportSet>'
    mappings = workflow.parse_ipg(xml, "WP_123.1")
    assert workflow.choose_representative_mapping(mappings).accver == "NZ_123.1"
    assert workflow.choose_representative_mapping(mappings, False).accver == "AB123.1"
    assert mappings[0].start == 9 and mappings[0].strand == "-"


def test_ipg_rejects_an_html_error_page():
    with pytest.raises(ValueError, match="not an IPG"):
        workflow.parse_ipg("<html><body>Service unavailable</body></html>", "WP_123.1")


def test_download_is_atomic_idempotent_and_records_truncation(tmp_path):
    class Client:
        calls = 0

        def fetch(self, **kwargs):
            self.calls += 1
            assert kwargs["seq_start"] == 1 and kwargs["seq_stop"] == 19
            assert kwargs["strand"] == 1
            return ">NC_123.1:1-12\nACGTACGTACGT\n"

    client = Client()
    mapping = workflow.CdsMapping("WP_123.1", "NC_123.1", 9, 3, "-", "5", "Test organism")
    result, status = workflow.download_neighborhood(mapping, 10, tmp_path, client)
    again, status_again = workflow.download_neighborhood(mapping, 10, tmp_path, client)
    assert result == again and status_again == "cached" and client.calls == 1
    assert result["actual_stop"] == 12 and result["truncated_right"] is True
    assert len(list(SeqIO.parse(result["path"], "fasta"))) == 1
    Path(result["path"]).write_text("truncated file")
    workflow.download_neighborhood(mapping, 10, tmp_path, client)
    assert client.calls == 2
    assert len(list(SeqIO.parse(result["path"], "fasta"))) == 1


def test_partial_download_cannot_produce_success_manifest(tmp_path):
    class Client:
        api_key = ""
        fail = True

        def fetch(self, **kwargs):
            if kwargs["db"] == "protein":
                return '<IPGReportSet><CDS accver="NC_123.1" start="3" stop="9" strand="+" taxid="5" org="B"/></IPGReportSet>'
            if self.fail:
                raise RuntimeError("simulated download failure")
            return ">nucleotide\nACGTACGTACGT\n"

    client = Client()
    tsv = tmp_path / "hits.tsv"
    tsv.write_text(BLAST_ROW)
    with pytest.raises(RuntimeError, match="downloads failed"):
        workflow.mine_neighborhoods(tsv, tmp_path / "downloads", client, flank=10, report=lambda _: None)
    manifest = next((tmp_path / "downloads").glob("*/manifest.json"))
    assert json.loads(manifest.read_text())["status"] == "failed"
    with pytest.raises(ValueError, match="incomplete"):
        workflow.manifest_files(manifest)
    client.fail = False
    workflow.mine_neighborhoods(tsv, tmp_path / "downloads", client, flank=10, report=lambda _: None)
    assert len(workflow.manifest_files(manifest)) == 1


def test_resource_budget_prevents_cpu_oversubscription():
    assert workflow.resource_budget(20, 4, available_cpus=6, memory_gb=100) == (1, 4)
    assert workflow.resource_budget(5, 1, available_cpus=16, memory_gb=14) == (2, 1)
    with pytest.raises(ValueError):
        workflow.resource_budget(0, 1)


def test_database_download_toggle_does_not_disable_prerequisites_or_pfam(tmp_path, monkeypatch):
    calls = []
    db = tmp_path / "databases"
    db.mkdir()
    (db / "pfam.hmm").write_text("fixture")
    monkeypatch.setattr(workflow, "run_logged", lambda command, *args, **kwargs: calls.append(command))
    workflow.prepare_antismash_databases(["micromamba", "run"], db, tmp_path / "logs", download=False)
    assert len(calls) == 1 and "--check-prereqs" in calls[0] and str(db) in calls[0]
    for mode in ("full", "filter-required"):
        args = workflow.antismash_arguments(db, mode)
        assert "--fullhmmer" in args and "--clusterhmmer" in args and str(db) in args


def test_run_logged_keeps_full_log_and_reports_failure(tmp_path):
    log = tmp_path / "long.log"
    with pytest.raises(RuntimeError, match="exit 2"):
        workflow.run_logged([sys.executable, "-c", "for i in range(1000): print(i)\nraise SystemExit(2)"], log)
    assert len(log.read_text().splitlines()) == 1000


@pytest.mark.parametrize("gap,expected", [(24999, True), (25000, True), (25001, False)])
def test_filter_gap_boundary(tmp_path, gap, expected):
    diagnosis = workflow.diagnose_region(region(tmp_path / "x.region001.gbk", gap=gap))
    assert diagnosis["error"] == ""
    assert diagnosis["pass"] is expected
    assert diagnosis["min_core_gap"] == gap


def test_duplicate_domain_annotations_do_not_double_count(tmp_path):
    diagnosis = workflow.diagnose_region(region(tmp_path / "x.region001.gbk", domains=1, duplicate_domains=True))
    assert diagnosis["pf00501"] == 1 and diagnosis["condensation"] == 1
    assert diagnosis["pass"] is False and diagnosis["domain_gate"] is False


def test_versionless_pfam_and_qualifier_names_are_supported(tmp_path):
    assert workflow.diagnose_region(region(tmp_path / "x.region001.gbk", versioned_pfams=False))["pass"] is True


def test_compound_core_location_and_distinct_record_domains(tmp_path):
    path = region(tmp_path / "x.region001.gbk")
    record = SeqIO.read(path, "genbank")
    record.features[2].qualifiers["core_location"] = ["join{[59000:60000](+), [26000:27000](+)}"]
    SeqIO.write(record, path, "genbank")
    assert workflow.diagnose_region(path)["min_core_gap"] == 25000
    path2 = region(tmp_path / "y.region001.gbk", domains=1)
    one = SeqIO.read(path2, "genbank")
    two = copy.deepcopy(one)
    two.id = "SYNTHETIC_2"
    SeqIO.write([one, two], path2, "genbank")
    assert workflow.diagnose_region(path2)["pf00501"] == 2


def test_invalid_genbank_is_reported_as_error_not_a_failed_biological_gate(tmp_path):
    path = tmp_path / "x.region001.gbk"
    path.write_text("truncated")
    diagnosis = workflow.diagnose_region(path)
    assert diagnosis["pass"] is False and diagnosis["error"]
    assert diagnosis["failed_rules"] == "parse_or_read_error"


def test_region_csv_quotes_and_explicit_column(tmp_path):
    path = tmp_path / "hits.csv"
    path.write_text('Description,Organism\n"contains, comma",x.region001\n')
    assert workflow.extract_region_names(path, 2) == ["x.region001"]
    with pytest.raises(ValueError, match="column 1"):
        workflow.extract_region_names(path, 1)


def test_filter_export_matches_diagnosis_and_preserves_previous_run(tmp_path):
    root, dest = tmp_path / "antismash", tmp_path / "final"
    good = region(root / "sample" / "good.region001.gbk", gap=25000)
    bad = region(root / "sample" / "bad.region001.gbk", gap=25001)
    table = tmp_path / "abspres.tsv"
    table.write_text("Organism\tScaffold\ngood.region001\tSYNTHETIC\nbad.region001\tSYNTHETIC\n")
    result = workflow.filter_regions(table, [good, bad], root, dest, copy_mode="gbk", report=lambda _: None)
    assert result["kept"] == [str(good)]
    assert len(result["exported"]) == 1 and dest in Path(result["exported"][0]).parents
    with Path(result["diagnosis"]).open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["path"] for row in rows if row["pass"] == "True"} == set(result["kept"])
    old_export = Path(result["exported"][0])
    empty = workflow.filter_regions(table, [good, bad], root, dest, gap_max=24999, copy_mode="gbk", report=lambda _: None)
    assert not empty["kept"] and not empty["exported"]
    assert Path(empty["diagnosis"]).is_file() and old_export.is_file()
    assert result["directory"] != empty["directory"]
    assert Path(empty["directory"], "kept_region_files.abs.txt").read_text() == ""
    # Selecting an earlier cached export also updates the latest pointer.
    workflow.filter_regions(table, [good, bad], root, dest, copy_mode="gbk", report=lambda _: None)
    assert json.loads((dest / "latest.json").read_text())["directory"] == result["directory"]


def test_filter_header_only_table_is_valid_zero_result(tmp_path):
    root = tmp_path / "antismash"
    path = region(root / "x.region001.gbk")
    table = tmp_path / "empty.tsv"
    table.write_text("Organism\tScaffold\tStart\tEnd\tScore\n")
    result = workflow.filter_regions(table, [path], root, tmp_path / "out")
    assert result["counts"] == {"matched": 0, "kept": 0}


def test_filter_refuses_ambiguous_names_and_destination_inside_inputs(tmp_path):
    root = tmp_path / "antismash"
    paths = [region(root / name / "x.region001.gbk") for name in ("a", "b")]
    table = tmp_path / "hits.tsv"
    table.write_text("Organism\tScaffold\nx.region001\tx\n")
    with pytest.raises(ValueError, match="uniquely matched"):
        workflow.filter_regions(table, paths, root, tmp_path / "output")
    with pytest.raises(ValueError, match="outside"):
        workflow.filter_regions(table, paths[:1], root, root / "output")
    assert all(path.exists() for path in paths)


@pytest.mark.parametrize("unique,identity,coverage,hits", [(3, 20, 50, 100), (1, 101, 50, 100), (1, 20, -1, 100), (1, 20, 50, 0)])
def test_cblaster_parameters_fail_before_expensive_tools(tmp_path, unique, identity, coverage, hits):
    with pytest.raises(ValueError):
        workflow.check_cblaster_query(query(tmp_path), unique, 3, identity, coverage, hits)


def test_database_validation_rejects_silent_partial_parse(tmp_path):
    db = tmp_path / "regions.sqlite3"
    with sqlite3.connect(db) as conn:
        conn.execute("CREATE TABLE feature(feature_type TEXT, organism TEXT, sequence TEXT)")
        conn.executemany("INSERT INTO feature VALUES(?,?,?)", [("scaffold", "a.region001", "AAA"), ("gene", "a.region001", "MMM")])
    workflow.validate_cblaster_database(db, [Path("a.region001.gbk")])
    with pytest.raises(ValueError, match="incomplete"):
        workflow.validate_cblaster_database(db, [Path("a.region001.gbk"), Path("b.region001.gbk")])


def test_pipeline_cache_invalidation_and_failed_antismash_retry(tmp_path, monkeypatch):
    """Exercise stage handoff with real FASTA/GenBank/SQLite files and a mocked tool boundary."""
    fasta = query(tmp_path)
    calls = []
    fail_once = [True]

    def fake_tool(command, log, **kwargs):
        calls.append(command)
        if command[0] == "fake-antismash":
            out = Path(command[command.index("--output-dir") + 1])
            out.mkdir(parents=True)
            if fail_once[0]:
                fail_once[0] = False
                (out / "partial.tmp").write_text("interrupted")
                raise RuntimeError("simulated tool failure")
            (out / "index.html").write_text("<h1>Fixture report</h1>")
            region(out / "x.region001.gbk")
        elif "makedb" in command:
            base = Path(command[command.index("-n") + 1])
            files = Path(command[-1]).read_text().splitlines()
            with sqlite3.connect(base.with_suffix(".sqlite3")) as conn:
                conn.execute("CREATE TABLE feature(feature_type TEXT, organism TEXT, sequence TEXT)")
                for path in files:
                    conn.executemany("INSERT INTO feature VALUES(?,?,?)", [("scaffold", Path(path).stem, "AAA"), ("gene", Path(path).stem, "MMM")])
            base.with_suffix(".dmnd").write_bytes(b"fixture diamond db")
            base.with_suffix(".fasta").write_text(">1\nMMM\n")
        elif "search" in command:
            for flag, value in (("--session_file", '{"organisms":[{}]}'),
                                ("--output", "Fixture"), ("--binary", "Organism\tScaffold\nx.region001\tSYNTHETIC\n")):
                Path(command[command.index(flag) + 1]).write_text(value)
        elif "-c" in command and "plot.save_html" in command[command.index("-c") + 1]:
            Path(command[-1]).write_text("<html>synthetic plot</html>")

    monkeypatch.setattr(workflow, "run_logged", fake_tool)
    versions, database = {"antismash": "8.fixture", "cblaster": "1.4.2"}, {"inventory": "fixture"}
    args = workflow.antismash_arguments(tmp_path / "db")
    common = ([fasta], tmp_path / "antismash", ["fake-antismash"], args, versions, database)
    with pytest.raises(RuntimeError, match="jobs failed"):
        workflow.run_antismash_batch(*common, report=lambda _: None)
    manifest = workflow.run_antismash_batch(*common, report=lambda _: None)
    assert any((tmp_path / "antismash").rglob("partial.tmp"))
    workflow.run_antismash_batch(*common, report=lambda _: None)
    assert len(calls) == 2  # One failed attempt, one success, then no tool invocation.
    regions = workflow.manifest_files(manifest)
    db = workflow.build_cblaster_db(regions, tmp_path / "cblaster-db", ["fake-cblaster"], versions)
    assert workflow.build_cblaster_db(regions, tmp_path / "cblaster-db", ["fake-cblaster"], versions) == db
    assert len(calls) == 3
    outputs = workflow.run_cblaster(fasta, db, tmp_path / "search", ["fake-cblaster"], versions, unique=1)
    workflow.run_cblaster(fasta, db, tmp_path / "search", ["fake-cblaster"], versions, unique=1)
    assert len(calls) == 5  # Search plus independent plot export, then both cached.
    result = workflow.filter_regions(outputs["binary"], regions, tmp_path / "antismash", tmp_path / "final", copy_mode="gbk")
    assert result["counts"]["kept"] == 1
    old = result["directory"]
    # A changed GenBank invalidates the DB and downstream filter without deleting old exports.
    region(regions[0], gap=25001)
    new_db = workflow.build_cblaster_db(regions, tmp_path / "cblaster-db", ["fake-cblaster"], versions)
    assert new_db != db and len(calls) == 6
    new = workflow.filter_regions(outputs["binary"], regions, tmp_path / "antismash", tmp_path / "final", copy_mode="gbk")
    assert new["counts"]["kept"] == 0 and Path(old).is_dir()
    # A changed anchor reruns antiSMASH even though an earlier index.html exists.
    fasta.write_text(">changed\nMKLL\n")
    workflow.run_antismash_batch(*common, report=lambda _: None)
    assert len(calls) == 7


def test_preview_server_preserves_working_directory_and_closes_previous(tmp_path):
    original = Path.cwd()
    (tmp_path / "index.html").write_text("fixture")
    first = workflow.serve_directory(tmp_path)
    second = None
    try:
        with urllib.request.urlopen(f"http://127.0.0.1:{first.server_port}/index.html", timeout=5) as response:
            assert response.read() == b"fixture"
        second = workflow.serve_directory(tmp_path, first)
        assert first.socket.fileno() == -1
        assert Path.cwd() == original
    finally:
        target = second or first
        target.shutdown()
        target.server_close()


def test_cblaster_plot_failure_retries_without_repeating_search(tmp_path, monkeypatch):
    fasta = query(tmp_path)
    db = tmp_path / "regions.dmnd"
    db.write_text("fixture database")
    db.with_suffix(".sqlite3").write_text("fixture fingerprint input")
    calls = {"search": 0, "plot": 0}

    def fake_tool(command, *args, **kwargs):
        if "search" in command:
            calls["search"] += 1
            for flag, text in (("--session_file", '{"organisms":[{}]}'), ("--output", "summary"),
                               ("--binary", "Organism\tScaffold\na.region001\ta\n")):
                Path(command[command.index(flag) + 1]).write_text(text)
        else:
            calls["plot"] += 1
            if calls["plot"] == 1:
                raise RuntimeError("simulated plot failure")
            Path(command[-1]).write_text("<html>plot</html>")

    monkeypatch.setattr(workflow, "run_logged", fake_tool)
    args = (fasta, db, tmp_path / "out", ["fake-cblaster"], {"cblaster": "1.4.2"})
    with pytest.raises(RuntimeError, match="plot failure"):
        workflow.run_cblaster(*args, unique=1)
    result = workflow.run_cblaster(*args, unique=1)
    assert calls == {"search": 1, "plot": 2}
    assert Path(result["plot"]).is_file()


def test_notebook_schema_python_cells_outputs_and_embedded_source():
    root = Path(__file__).resolve().parents[1]
    notebook = nbformat.read(root / "workflow.ipynb", as_version=4)
    nbformat.validate(notebook)
    for cell in notebook.cells:
        if cell.cell_type == "code":
            ast.parse(cell.source)
            assert cell.outputs == [] and cell.execution_count is None
    assert "widgets" not in notebook.metadata
    subprocess.run([sys.executable, str(root / "scripts/sync_notebook.py"), "--check"], check=True)


def test_notebook_input_cells_execute_without_colab_or_remote_calls(tmp_path, monkeypatch):
    root = Path(__file__).resolve().parents[1]
    nb = nbformat.read(root / "workflow.ipynb", as_version=4)
    cells = {cell.id: cell.source for cell in nb.cells if cell.cell_type == "code"}
    namespace = {"__name__": "__main__"}
    monkeypatch.setenv("NCBI_EMAIL", "test@example.org")
    monkeypatch.setenv("NCBI_API_KEY", "")

    def run_cell(identifier, **overrides):
        tree = ast.parse(cells[identifier])
        for node in tree.body:
            if isinstance(node, ast.Assign) and len(node.targets) == 1:
                target = node.targets[0]
                if isinstance(target, ast.Name) and target.id in overrides:
                    node.value = ast.Constant(overrides[target.id])
        exec(compile(ast.fix_missing_locations(tree), identifier, "exec"), namespace)

    run_cell("workflow-support")
    run_cell("setup", WORK_DIR=str(tmp_path / "results"))
    run_cell("anchor", query_sequence=">q\nMACDEFGHIK\n")
    markers = tmp_path / "markers.faa"
    markers.write_text(">one\nMACD\n>two\nMKKK\n>three\nMNLL\n")
    run_cell("markers", QUERY_FASTA=str(markers))
    tsv = tmp_path / "prior.tsv"
    tsv.write_text(BLAST_ROW)
    namespace["run_blast"] = lambda *args, **kwargs: pytest.fail("Imported TSV must bypass remote BLAST")
    run_cell("blast", BLAST_TSV_INPUT=str(tsv))
    assert namespace["BLAST_TSV"].read_text() == BLAST_ROW
    assert namespace["MARKER_QUERY"].is_file()
    assert namespace["FNA_MANIFEST"] is None
    runtime = json.loads((tmp_path / "results/runtime.json").read_text())
    assert "api_key" not in runtime and "email" not in runtime
