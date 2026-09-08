"""Reusable, offline-testable helpers embedded into workflow.ipynb.

This module performs no installation, network access, or filesystem writes at
import time. Run scripts/sync_notebook.py after editing it.
"""

import csv
import hashlib
import html
import io
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import threading
import time
import uuid
import xml.etree.ElementTree as ET
from collections import Counter, deque
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

WORKFLOW_VERSION = "2.0.0"


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def workflow_source_sha():
    embedded = globals().get("_WORKFLOW_SOURCE_SHA")
    return embedded if embedded else sha256_file(__file__)


def atomic_write(path, content):
    """Replace a complete file; a failed write never leaves a partial target."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = "wb" if isinstance(content, bytes) else "w"
    kwargs = {} if mode == "wb" else {"encoding": "utf-8", "newline": ""}
    name = None
    try:
        with tempfile.NamedTemporaryFile(mode, dir=path.parent, delete=False, **kwargs) as handle:
            name = handle.name
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(name, path)
    finally:
        if name and Path(name).exists():
            Path(name).unlink()


def write_json(path, value):
    atomic_write(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def read_json(path, default=None):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except FileNotFoundError:
        return default


def write_table(path, rows, columns):
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, columns, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    atomic_write(path, buffer.getvalue())


def write_list(path, values):
    values = list(map(str, values))
    atomic_write(path, "\n".join(values) + ("\n" if values else ""))


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def sha256_file(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def file_signature(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": sha256_file(path), "size": path.stat().st_size}


def cache_valid(marker, fingerprint):
    try:
        state = read_json(marker, {})
        return (
            state.get("fingerprint") == fingerprint
            and state.get("status") == "complete"
            and bool(state.get("outputs"))
            and all(file_signature(item["path"]) == item for item in state["outputs"])
        )
    except (OSError, ValueError, KeyError, TypeError):
        return False


def complete_cache(marker, fingerprint, outputs, **metadata):
    state = dict(metadata, fingerprint=fingerprint, status="complete", completed_at=utc_now(),
                 outputs=[file_signature(path) for path in outputs], workflow_version=WORKFLOW_VERSION)
    write_json(marker, state)
    return state


def safe_name(value, length=100):
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("._")
    return value[:length] or "query"


def resolve_credentials(email="", api_key="", environ=None, secret_getter=None):
    """Form field > environment > optional Colab Secrets; never print secrets."""
    env = os.environ if environ is None else environ

    def resolve(value, key):
        if str(value or "").strip():
            return str(value).strip()
        if str(env.get(key, "")).strip():
            return str(env[key]).strip()
        if secret_getter is not None:
            try:
                return str(secret_getter(key) or "").strip()
            except Exception:  # Secrets are optional, including access not granted.
                pass
        return ""

    email = resolve(email, "NCBI_EMAIL")
    api_key = resolve(api_key, "NCBI_API_KEY")
    if not email or "@" not in email:
        raise ValueError("Set a contact email in NCBI_EMAIL (form, environment, or Colab Secrets).")
    return email, api_key


def protein_records(text, default_id="query"):
    """Validate all records before writing. Preserve sequence residues and IDs."""
    lines = [line.strip() for line in str(text).lstrip("\ufeff").splitlines() if line.strip()]
    if not lines:
        raise ValueError("Protein input is empty.")
    records = []
    if any(line.startswith(">") for line in lines):
        if not lines[0].startswith(">"):
            raise ValueError("Sequence text occurs before the first FASTA header.")
        header, seq = None, []
        for line in lines:
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header, seq = line[1:].strip(), []
            else:
                seq.append(re.sub(r"\s+", "", line).upper())
        records.append((header, "".join(seq)))
    else:
        records = [(safe_name(default_id), re.sub(r"\s+", "", "".join(lines)).upper())]
    seen = set()
    allowed = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
    for header, sequence in records:
        identifier = header.split()[0] if header else ""
        if not identifier or ">" in identifier:
            raise ValueError("Every FASTA record needs a non-empty identifier.")
        if identifier in seen:
            raise ValueError(f"Duplicate FASTA identifier: {identifier}")
        seen.add(identifier)
        if not sequence:
            raise ValueError(f"Empty protein sequence: {identifier}")
        invalid = sorted(set(sequence) - allowed)
        if invalid:
            raise ValueError(f"Invalid protein characters in {identifier}: {invalid}. Remove gaps/stop symbols explicitly.")
    return records


def write_protein_fasta(text, path, default_id="query", width=80):
    if not isinstance(width, int) or width < 1:
        raise ValueError("FASTA wrap width must be a positive integer.")
    records = protein_records(text, default_id)
    output = "".join(
        f">{header}\n" + "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width)) + "\n"
        for header, sequence in records
    )
    atomic_write(path, output)
    return records


class RateLimiter:
    """Thread-safe request spacing, including retries."""

    def __init__(self, interval, clock=time.monotonic, sleep=time.sleep):
        self.interval, self.clock, self.sleep = interval, clock, sleep
        self.lock, self.next_time = threading.Lock(), 0.0

    def wait(self):
        with self.lock:
            delay = self.next_time - self.clock()
            if delay > 0:
                self.sleep(delay)
            self.next_time = self.clock() + self.interval


def bounded_map(function, items, workers):
    """Keep at most 2 * workers tasks in memory; yield completed results."""
    if workers < 1:
        raise ValueError("workers must be positive")
    iterator = iter(items)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        pending = set()
        exhausted = False
        while pending or not exhausted:
            while not exhausted and len(pending) < 2 * workers:
                try:
                    item = next(iterator)
                except StopIteration:
                    exhausted = True
                else:
                    pending.add(executor.submit(function, item))
            if pending:
                done, pending = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    yield future.result()


def parse_blast_tabular(raw):
    """Reject malformed rows, login/error pages, and unexpected empty responses."""
    pre = re.search(r"<pre\b[^>]*>(.*?)</pre>", raw, re.I | re.S)
    text = html.unescape(pre.group(1) if pre else raw)
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 12 or not fields[0] or not fields[1]:
            raise ValueError("BLAST response is not a 12-column tabular result.")
        try:
            numbers = [float(value) for value in fields[2:]]
            if not all(math.isfinite(value) for value in numbers):
                raise ValueError
            if not 0 <= numbers[0] <= 100 or numbers[-2] < 0:
                raise ValueError
            for value in fields[3:10]:
                int(value)
        except ValueError as exc:
            raise ValueError("Invalid numeric field in BLAST tabular result.") from exc
        rows.append(fields)
    if not rows:
        raise ValueError("No tabular hits in BLAST response; only SearchInfo may confirm no hits.")
    return rows


def run_blast(query_fasta, output_root, email, database="refseq_protein", hit_limit=10,
              expect="1e-5", poll_seconds=60, max_wait_hours=6, force_new=False,
              session=None, clock=time.time, sleep=time.sleep, report=print):
    """Resume one matching RID. A submission is never automatically retried."""
    import requests

    if hit_limit < 1 or not math.isfinite(float(expect)) or float(expect) <= 0 or not math.isfinite(max_wait_hours) or max_wait_hours <= 0:
        raise ValueError("BLAST hit limit, expectation and waiting time must be positive.")
    protein_records(Path(query_fasta).read_text())
    config = {"query_sha256": sha256_file(query_fasta), "database": database,
              "program": "blastp", "hit_limit": hit_limit, "expect": str(expect)}
    fingerprint = digest(config)
    folder = Path(output_root) / fingerprint[:16]
    folder.mkdir(parents=True, exist_ok=True)
    checkpoint = folder / "job.json"
    server_clock = Path(output_root) / "last_contact.json"
    output = folder / "hits.tsv"
    state = read_json(checkpoint, {})
    if state and state.get("fingerprint") != fingerprint:
        raise RuntimeError("BLAST checkpoint does not match this query/configuration.")
    if state.get("status") == "complete" and not force_new:
        if state.get("output_sha256") == (sha256_file(output) if output.exists() else None):
            report(f"BLAST: reuse verified result ({state.get('hit_rows', 0)} rows): {output}")
            return output
    if force_new and state:
        # Retain the old RID and results for audit, including ambiguous submissions.
        archive = folder / ("previous-" + uuid.uuid4().hex[:12])
        archive.mkdir()
        for path in list(folder.iterdir()):
            if path.is_file():
                shutil.move(str(path), archive / path.name)
        state = {}
    if state.get("status") in {"submitting", "submission_unknown"} and not state.get("rid"):
        raise RuntimeError("Previous BLAST submission has no confirmed RID. Check job.json; set FORCE_NEW_BLAST only if resubmission is intended.")
    if state.get("status") in {"FAILED", "UNKNOWN"}:
        raise RuntimeError("The saved RID failed or expired. Set FORCE_NEW_BLAST to submit a new search.")
    session = session or requests.Session()
    url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    poll_seconds = max(60.0, float(poll_seconds))
    deadline = clock() + float(max_wait_hours) * 3600

    def save():
        write_json(checkpoint, state)

    def request(method, parameters, is_poll=False):
        # Persist wall-clock times so restarting the notebook cannot bypass limits.
        for attempt in range(4 if method == "GET" else 1):
            contact = read_json(server_clock, {})
            earliest = max(state.get("last_request_at", 0), contact.get("at", 0)) + 10
            if is_poll:
                earliest = max(earliest, state.get("last_poll_at", 0) + poll_seconds)
            delay = max(0, earliest - clock())
            if clock() + delay > deadline:
                raise TimeoutError(f"BLAST waiting budget reached; rerun to resume RID {state.get('rid', '')}.")
            if delay:
                sleep(delay)
            state["last_request_at"] = clock()
            write_json(server_clock, {"at": clock()})
            if is_poll:
                state["last_poll_at"] = clock()
            save()
            try:
                kwargs = {"params" if method == "GET" else "data": parameters, "timeout": (20, 300)}
                response = session.request(method, url, **kwargs)
                if response.status_code == 429 or response.status_code >= 500:
                    if method == "GET" and attempt < 3:
                        try:
                            retry_after = float(response.headers.get("Retry-After", 0))
                        except (TypeError, ValueError):
                            retry_after = 0
                        delay = max(10 * (2 ** attempt), retry_after)
                        if clock() + delay > deadline:
                            raise TimeoutError("BLAST retry would exceed waiting budget; rerun to resume.")
                        sleep(delay)
                        continue
                response.raise_for_status()
                return response.text
            except (requests.Timeout, requests.ConnectionError):
                if method != "GET" or attempt == 3:
                    raise RuntimeError("BLAST connection failed; the checkpoint was preserved. Rerun to resume a confirmed RID.") from None
        raise RuntimeError("BLAST request failed after retries.")

    if not state.get("rid"):
        state.update(fingerprint=fingerprint, config=config, status="submitting", started_at=utc_now())
        save()  # Record intent before POST so a disconnect cannot cause a duplicate job.
        try:
            raw = request("POST", {"CMD": "Put", "PROGRAM": "blastp", "DATABASE": database,
                "QUERY": Path(query_fasta).read_text(), "HITLIST_SIZE": str(hit_limit),
                "EXPECT": str(expect), "email": email, "tool": "symlab_dxpBGC_workflow"})
            atomic_write(folder / "submission.txt", raw)
            match = re.search(r"\bRID\s*=\s*([A-Z0-9-]+)", raw)
            if not match:
                raise RuntimeError("BLAST did not return a RID; inspect submission.txt before resubmitting.")
            state.update(rid=match.group(1), status="WAITING")
            rtoe = re.search(r"\bRTOE\s*=\s*(\d+)", raw)
            state["last_poll_at"] = clock() + max(0, int(rtoe.group(1)) - poll_seconds) if rtoe else clock()
            save()
        except Exception:
            if not state.get("rid"):
                state["status"] = "submission_unknown"
                save()
            raise
    report(f"BLAST RID {state['rid']} — rerun this cell to resume; checkpoint: {checkpoint}")
    while clock() < deadline:
        raw = request("GET", {"CMD": "Get", "RID": state["rid"], "FORMAT_OBJECT": "SearchInfo",
                              "email": email, "tool": "symlab_dxpBGC_workflow"}, is_poll=True)
        atomic_write(folder / "status.txt", raw)
        match = re.search(r"Status\s*=\s*(\w+)", raw)
        if not match:
            raise RuntimeError("Unrecognized BLAST status; inspect status.txt, then rerun to resume.")
        state["status"] = match.group(1).upper()
        save()
        if state["status"] in {"FAILED", "UNKNOWN"}:
            raise RuntimeError(f"BLAST RID {state['rid']} is {state['status']}; checkpoint preserved.")
        if state["status"] != "READY":
            report(f"BLAST: {state['status']}")
            continue
        no_hits = bool(re.search(r"ThereAreHits\s*=\s*no", raw, re.I))
        if no_hits:
            rows = []
        elif re.search(r"ThereAreHits\s*=\s*yes", raw, re.I):
            result = request("GET", {"CMD": "Get", "RID": state["rid"], "FORMAT_TYPE": "Text",
                "ALIGNMENT_VIEW": "Tabular", "DESCRIPTIONS": str(hit_limit), "ALIGNMENTS": str(hit_limit),
                "email": email, "tool": "symlab_dxpBGC_workflow"})
            atomic_write(folder / "results.raw.txt", result)
            rows = parse_blast_tabular(result)
        else:
            raise RuntimeError("BLAST READY response has no hit indicator; checkpoint preserved.")
        atomic_write(output, "".join("\t".join(row) + "\n" for row in rows))
        state.update(status="complete", hit_rows=len(rows), no_hits=no_hits,
                     completed_at=utc_now(), output_sha256=sha256_file(output))
        save()
        report(f"BLAST complete: {len(rows)} alignment rows → {output}")
        return output
    raise TimeoutError(f"BLAST waiting budget reached; rerun to resume RID {state['rid']}.")


class EntrezClient:
    """Small HTTP client with explicit timeouts and one shared request limiter."""

    def __init__(self, email, api_key="", timeout=(15, 90), attempts=4):
        self.email, self.api_key = email, api_key
        self.timeout, self.attempts = timeout, attempts
        self.limiter = RateLimiter(0.11 if api_key else 0.35)
        self.local = threading.local()

    def fetch(self, **parameters):
        import requests

        if not getattr(self.local, "session", None):
            self.local.session = requests.Session()
        payload = dict(parameters, email=self.email, tool="symlab_dxpBGC_workflow")
        if self.api_key:
            payload["api_key"] = self.api_key
        for attempt in range(self.attempts):
            self.limiter.wait()
            try:
                # POST keeps credentials out of URLs and exception request URLs.
                response = self.local.session.post("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                                                    data=payload, timeout=self.timeout)
                if response.status_code == 429 or response.status_code >= 500:
                    if attempt + 1 < self.attempts:
                        try:
                            retry_after = float(response.headers.get("Retry-After", 0))
                        except (TypeError, ValueError):
                            retry_after = 0
                        time.sleep(max(2 ** attempt, retry_after))
                        continue
                response.raise_for_status()
                if re.search(r"<ERROR>|Error occurred:", response.text, re.I):
                    raise RuntimeError("NCBI returned an error response for this accession.")
                return response.text
            except (requests.Timeout, requests.ConnectionError):
                if attempt + 1 == self.attempts:
                    break
                time.sleep(2 ** attempt)
        raise RuntimeError(f"NCBI request failed after {self.attempts} attempts.")


@dataclass(frozen=True)
class CdsMapping:
    wp: str
    accver: str
    start: int
    stop: int
    strand: str
    taxid: str = ""
    org: str = ""


def parse_ipg(xml_text, wp):
    root = ET.fromstring(xml_text)
    if root.tag not in {"IPGReportSet", "IPGReport"}:
        raise ValueError("NCBI response is not an IPG XML report.")
    errors = list(root.iter("ERROR"))
    if errors:
        raise ValueError("NCBI IPG report contains an error.")
    mappings = []
    for cds in root.iter("CDS"):
        fields = cds.attrib
        if not all(fields.get(key) for key in ("accver", "start", "stop")):
            continue
        start, stop = int(fields["start"]), int(fields["stop"])
        if min(start, stop) < 1:
            raise ValueError("Invalid IPG coordinates; expected one-based inclusive values.")
        mappings.append(CdsMapping(wp, fields["accver"], start, stop,
                                   fields.get("strand", "?"), fields.get("taxid", ""), fields.get("org", "")))
    return mappings


def choose_representative_mapping(mappings, prefer_refseq=True):
    """Keep the original ordered, one-representative-per-WP selection policy."""
    candidates = mappings
    if prefer_refseq:
        candidates = [item for item in mappings if item.accver.startswith(("NC_", "NZ_"))] or mappings
    return next((item for item in candidates if item.taxid.strip()), candidates[0] if candidates else None)


def download_neighborhood(mapping, flank, folder, client):
    from Bio import SeqIO

    start = max(1, min(mapping.start, mapping.stop) - flank)
    stop = max(mapping.start, mapping.stop) + flank
    config = dict(asdict(mapping), requested_start=start, requested_stop=stop, sequence_orientation="genomic-forward")
    fingerprint = digest(config)
    path = Path(folder) / f"{safe_name(mapping.accver)}_{start}_{stop}_{fingerprint[:12]}.fna"
    marker = path.with_suffix(".json")
    if cache_valid(marker, fingerprint):
        return read_json(marker)["record"], "cached"
    raw = client.fetch(db="nuccore", id=mapping.accver, rettype="fasta", retmode="text",
                       seq_start=start, seq_stop=stop, strand=1)
    record = SeqIO.read(io.StringIO(raw), "fasta")
    sequence = str(record.seq).upper()
    if not sequence or set(sequence) - set("ACGTRYSWKMBDHVN"):
        raise ValueError("NCBI returned an empty or invalid nucleotide sequence.")
    if len(sequence) > stop - start + 1:
        raise ValueError("Returned sequence is longer than the requested interval.")
    actual_stop = start + len(sequence) - 1
    if actual_stop < max(mapping.start, mapping.stop):
        raise ValueError("Returned interval does not contain the complete anchor CDS.")
    identifier = f"{safe_name(mapping.accver)}_{start}_{actual_stop}_{fingerprint[:10]}"
    description = f"WP={mapping.wp} taxid={mapping.taxid or 'NA'} strand=+ cds_strand={mapping.strand} {mapping.org}"
    fasta = f">{identifier} {description.strip()}\n" + "\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + "\n"
    atomic_write(path, fasta)
    result = dict(config, path=str(path.resolve()), actual_start=start, actual_stop=actual_stop,
                  length=len(sequence), truncated_right=actual_stop < stop)
    complete_cache(marker, fingerprint, [path], record=result)
    return result, "downloaded"


def mine_neighborhoods(blast_tsv, output_root, client, flank=100000, threads=10, top_taxids=1,
                       prefer_refseq=True, refresh_ipg=False, report=print):
    if flank < 0 or threads < 1 or top_taxids < 0:
        raise ValueError("FLANK/TOP_N_UNIQ_TAXID must be nonnegative and THREADS positive.")
    raw = Path(blast_tsv).read_text()
    if not raw.strip():
        raise ValueError("BLAST completed without hits; there are no neighborhoods to download.")
    rows = parse_blast_tabular(raw)
    root = Path(output_root)
    config = {"blast_sha256": sha256_file(blast_tsv), "flank": flank,
              "top_taxids": top_taxids, "prefer_refseq": prefer_refseq}
    folder = root / digest(config)[:16]
    folder.mkdir(parents=True, exist_ok=True)
    ipg_cache = root / "ipg-cache"
    picked, audit, seen_wp, seen_taxid = [], [], set(), set()
    for row in rows:
        hit = re.search(r"\b(WP_\d+\.\d+)\b", row[1])
        if not hit:
            audit.append({"subject": row[1], "status": "unsupported_accession", "error": "Only WP accessions are mapped by this workflow."})
            continue
        wp = hit.group(1)
        if wp in seen_wp:
            continue
        seen_wp.add(wp)
        item = {"subject": wp, "status": "", "error": ""}
        try:
            xml_path = ipg_cache / (wp + ".xml")
            if refresh_ipg or not xml_path.exists():
                xml = client.fetch(db="protein", id=wp, rettype="ipg", retmode="xml")
                mappings = parse_ipg(xml, wp)
                atomic_write(xml_path, xml)
            else:
                mappings = parse_ipg(xml_path.read_text(), wp)
            chosen = choose_representative_mapping(mappings, prefer_refseq)
            if not chosen:
                item["status"] = "no_cds_mapping"
            else:
                # Missing taxid retains the original per-WP fallback, explicitly reported.
                key = chosen.taxid or wp
                if key in seen_taxid:
                    item["status"] = "duplicate_taxid"
                else:
                    seen_taxid.add(key)
                    picked.append(chosen)
                    item.update(status="selected" if chosen.taxid else "selected_missing_taxid", **asdict(chosen))
        except Exception as exc:
            item.update(status="ipg_error", error=str(exc))
        audit.append(item)
        if top_taxids and len(picked) >= top_taxids:
            break
    write_table(folder / "selection.tsv", audit, ["subject", "status", "wp", "accver", "start", "stop", "strand", "taxid", "org", "error"])
    write_json(folder / "selection.json", {"config": config, "selected": [asdict(item) for item in picked],
                                           "audit": audit, "created_at": utc_now()})
    if any(item["status"] == "ipg_error" for item in audit):
        raise RuntimeError(f"IPG mapping errors are listed in {folder / 'selection.tsv'}. Rerun to retry; incomplete selection is not used downstream.")
    if not picked:
        raise ValueError(f"No eligible WP/CDS mappings; inspect {folder / 'selection.tsv'}.")
    records, failures, counts = [], [], Counter()

    def download(item):
        try:
            record, status = download_neighborhood(item, flank, folder / "fasta", client)
            return record, status, None
        except Exception as exc:
            return None, "failed", dict(asdict(item), error=str(exc))

    for record, status, error in bounded_map(download, picked, min(threads, 8 if client.api_key else 3)):
        counts[status] += 1
        if record:
            records.append(record)
        if error:
            failures.append(error)
        report(f"Neighborhoods: {sum(counts.values())}/{len(picked)} — {dict(counts)}")
    records.sort(key=lambda item: item["path"])
    write_table(folder / "download_errors.tsv", failures, ["wp", "accver", "start", "stop", "error"])
    state = {"status": "failed" if failures else "complete", "config": config, "records": records,
             "outputs": [file_signature(item["path"]) for item in records], "updated_at": utc_now()}
    write_json(folder / "manifest.json", state)
    if failures:
        raise RuntimeError(f"{len(failures)} neighborhood downloads failed. Rerun to retry; completed files will be verified and reused.")
    return folder / "manifest.json"


def manifest_files(path):
    state = read_json(path)
    if not state or state.get("status") != "complete":
        raise ValueError(f"The upstream step is incomplete: {path}")
    signatures = state.get("outputs", [])
    if not signatures and state.get("regions") == []:
        raise ValueError("antiSMASH completed but detected no regions. Review its reports before running cblaster.")
    if not signatures or any(file_signature(item["path"]) != item for item in signatures):
        raise ValueError(f"An upstream file is missing or changed. Rerun that step: {path}")
    return [Path(item["path"]) for item in signatures]


def run_logged(command, log_path, stream=False, env=None, cwd=None, report=print):
    """Stream to disk, retain only a short error tail, and stop owned children on interrupt."""
    import signal

    path = Path(log_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tail = deque(maxlen=30)
    with path.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(list(map(str, command)), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, errors="replace", env=env, cwd=cwd, bufsize=1, start_new_session=True)
        try:
            for line in process.stdout:
                log.write(line)
                log.flush()
                tail.append(line.rstrip())
                if stream:
                    report(line.rstrip())
            code = process.wait()
        except BaseException:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
            raise
        finally:
            process.stdout.close()
    if code:
        raise RuntimeError(f"Command failed (exit {code}); log: {path}\n" + "\n".join(tail))
    return "\n".join(tail)


def cpu_count():
    return len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)


def available_memory_gb():
    try:
        match = re.search(r"MemAvailable:\s+(\d+)", Path("/proc/meminfo").read_text())
        return int(match.group(1)) / (1024 ** 2) if match else None
    except OSError:
        return None


def resource_budget(max_jobs, cpus_per_job, memory_per_job_gb=6, available_cpus=None, memory_gb=None):
    if max_jobs < 1 or cpus_per_job < 1 or memory_per_job_gb <= 0:
        raise ValueError("Job count, CPUs and estimated memory per job must be positive.")
    available_cpus = cpu_count() if available_cpus is None else available_cpus
    memory_gb = available_memory_gb() if memory_gb is None else memory_gb
    cpus = min(cpus_per_job, max(1, available_cpus))
    jobs = min(max_jobs, max(1, available_cpus // cpus))
    if memory_gb is not None:
        jobs = min(jobs, max(1, int(max(0, memory_gb - 2) // memory_per_job_gb)))
    return jobs, cpus


def database_inventory(folder):
    """Fast database invalidation using relative path, size and mtime, not multi-GB rehashing."""
    folder = Path(folder).resolve()
    records = []
    for path in sorted(folder.rglob("*")):
        if path.is_file():
            stat = path.stat()
            records.append([str(path.relative_to(folder)), stat.st_size, stat.st_mtime_ns])
    if not records:
        raise ValueError(f"Database folder is empty: {folder}")
    return {"path": str(folder), "inventory_sha256": digest(records), "files": len(records),
            "bytes": sum(item[1] for item in records)}


def ensure_toolchain(tool_root, log_root, report=print):
    """Install only missing tools in a dedicated Linux environment; record the exact solution."""
    import platform
    import tarfile
    import requests

    if platform.system() != "Linux" or platform.machine() not in {"x86_64", "amd64"}:
        raise RuntimeError("Automatic tool installation targets Linux x86_64 (Colab). On other platforms use a Linux runtime/container.")
    root, logs = Path(tool_root).resolve(), Path(log_root).resolve()
    mm, env = root / "bin/micromamba", root / "envs/dxpbgc-v2"
    root.mkdir(parents=True, exist_ok=True)
    if not mm.is_file():
        report("Installing micromamba from its official distribution…")
        response = requests.get("https://micro.mamba.pm/api/micromamba/linux-64/latest", timeout=(20, 180))
        response.raise_for_status()
        with tarfile.open(fileobj=io.BytesIO(response.content), mode="r:bz2") as archive:
            member = archive.extractfile("bin/micromamba")
            if member is None:
                raise RuntimeError("micromamba archive does not contain bin/micromamba")
            atomic_write(mm, member.read())
        mm.chmod(0o755)
    os.environ["MAMBA_ROOT_PREFIX"] = str(root)
    binaries = ["antismash", "cblaster", "diamond", "blastp", "makeblastdb", "hmmscan", "hmmsearch",
                "hmmpress", "hmmpfam2", "fasttree", "prodigal"]
    if any(not (env / "bin" / name).is_file() for name in binaries):
        action = "install" if (env / "conda-meta").is_dir() else "create"
        report("Preparing antiSMASH 8 and cblaster 1.4.2; dependency solving runs only when tools are missing.")
        run_logged([str(mm), action, "-y", "-p", str(env), "-c", "conda-forge", "-c", "bioconda",
                    "python=3.11", "antismash=8", "cblaster=1.4.2", "hmmer", "hmmer2", "diamond",
                    "blast", "fasttree", "prodigal"], logs / "tool-install.log", stream=True, report=report)
    prefix = [str(mm), "run", "-p", str(env)]
    versions = {}
    for name, flags in (("antismash", ["--version"]), ("cblaster", ["--version"]),
                        ("diamond", ["version"]), ("blastp", ["-version"]), ("prodigal", ["-v"])):
        value = subprocess.check_output([*prefix, name, *flags], text=True, stderr=subprocess.STDOUT).strip()
        versions[name] = value
    if not re.search(r"\b8\.", versions["antismash"]):
        raise RuntimeError("This workflow requires antiSMASH 8. Use the dedicated tool environment.")
    explicit = subprocess.check_output([str(mm), "list", "-p", str(env), "--explicit"], text=True)
    atomic_write(logs / "environment-explicit.txt", explicit)
    versions["environment_sha256"] = hashlib.sha256(explicit.encode()).hexdigest()
    write_json(logs / "tool-versions.json", versions)
    return prefix, versions


def prepare_antismash_databases(command_prefix, database_dir, log_root, download=True, analysis_args=None):
    folder, logs = Path(database_dir).resolve(), Path(log_root).resolve()
    if download:
        folder.mkdir(parents=True, exist_ok=True)
        help_text = subprocess.check_output([*command_prefix, "download-antismash-databases", "--help"], text=True)
        flag = next((flag for flag in ("--database-dir", "--install-dir") if flag in help_text), None)
        if not flag:
            raise RuntimeError("Unrecognized database downloader CLI; specify a prepared DBDIR and set DOWNLOAD_DB=False.")
        run_logged([*command_prefix, "download-antismash-databases", flag, str(folder)],
                   logs / "database-download.log", stream=True)
    if not folder.is_dir():
        raise ValueError("Database directory does not exist. Enable DOWNLOAD_DB for the first run.")
    # No deletion or symlinking of package-managed data. Surface a genuine installation failure.
    analysis_args = antismash_arguments(folder) if analysis_args is None else analysis_args
    run_logged([*command_prefix, "antismash", "--check-prereqs", *analysis_args],
               logs / "prerequisites.log", stream=True)
    inventory = database_inventory(folder)
    write_json(logs / "database-inventory.json", inventory)
    return inventory


def antismash_arguments(database_dir, annotation_mode="full"):
    # Database download/reuse is deliberately independent of analysis modules.
    args = ["--taxon", "bacteria", "--genefinding-tool", "prodigal", "--allow-long-headers",
            "--databases", str(Path(database_dir).resolve()), "--fullhmmer", "--clusterhmmer"]
    if annotation_mode == "full":
        args += ["--tfbs", "--cb-general", "--cb-knownclusters", "--cb-subclusters", "--rre",
                 "--cc-mibig", "--tigrfam", "--asf", "--pfam2go", "--smcog-trees"]
    elif annotation_mode != "filter-required":
        raise ValueError("ANNOTATION_MODE must be full or filter-required.")
    return args


def run_antismash_batch(inputs, output_root, command_prefix, common_args, versions, database,
                        max_jobs=2, cpus_per_job=1, memory_per_job_gb=6, stream=False, report=print):
    inputs = sorted(set(Path(path).resolve() for path in inputs))
    if not inputs:
        raise ValueError("No FASTA inputs supplied to antiSMASH.")
    jobs, cpus = resource_budget(max_jobs, cpus_per_job, memory_per_job_gb)
    report(f"antiSMASH resources: {jobs} concurrent jobs × {cpus} CPU(s); memory budget is an estimate.")
    root = Path(output_root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    results = []
    env = dict(os.environ, MPLBACKEND="agg")

    def run_one(path):
        config = {"input": file_signature(path), "arguments": common_args, "cpus": cpus,
                  "versions": versions, "database": database, "workflow_version": WORKFLOW_VERSION,
                  "workflow_source_sha256": workflow_source_sha()}
        fingerprint = digest(config)
        folder = root / (safe_name(path.stem, 85) + "-" + fingerprint[:16])
        folder.mkdir(parents=True, exist_ok=True)
        marker = folder / "complete.json"
        if cache_valid(marker, fingerprint):
            state = read_json(marker)
            report(f"antiSMASH cached: {path.name}")
            return {"input": str(path), "status": "cached", "regions": state["regions"],
                    "report": state["report"], "marker": str(marker)}
        # Fresh attempt directory avoids both overwriting old results and reusing a partial run.
        attempt = folder / ("attempt-" + uuid.uuid4().hex[:12])
        log = folder / (attempt.name + ".log")
        try:
            run_logged([*command_prefix, *common_args, "-c", str(cpus), "--output-dir", str(attempt), str(path)],
                       log, stream=stream, env=env, report=report)
            index = attempt / "index.html"
            if not index.is_file():
                raise RuntimeError("antiSMASH exited without an HTML report.")
            regions = sorted(p for p in attempt.rglob("*.gbk") if re.search(r"\.region\d+\.gbk$", p.name))
            complete_cache(marker, fingerprint, [index, *regions], config=config,
                           regions=[str(p) for p in regions], report=str(index))
            report(f"antiSMASH complete: {path.name} ({len(regions)} regions)")
            return {"input": str(path), "status": "complete", "regions": list(map(str, regions)),
                    "report": str(index), "marker": str(marker)}
        except Exception as exc:
            return {"input": str(path), "status": "failed", "regions": [], "error": str(exc), "log": str(log)}

    for result in bounded_map(run_one, inputs, jobs):
        results.append(result)
    results.sort(key=lambda item: item["input"])
    failures = [item for item in results if item["status"] == "failed"]
    batch_key = digest([file_signature(p) for p in inputs] + [common_args, versions, database])[:16]
    manifest = root / ("batch-" + batch_key + ".json")
    regions = sorted({path for item in results for path in item["regions"]})
    write_json(manifest, {"status": "failed" if failures else "complete", "results": results,
                         "regions": regions, "outputs": [file_signature(path) for path in regions],
                         "updated_at": utc_now()})
    write_table(root / "last_batch.tsv", results, ["input", "status", "report", "log", "error"])
    if failures:
        raise RuntimeError(f"{len(failures)} antiSMASH jobs failed. See {root / 'last_batch.tsv'}; rerun to retry only failed jobs.")
    return manifest


def check_cblaster_query(query_fasta, unique, minimum_hits, identity, coverage, hit_limit):
    records = protein_records(Path(query_fasta).read_text())
    if not 1 <= unique <= len(records):
        raise ValueError(f"U={unique} requires at least {unique} unique marker sequences; query contains {len(records)}. Supply the marker set or deliberately change U.")
    if minimum_hits < 1 or not 0 <= identity <= 100 or not 0 <= coverage <= 100 or hit_limit < 1:
        raise ValueError("MH/HS must be positive; MI/MC must be between 0 and 100.")
    return records


def validate_cblaster_database(database, region_paths):
    """cblaster 1.4 may catch parser errors without failing its CLI: verify every source."""
    import sqlite3

    with sqlite3.connect(f"file:{Path(database).resolve()}?mode=ro", uri=True) as connection:
        if connection.execute("PRAGMA quick_check").fetchone()[0] != "ok":
            raise ValueError("cblaster SQLite integrity check failed.")
        organisms = {row[0] for row in connection.execute("SELECT DISTINCT organism FROM feature WHERE feature_type='scaffold'")}
        gene_count = connection.execute("SELECT count(*) FROM feature WHERE feature_type='gene' AND sequence IS NOT NULL").fetchone()[0]
    expected = {Path(path).stem for path in region_paths}
    if organisms != expected or not gene_count:
        raise ValueError(f"cblaster database is incomplete: {len(organisms)}/{len(expected)} source regions, {gene_count} translated genes.")


def build_cblaster_db(region_paths, output_root, command_prefix, versions, cpus=2, batch_size=32, report=print):
    paths = sorted(set(Path(path).resolve() for path in region_paths))
    if not paths or cpus < 1 or batch_size < 1:
        raise ValueError("Provide region GBKs, positive CPUs and a positive database batch size.")
    names = [path.stem for path in paths]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate region filenames would merge different cblaster organisms. Use unique FASTA record IDs and regenerate antiSMASH outputs.")
    config = {"inputs": [file_signature(path) for path in paths], "versions": versions,
              "workflow_version": WORKFLOW_VERSION, "workflow_source_sha256": workflow_source_sha()}
    fingerprint = digest(config)
    folder = Path(output_root).resolve() / fingerprint[:16]
    folder.mkdir(parents=True, exist_ok=True)
    marker = folder / "complete.json"
    if cache_valid(marker, fingerprint):
        state = read_json(marker)
        report(f"cblaster: reuse verified database ({len(paths)} regions)")
        return Path(state["database"])
    attempt = folder / ("attempt-" + uuid.uuid4().hex[:12])
    attempt.mkdir()
    # Native file-list input avoids shell argument limits and preserves original region names.
    file_list = attempt / "regions.txt"
    write_list(file_list, paths)
    database = attempt / "regions"
    run_logged([*command_prefix, "makedb", "-n", str(database), "--cpus", str(min(cpus, cpu_count())),
                "--batch", str(batch_size), str(file_list)], attempt / "makedb.log", report=report)
    outputs = [database.with_suffix(ext) for ext in (".dmnd", ".sqlite3", ".fasta")]
    if any(not path.is_file() or path.stat().st_size == 0 for path in outputs):
        raise RuntimeError(f"cblaster database outputs are missing or empty: {attempt}")
    validate_cblaster_database(database.with_suffix(".sqlite3"), paths)
    complete_cache(marker, fingerprint, outputs, config=config, database=str(database.with_suffix(".dmnd")))
    return database.with_suffix(".dmnd")


def run_cblaster(query, database, output_root, command_prefix, versions, unique=3, minimum_hits=3,
                  identity=20, coverage=50, hit_limit=10000, cpus=2, report=print, plot_python_prefix=None):
    import sys

    check_cblaster_query(query, unique, minimum_hits, identity, coverage, hit_limit)
    if cpus < 1:
        raise ValueError("cblaster CPU count must be positive.")
    database = Path(database).resolve()
    config = {"query": file_signature(query), "databases": [file_signature(database), file_signature(database.with_suffix(".sqlite3"))],
              "U": unique, "MH": minimum_hits, "MI": identity, "MC": coverage, "HS": hit_limit,
              "versions": versions, "workflow_version": WORKFLOW_VERSION, "workflow_source_sha256": workflow_source_sha()}
    fingerprint = digest(config)
    folder = Path(output_root).resolve() / fingerprint[:16]
    folder.mkdir(parents=True, exist_ok=True)
    marker = folder / "complete.json"
    if cache_valid(marker, fingerprint):
        report("cblaster: reuse verified search results")
        return read_json(marker)["files"]
    search_marker = folder / "search-complete.json"
    if cache_valid(search_marker, fingerprint):
        files = read_json(search_marker)["files"]
        attempt = Path(files["session"]).parent
        report("cblaster: search is complete; retrying only the plot export")
    else:
        attempt = folder / ("attempt-" + uuid.uuid4().hex[:12])
        attempt.mkdir()
        files = {name: str(attempt / filename) for name, filename in {
            "plot": "plot.html", "session": "session.json", "summary": "summary.tsv", "binary": "abspres.tsv"}.items()}
        args = [*command_prefix, "search", "-qf", str(Path(query).resolve()), "--mode", "local", "-db", str(database),
            "-u", str(unique), "-mh", str(minimum_hits), "-mi", str(identity), "-mc", str(coverage), "-hs", str(hit_limit),
            "--cpus", str(min(cpus, cpu_count())), "--session_file", files["session"],
            "--output", files["summary"], "--output_delimiter", "\t", "--output_decimals", "2",
            "--binary", files["binary"], "--binary_delimiter", "\t"]
        run_logged(args, attempt / "search.log", report=report)
    session = read_json(files["session"])
    if not isinstance(session, dict) or "organisms" not in session:
        raise RuntimeError("cblaster did not produce a valid session JSON; inspect search.log.")
    if not session["organisms"]:
        for key in ("binary", "summary"):
            if not Path(files[key]).exists():
                atomic_write(files[key], "Organism\tScaffold\tStart\tEnd\tScore\n")
    data_outputs = [files[key] for key in ("session", "summary", "binary")]
    if any(not Path(path).is_file() for path in data_outputs):
        raise RuntimeError("cblaster search output is incomplete; see search.log.")
    complete_cache(search_marker, fingerprint, data_outputs, config=config, files=files)
    if not session["organisms"]:
        atomic_write(files["plot"], "<!doctype html><title>cblaster</title><p>No clusters met the search thresholds.</p>")
    else:
        # cblaster 1.4.2's clustering array is not JSON-serializable when >1 result.
        # Convert NumPy values in rendering data; leave saved search results untouched.
        render_script = """import json, sys
from cblaster.classes import Session
from cblaster import plot
with open(sys.argv[1]) as handle:
    session = Session.from_json(handle)
data = plot.get_data(session)
data = json.loads(json.dumps(data, default=lambda value: value.tolist()))
plot.save_html(data, sys.argv[2])
"""
        run_logged([*(plot_python_prefix or [sys.executable]), "-c", render_script, files["session"], files["plot"]],
                   attempt / "plot.log", report=report)
    if any(not Path(path).is_file() for path in files.values()):
        raise RuntimeError("cblaster output is incomplete; see search.log.")
    complete_cache(marker, fingerprint, files.values(), config=config, files=files)
    return files


def normalize_region_base(value):
    value = str(value).lstrip("\ufeff").strip().strip("\"'").split("#", 1)[0].strip()
    value = re.split(r"[/\\]", value)[-1]
    value = re.sub(r"\.(gbk|gbff|gb)$", "", value, flags=re.I)
    return value if re.fullmatch(r"[^\s,;\"']+\.region\d+", value) else None


def extract_region_names(input_path, column=1):
    if column < 1:
        raise ValueError("COL is one-based and must be positive.")
    text = Path(input_path).read_text(encoding="utf-8-sig")
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        raise ValueError("The cblaster table is empty; a header-only table is required for a valid zero-hit result.")
    delimiter = "\t" if "\t" in lines[0] else ","
    rows = list(csv.reader(io.StringIO(text), delimiter=delimiter))
    names, invalid = set(), []
    for number, row in enumerate(rows, 1):
        if not row or not any(cell.strip() for cell in row):
            continue
        if len(row) < column:
            invalid.append(number)
            continue
        value = row[column - 1].strip()
        if number == 1 and value.lower() in {"organism", "scaffold", "name", "filename"}:
            continue
        name = normalize_region_base(value)
        if name:
            names.add(name)
        else:
            invalid.append(number)
    if invalid:
        raise ValueError(f"Region names not recognized in column {column}, rows {invalid[:10]}. Set COL explicitly; no other-column fallback is used.")
    return sorted(names)


def gap_half_open(a_start, a_end, b_start, b_end):
    return max(0, max(a_start, b_start) - min(a_end, b_end))


def diagnose_region(path, gap_max=25000):
    """One implementation for selection and diagnosis; original biological gates retained."""
    from Bio import SeqIO

    if gap_max < 0:
        raise ValueError("PROTO_GAP_MAX must be nonnegative.")
    result = {"path": str(Path(path).resolve()), "pass": False, "error": "", "gap_max": gap_max,
              "has_pks_at": False, "has_ketoacyl_synthase": False, "has_exact_nrps": False,
              "pf00501": 0, "amp_binding": 0, "condensation": 0, "min_core_gap": "",
              "domain_gate": False, "core_gap_gate": False, "failed_rules": ""}
    try:
        text = Path(path).read_text(encoding="utf-8")
        if not text.rstrip().endswith("//"):
            raise ValueError("Incomplete GenBank file: missing record terminator.")
        records = list(SeqIO.parse(io.StringIO(text), "genbank"))
        if not records:
            raise ValueError("No GenBank records found.")
        result["has_pks_at"] = bool(re.search(r"PKS_AT", text, re.I))
        result["has_ketoacyl_synthase"] = bool(re.search(r"ketoacyl synthase", text, re.I))
        domains, clusters = set(), []
        for index, record in enumerate(records):
            for feature in record.features:
                q = feature.qualifiers
                products = q.get("product", [])
                result["has_exact_nrps"] |= "NRPS" in products
                if feature.type == "protocluster":
                    core = " ".join(q.get("core_location", []))
                    spans = [(int(a), int(b)) for a, b in re.findall(r"\[(\d+):(\d+)\]\([+-]\)", core)]
                    if core and (not spans or any(a > b for a, b in spans)):
                        raise ValueError(f"Unsupported protocluster core_location: {core}")
                    for product in products:
                        clusters.append((index, product.strip().upper(), spans))
                if feature.type != "PFAM_domain":
                    continue
                pfams = tuple(sorted({match.group(1) for value in q.get("db_xref", [])
                                      if (match := re.fullmatch(r"(PF\d+)(?:\.\d+)?", value))}))
                domain = (q.get("aSDomain") or [""])[0]
                # Record identity prevents domains on distinct contigs from collapsing.
                key = (index, str(feature.location), pfams, domain,
                       tuple(q.get("protein_start", [])), tuple(q.get("protein_end", [])))
                domains.add(key)
        for _, _, pfams, domain, _, _ in domains:
            result["pf00501"] += "PF00501" in pfams
            result["amp_binding"] += domain.lower() == "amp-binding"
            result["condensation"] += "PF00668" in pfams or domain.lower() in {"condensation_lcl", "condensation_dcl"}
        gaps = [gap_half_open(a, b, c, d)
                for index_a, product_a, spans_a in clusters if product_a == "NRPS"
                for index_b, product_b, spans_b in clusters if "HGLE-KS" in product_b and index_a == index_b
                for a, b in spans_a for c, d in spans_b]
        result["min_core_gap"] = min(gaps) if gaps else ""
        result["domain_gate"] = (result["pf00501"] >= 2 or result["amp_binding"] >= 2) and result["condensation"] >= 2
        result["core_gap_gate"] = bool(gaps) and min(gaps) <= gap_max
        rules = ["has_pks_at", "has_ketoacyl_synthase", "domain_gate", "core_gap_gate", "has_exact_nrps"]
        result["failed_rules"] = ";".join(rule for rule in rules if not result[rule])
        result["pass"] = not result["failed_rules"]
    except Exception as exc:
        result.update(error=str(exc), failed_rules="parse_or_read_error", **{"pass": False})
    return result


DIAGNOSIS_COLUMNS = ["path", "pass", "has_pks_at", "has_ketoacyl_synthase", "pf00501", "amp_binding",
                     "condensation", "domain_gate", "min_core_gap", "gap_max", "core_gap_gate",
                     "has_exact_nrps", "failed_rules", "error"]


def filter_regions(input_table, region_paths, source_root, destination, gap_max=25000,
                   column=1, copy_mode="dir", dry_run=False, report=print):
    if copy_mode not in {"dir", "gbk"} or gap_max < 0:
        raise ValueError("COPY_MODE must be dir or gbk; PROTO_GAP_MAX must be nonnegative.")
    root, dest = Path(source_root).resolve(), Path(destination).resolve()
    if dest == root or root in dest.parents:
        raise ValueError("DEST must be outside the antiSMASH input tree.")
    paths = sorted(set(Path(path).resolve() for path in region_paths))
    if any(root not in path.parents for path in paths):
        raise ValueError("Every region file must be inside ROOT.")
    names = extract_region_names(input_table, column)
    config = {"table": file_signature(input_table), "regions": [file_signature(path) for path in paths],
              "gap_max": gap_max, "column": column, "copy_mode": copy_mode, "dry_run": dry_run,
              "workflow_version": WORKFLOW_VERSION, "workflow_source_sha256": workflow_source_sha()}
    fingerprint = digest(config)
    folder = dest / fingerprint[:16]
    folder.mkdir(parents=True, exist_ok=True)
    marker = folder / "complete.json"
    if cache_valid(marker, fingerprint):
        state = read_json(marker)
        write_json(dest / "latest.json", {"manifest": str(marker), "directory": str(folder)})
        report(f"Filter: reuse verified export ({len(state['kept'])} retained regions)")
        return state
    index = {}
    for path in paths:
        index.setdefault(path.stem, []).append(path)
    matching_errors = [{"name": name, "status": "unmatched" if name not in index else "ambiguous",
                        "paths": ";".join(map(str, index.get(name, [])))}
                       for name in names if len(index.get(name, [])) != 1]
    write_table(folder / "matching_errors.tsv", matching_errors, ["name", "status", "paths"])
    if matching_errors:
        raise ValueError(f"Some region names could not be uniquely matched. Inspect {folder / 'matching_errors.tsv'}.")
    matched = [index[name][0] for name in names]
    diagnosis = [diagnose_region(path, gap_max) for path in matched]
    write_table(folder / "filter_diagnosis.tsv", diagnosis, DIAGNOSIS_COLUMNS)
    write_list(folder / "matched_region_files.abs.txt", matched)
    if any(row["error"] for row in diagnosis):
        raise ValueError(f"GenBank read/parse errors are recorded in {folder / 'filter_diagnosis.tsv'}; export stopped.")
    kept = [Path(row["path"]) for row in diagnosis if row["pass"]]
    targets = sorted({path.parent if copy_mode == "dir" else path for path in kept})
    exported = []
    outputs = [folder / "filter_diagnosis.tsv", folder / "matching_errors.tsv", folder / "matched_region_files.abs.txt"]
    if kept and not dry_run:
        # Each export attempt is isolated. Existing user files and earlier results are preserved.
        export = folder / ("export-" + uuid.uuid4().hex[:12])
        for source in targets:
            target = export / source.relative_to(root)
            target.parent.mkdir(parents=True, exist_ok=True)
            if source.is_dir():
                shutil.copytree(source, target)
                outputs.extend(path for path in target.rglob("*") if path.is_file())
            else:
                shutil.copy2(source, target)
                outputs.append(target)
        exported = [export / path.relative_to(root) for path in kept]
    write_list(folder / "kept_region_files.abs.txt", exported if not dry_run else kept)
    write_list(folder / "kept_source_files.abs.txt", kept)
    outputs += [folder / "kept_region_files.abs.txt", folder / "kept_source_files.abs.txt"]
    state = complete_cache(marker, fingerprint, outputs, config=config, kept=list(map(str, kept)),
                           exported=list(map(str, exported)), diagnosis=str(folder / "filter_diagnosis.tsv"),
                           directory=str(folder), counts={"matched": len(matched), "kept": len(kept)})
    write_json(dest / "latest.json", {"manifest": str(marker), "directory": str(folder)})
    report(f"Filter: {len(kept)}/{len(matched)} retained. Diagnosis and export: {folder}")
    if copy_mode == "dir" and kept:
        report("COPY_MODE=dir includes complete antiSMASH reports and their other regions; kept_region_files.abs.txt identifies only the passing GBKs.")
    return state


def serve_directory(directory, previous_server=None):
    """Notebook preview without changing process cwd or racing to allocate a port."""
    from functools import partial
    from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

    if previous_server is not None:
        previous_server.shutdown()
        previous_server.server_close()

    class QuietHandler(SimpleHTTPRequestHandler):
        def log_message(self, *args):
            pass

    handler = partial(QuietHandler, directory=str(Path(directory).resolve()))
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server
