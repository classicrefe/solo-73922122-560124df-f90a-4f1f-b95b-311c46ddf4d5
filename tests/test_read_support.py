"""Read-support endpoint: pileup counting, validation, concurrency."""
import os
import threading

import pytest
from conftest import (
    SAMPLE,
    make_bam,
    make_config,
    post_read_support,
    read_seq,
    upload,
    vcf,
)

from app.main import create_app

SNV = vcf([(105, "G", "T")])
SITE0 = 104  # 0-based position of the SNV at chr1:105


def _read(name, base, start=80, **kw):
    d = {"name": name, "start": start, "seq": read_seq(base, SITE0 - start)}
    d.update(kw)
    return d


# BQ strings with a low-quality base exactly at the site offset.
LOWBQ_24 = "I" * 24 + "+" + "I" * 25   # Q10 at offset 24 (read starting at 80)
LOWBQ_4 = "I" * 4 + "+" + "I" * 45     # Q10 at offset 4 (read starting at 100)
BQ20 = "5" * 50                        # Q20 everywhere


def site_reads():
    """Reads at chr1:105 exercising every counting rule."""
    return [
        _read("ref1", "G"),
        _read("ref2", "G"),
        _read("alt1", "T"),
        _read("oth1", "C"),
        _read("nbase", "N"),                       # N supports nothing
        # boundary qualities: exactly 20 qualifies
        _read("mq20", "G", mapq=20),
        _read("bq20", "G", qual=BQ20),
        # pair, both ends qualified, agree on T -> one ALT molecule
        _read("pairOK", "T", flag=99),
        _read("pairOK", "T", start=100, flag=147),
        # pair, both ends qualified, disagree -> conflict, no support
        _read("pairX", "T", flag=99),
        _read("pairX", "A", start=100, flag=147),
        # pair, mate below base-quality threshold -> qualified end only
        _read("pairQ", "G", flag=99),
        _read("pairQ", "A", start=100, flag=147, qual=LOWBQ_4),
        # pair, mate deleted at the site -> based end only
        _read("pairD", "T", flag=99),
        _read("pairD", "A", start=100, flag=147, cigar=[(0, 4), (2, 1), (0, 46)]),
        # excluded categories, all carrying ALT
        _read("dup", "T", flag=0x400),
        _read("sec", "T", flag=0x100),
        _read("supp", "T", flag=0x800),
        _read("lowmq", "T", mapq=19),
        _read("lowbq", "T", qual=LOWBQ_24),
        _read("delr", "A", cigar=[(0, 24), (2, 1), (0, 26)]),
        _read("skipr", "A", cigar=[(0, 24), (3, 1), (0, 26)]),
    ]


def test_snv_support_counts_pairs_and_filters(client, tmp_path):
    assert upload(client, SNV, base_version=0).status_code == 201
    bam = make_bam(tmp_path / "site.bam", site_reads())
    r = post_read_support(client, bam, "chr1", 100, 110)
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["version"] == 1
    assert body["chrom"] == "chr1"
    assert body["n_records"] == 1
    res = body["results"][0]
    assert res["record"] == "chr1\t105\t.\tG\tT\t.\t.\t."
    assert res["applicable"] is True
    # ref: ref1, ref2, mq20, bq20, pairQ (mate low-BQ) = 5
    # alt: alt1, pairOK (both ends, one molecule), pairD (mate deleted) = 3
    # other: oth1 (C) = 1; pairX conflicts and counts for nothing
    assert res["support"] == {"ref": 5, "alt": 3, "other": 1}
    assert res["other_bases"] == {"C": 1}
    assert res["conflicts"] == [{"qname": "pairX", "bases": ["A", "T"]}]


def test_non_snv_records_marked_not_applicable(client, tmp_path):
    # SAMPLE: del@100 (ACGT->A), snp@105 (G->T), ins@110 (C->CT)
    assert upload(client, SAMPLE, base_version=0).status_code == 201
    bam = make_bam(tmp_path / "m.bam", [_read("alt1", "T")])
    r = post_read_support(client, bam, "chr1", 100, 110)
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["n_records"] == 3
    by_pos = {e["pos"]: e for e in body["results"]}
    assert by_pos[100]["applicable"] is False
    assert by_pos[100]["reason"] == "not a single-base substitution"
    assert by_pos[110]["applicable"] is False
    assert by_pos[105]["applicable"] is True
    assert by_pos[105]["support"] == {"ref": 0, "alt": 1, "other": 0}
    # every entry carries the original published record line
    assert by_pos[100]["record"] == "chr1\t100\t.\tACGT\tA\t.\t.\t."
    assert by_pos[110]["record"] == "chr1\t110\t.\tC\tCT\t.\t.\t."


def test_zero_coverage_site(client, tmp_path):
    upload(client, SNV, base_version=0)
    bam = make_bam(tmp_path / "z.bam",
                   [{"name": "far", "start": 5000, "seq": "A" * 50}])
    r = post_read_support(client, bam, "chr1", 100, 110)
    assert r.status_code == 200
    res = r.json()["results"][0]
    assert res["support"] == {"ref": 0, "alt": 0, "other": 0}
    assert res["other_bases"] == {}
    assert res["conflicts"] == []


def test_interval_without_records(client, tmp_path):
    upload(client, SNV, base_version=0)
    bam = make_bam(tmp_path / "e.bam", [_read("ref1", "G")])
    r = post_read_support(client, bam, "chr1", 500, 600)
    assert r.status_code == 200
    assert r.json()["n_records"] == 0
    assert r.json()["results"] == []


def test_response_carries_pinned_version_and_record(client, tmp_path):
    upload(client, SNV, base_version=0)
    upload(client, vcf([(105, "G", "T"), (300, "A", "C")]), base_version=1)
    bam = make_bam(tmp_path / "v.bam", [_read("ref1", "G")])
    r = post_read_support(client, bam, "chr1", 100, 110)
    body = r.json()
    assert body["version"] == 2
    assert body["results"][0]["record"] == "chr1\t105\t.\tG\tT\t.\t.\t."


# ------------------------------------------------------------------ validation


def test_interval_span_limit(client):
    upload(client, SNV, base_version=0)
    r = client.post(
        "/read-support",
        params={"chrom": "chr1", "start": 1, "end": 2001},
        files={"bam": ("r.bam", b"x"), "bai": ("r.bam.bai", b"y")},
    )
    assert r.status_code == 400
    assert "2001" in r.json()["detail"]


def test_interval_bounds_validated(client):
    upload(client, SNV, base_version=0)
    files = {"bam": ("r.bam", b"x"), "bai": ("r.bam.bai", b"y")}
    r = client.post("/read-support",
                    params={"chrom": "chr1", "start": 10, "end": 5}, files=files)
    assert r.status_code == 400
    r = client.post("/read-support",
                    params={"chrom": "chr1", "start": 0, "end": 5}, files=files)
    assert r.status_code == 422


def test_more_than_20_records_rejected(client, tmp_path):
    body = vcf([(100 + i, "G", "T") for i in range(21)])
    assert upload(client, body, base_version=0).status_code == 201
    bam = make_bam(tmp_path / "x.bam", [_read("ref1", "G")])
    r = post_read_support(client, bam, "chr1", 1, 2000)
    assert r.status_code == 400
    assert "21" in r.json()["detail"]


def test_exactly_20_records_allowed(client, tmp_path):
    body = vcf([(100 + i, "G", "T") for i in range(20)])
    assert upload(client, body, base_version=0).status_code == 201
    bam = make_bam(tmp_path / "x.bam", [_read("ref1", "G")])
    r = post_read_support(client, bam, "chr1", 1, 2000)
    assert r.status_code == 200, r.text
    assert r.json()["n_records"] == 20


def test_request_chromosome_mismatch(client, tmp_path):
    upload(client, SNV, base_version=0)
    bam = make_bam(tmp_path / "c.bam", [_read("ref1", "G")])
    r = post_read_support(client, bam, "chr2", 100, 110)
    assert r.status_code == 400
    assert "chr1" in r.json()["detail"]


def test_bam_missing_chromosome(client, tmp_path):
    upload(client, SNV, base_version=0)
    bam = make_bam(tmp_path / "c2.bam", [_read("ref1", "G")], chrom="chr2")
    r = post_read_support(client, bam, "chr1", 100, 110)
    assert r.status_code == 400
    assert "chr1" in r.json()["detail"]


def test_unreadable_bam_rejected(client):
    upload(client, SNV, base_version=0)
    r = client.post(
        "/read-support",
        params={"chrom": "chr1", "start": 100, "end": 110},
        files={"bam": ("r.bam", b"this is not a bam"),
               "bai": ("r.bam.bai", b"nor an index")},
    )
    assert r.status_code == 400
    assert "BAM" in r.json()["detail"]


def test_unreadable_bai_rejected(client, tmp_path):
    upload(client, SNV, base_version=0)
    bam = make_bam(tmp_path / "g.bam", [_read("ref1", "G")])
    with open(bam, "rb") as fh:
        bam_bytes = fh.read()
    r = client.post(
        "/read-support",
        params={"chrom": "chr1", "start": 100, "end": 110},
        files={"bam": ("r.bam", bam_bytes), "bai": ("r.bam.bai", b"garbage")},
    )
    assert r.status_code == 400
    assert "BAI" in r.json()["detail"]


def test_bad_bam_rejected_even_without_records(client):
    """BAM/BAI validity is enforced even when no record needs a pileup."""
    upload(client, SNV, base_version=0)
    r = client.post(
        "/read-support",
        params={"chrom": "chr1", "start": 500, "end": 600},
        files={"bam": ("r.bam", b"junk"), "bai": ("r.bam.bai", b"junk")},
    )
    assert r.status_code == 400


def test_unsorted_bam_rejected(client, tmp_path):
    upload(client, SNV, base_version=0)
    bam = make_bam(tmp_path / "u.bam", [_read("ref1", "G")], so="queryname")
    r = post_read_support(client, bam, "chr1", 100, 110)
    assert r.status_code == 400
    assert "sorted" in r.json()["detail"]


def test_read_support_without_published_vcf_is_404(client, tmp_path):
    bam = make_bam(tmp_path / "n.bam", [_read("ref1", "G")])
    r = post_read_support(client, bam, "chr1", 1, 10)
    assert r.status_code == 404


def test_reads_upload_size_limit(tmp_path):
    from fastapi.testclient import TestClient

    bam = make_bam(tmp_path / "s.bam", [_read("ref1", "G")])
    total = os.path.getsize(bam) + os.path.getsize(bam + ".bai")

    # one byte under the combined size -> 413
    app1 = create_app(make_config(tmp_path, instance_id="lim-1",
                                  max_reads_upload_bytes=total - 1))
    with TestClient(app1) as c:
        assert upload(c, SNV, base_version=0).status_code == 201
        r = post_read_support(c, bam, "chr1", 100, 110)
        assert r.status_code == 413

    # exactly at the combined size -> accepted
    app2 = create_app(make_config(tmp_path / "second", instance_id="lim-2",
                                  max_reads_upload_bytes=total))
    with TestClient(app2) as c:
        assert upload(c, SNV, base_version=0).status_code == 201
        r = post_read_support(c, bam, "chr1", 100, 110)
        assert r.status_code == 200, r.text


# ------------------------------------------------------------------ concurrency


def test_slow_read_support_blocks_neither_upload_nor_query(tmp_path, monkeypatch):
    """A running statistics job pins its version but must not serialize
    corrections (uploads) or the original interval queries behind it."""
    from fastapi.testclient import TestClient

    from app import main as main_mod
    from app import readsupport

    app1 = create_app(make_config(tmp_path, instance_id="inst-1"))
    app2 = create_app(make_config(tmp_path, instance_id="inst-2"))

    started = threading.Event()
    release = threading.Event()
    real_collect = readsupport.collect_support

    def slow_collect(bam_path, chrom, records):
        started.set()
        release.wait(timeout=60)
        return real_collect(bam_path, chrom, records)

    monkeypatch.setattr(main_mod, "collect_support", slow_collect)

    with TestClient(app1) as c1, TestClient(app2) as c2:
        assert upload(c1, SNV, base_version=0).status_code == 201
        bam = make_bam(tmp_path / "conc.bam", [_read("alt1", "T")])
        holder = {}

        def call():
            holder["resp"] = post_read_support(c1, bam, "chr1", 100, 110)

        t = threading.Thread(target=call)
        t.start()
        try:
            assert started.wait(timeout=30), "read-support never reached the pileup"
            # correcting the variant list proceeds while statistics run
            r = upload(c2, vcf([(500, "A", "G")]), base_version=1)
            assert r.status_code == 201, r.text
            assert r.json()["version"] == 2
            # the original interval query is served, on both instances
            q2 = c2.get("/variants", params={"chrom": "chr1", "start": 1, "end": 1000})
            assert q2.status_code == 200 and q2.json()["version"] == 2
            q1 = c1.get("/variants", params={"chrom": "chr1", "start": 1, "end": 1000})
            assert q1.status_code == 200 and q1.json()["version"] == 2
        finally:
            release.set()
        t.join(timeout=30)
        assert not t.is_alive(), "read-support request did not finish"

    resp = holder["resp"]
    assert resp.status_code == 200, resp.text
    body = resp.json()
    # the statistics still describe the version pinned when the request began
    assert body["version"] == 1
    assert body["results"][0]["record"] == "chr1\t105\t.\tG\tT\t.\t.\t."
    assert body["results"][0]["support"] == {"ref": 0, "alt": 1, "other": 0}


def test_openapi_documents_read_support(client):
    assert client.get("/docs").status_code == 200
    spec = client.get("/openapi.json").json()
    assert "/read-support" in spec["paths"]
    post = spec["paths"]["/read-support"]["post"]
    params = {p["name"] for p in post["parameters"]}
    assert {"chrom", "start", "end"} <= params
