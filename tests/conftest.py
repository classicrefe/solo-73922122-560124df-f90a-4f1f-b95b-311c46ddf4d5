"""Shared fixtures/helpers for the test suite."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pysam
import pytest

from app.main import Config, create_app

META = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n"
CHROM_LINE = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
HEADER = META + CHROM_LINE


def vcf(records, header=HEADER):
    """Build a VCF body from (pos, ref, alt) tuples or raw strings."""
    lines = []
    for r in records:
        if isinstance(r, str):
            lines.append(r)
        else:
            pos, ref, alt = r
            lines.append(f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.")
    return header + "\n".join(lines) + "\n"


# chr1:100 deletion spanning 100-103, SNP at 105, insertion at 110, del 200-204
SAMPLE = vcf([
    (100, "ACGT", "A"),
    (105, "G", "T"),
    (110, "C", "CT"),
    (200, "TAAAA", "T"),
])


def make_bam(path, reads, chrom="chr1", length=1_000_000, so="coordinate"):
    """Write a coordinate-sorted, indexed BAM from simple read dicts.

    Each read: {"name", "start" (0-based), "seq", optional "flag", "mapq"
    (default 60), "qual" (default all Q40), "cigar" (default all-M)}.
    Returns the BAM path; the .bai sits next to it.
    """
    header = {"HD": {"VN": "1.6", "SO": so},
              "SQ": [{"SN": chrom, "LN": length}]}
    segs = []
    for r in reads:
        a = pysam.AlignedSegment()
        a.query_name = r["name"]
        a.flag = r.get("flag", 0)
        a.reference_id = 0
        a.reference_start = r["start"]
        a.mapping_quality = r.get("mapq", 60)
        a.cigar = r.get("cigar") or [(0, len(r["seq"]))]
        a.query_sequence = r["seq"]
        qual = r.get("qual")
        a.query_qualities = pysam.qualitystring_to_array(
            qual if qual is not None else "I" * len(r["seq"])
        )
        segs.append(a)
    path = str(path)
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for a in sorted(segs, key=lambda x: (x.reference_start, x.query_name)):
            out.write(a)
    pysam.index(path)
    return path


def read_seq(base_at_site, site_offset, length=50, fill="A"):
    """A read sequence whose base at the pileup site is `base_at_site`."""
    seq = [fill] * length
    seq[site_offset] = base_at_site
    return "".join(seq)


def post_read_support(client, bam_path, chrom, start, end):
    with open(bam_path, "rb") as fb, open(bam_path + ".bai", "rb") as fi:
        return client.post(
            "/read-support",
            params={"chrom": chrom, "start": start, "end": end},
            files={"bam": ("reads.bam", fb), "bai": ("reads.bam.bai", fi)},
        )


def make_config(tmp_path, **kw):
    kw.setdefault("instance_id", "test-inst")
    kw.setdefault("reaper_interval", 0)    # deterministic: tests call reap_once()
    kw.setdefault("heartbeat_interval", 0)
    return Config(data_dir=str(tmp_path / "data"), **kw)


@pytest.fixture()
def app(tmp_path):
    return create_app(make_config(tmp_path))


@pytest.fixture()
def client(app):
    from fastapi.testclient import TestClient

    with TestClient(app) as c:
        yield c


def upload(client, body, base_version):
    return client.post(
        "/vcf",
        params={"base_version": base_version},
        files={"file": ("variants.vcf", body if isinstance(body, bytes) else body.encode())},
    )
