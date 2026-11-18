"""Read-support statistics: pileup-based evidence for published SNV records.

For each single-base-substitution record intersecting the requested interval,
count how many independent molecules in the uploaded BAM support the REF
allele, the ALT allele, or any other base.

Counting rules (enforced explicitly in Python, not left to htslib defaults):

- a read end qualifies only when its mapping quality is >= 20 and the base
  quality at the site is >= 20;
- duplicate, secondary and supplementary alignments never count;
- positions where the read carries a deletion or a reference skip (CIGAR
  D/N at the site) never count;
- when both ends of a read pair cover the site, the pair is one molecule:
  agreement counts once; qualified disagreement is reported as a conflict
  and counts for nothing.
"""
from __future__ import annotations

import pysam

MIN_MAPPING_QUALITY = 20
MIN_BASE_QUALITY = 20
MAX_INTERVAL_SPAN = 2000          # closed-interval length limit, in bases
MAX_RECORDS = 20                  # max VCF records intersecting the interval

# BAM flag bits that disqualify a read: UNMAP | SECONDARY | SUPPLEMENTARY | DUP
EXCLUDED_FLAG_MASK = 0x004 | 0x100 | 0x800 | 0x400

_VALID_BASES = frozenset("ACGT")


class BamError(ValueError):
    """The uploaded BAM/BAI cannot be read or does not match the chromosome."""


def _qualified_base(pileup_read) -> str | None:
    """Base observed at the pileup column, or None when it must not count."""
    aln = pileup_read.alignment
    if aln.flag & EXCLUDED_FLAG_MASK:
        return None
    if aln.mapping_quality < MIN_MAPPING_QUALITY:
        return None
    if pileup_read.is_del or pileup_read.is_refskip:
        return None
    qp = pileup_read.query_position
    if qp is None:
        return None
    quals = aln.query_qualities
    if quals is None or quals[qp] < MIN_BASE_QUALITY:
        return None
    base = aln.query_sequence[qp].upper()
    # An N (or any non-ACGT call) supports no allele; skip it.
    return base if base in _VALID_BASES else None


def _snv_support(bam, chrom: str, pos1: int, ref: str, alt: str) -> dict:
    """Count ref/alt/other molecules at a 1-based single-base position."""
    pairs: dict[str, list[str]] = {}   # query name -> qualified base calls
    for col in bam.pileup(
        chrom,
        pos1 - 1,
        pos1,
        truncate=True,
        max_depth=1_000_000,
        stepper="samtools",
        # All filtering is done in _qualified_base; keep htslib permissive
        # (and disable BAQ so raw base qualities are compared).
        min_mapping_quality=0,
        min_base_quality=0,
        compute_baq=False,
        ignore_overlaps=False,
        ignore_orphans=False,
        flag_filter=0,
    ):
        if col.reference_pos != pos1 - 1:
            continue
        for pr in col.pileups:
            base = _qualified_base(pr)
            if base is not None:
                pairs.setdefault(pr.alignment.query_name, []).append(base)

    n_ref = n_alt = 0
    other: dict[str, int] = {}
    conflicts = []
    for qname, bases in pairs.items():
        if len(bases) > 1:
            # Both ends of one molecule cover the site.
            if any(b != bases[0] for b in bases):
                conflicts.append({"qname": qname, "bases": sorted(bases)})
                continue
            bases = [bases[0]]          # ends agree: count the molecule once
        base = bases[0]
        if base == ref:
            n_ref += 1
        elif base == alt:
            n_alt += 1
        else:
            other[base] = other.get(base, 0) + 1
    conflicts.sort(key=lambda c: c["qname"])
    return {
        "support": {"ref": n_ref, "alt": n_alt, "other": sum(other.values())},
        "other_bases": dict(sorted(other.items())),
        "conflicts": conflicts,
    }


def collect_support(bam_path: str, chrom: str, records: list) -> list:
    """Compute per-record read support for the published records.

    `records` is a list of (pos, ref, alt, raw_line) tuples. Non-SNV records
    are marked not applicable. Raises BamError when the BAM/BAI cannot be
    read, is not coordinate-sorted, or lacks the chromosome.
    """
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except (OSError, ValueError) as exc:
        raise BamError(f"BAM cannot be read: {exc}") from exc
    try:
        so = bam.header.get("HD", {}).get("SO")
        if so is not None and so.lower() != "coordinate":
            raise BamError(f"BAM is not coordinate-sorted (header SO={so!r})")
        if chrom not in bam.references:
            refs = ", ".join(bam.references) or "none"
            raise BamError(f"chromosome '{chrom}' not in BAM references ({refs})")
        # Force the index to load so a missing/unreadable BAI fails the whole
        # request here, even when no SNV record needs a pileup.
        try:
            for _ in bam.fetch(chrom, 0, 1):
                break
        except (OSError, ValueError) as exc:
            raise BamError(f"BAI index cannot be read: {exc}") from exc

        results = []
        for pos, ref, alt, line in records:
            entry = {"record": line, "pos": pos, "ref": ref, "alt": alt}
            if len(ref) == 1 and len(alt) == 1:
                entry["applicable"] = True
                entry.update(_snv_support(bam, chrom, pos, ref, alt))
            else:
                entry["applicable"] = False
                entry["reason"] = "not a single-base substitution"
            results.append(entry)
        return results
    except BamError:
        raise
    except (OSError, ValueError) as exc:
        raise BamError(f"BAM cannot be read: {exc}") from exc
    finally:
        bam.close()
