# Versioned VCF Variant Store

FastAPI + pysam + SQLite service for uploading single-chromosome VCF files and
querying them by 1-based closed chromosome intervals, with atomic versioned
replacement and safe concurrent access from two app instances. A
`/read-support` endpoint cross-checks published single-base substitutions
against an uploaded BAM to flag suspected false positives. Interactive API
docs (Swagger UI) are served at `/docs` on every instance.

## Run

```bash
docker compose up --build
```

Two instances come up sharing one data volume:

- instance 1: http://localhost:8001/docs
- instance 2: http://localhost:8002/docs

Both serve the same published data; uploads may go to either one.

## API

### `POST /vcf?base_version=<n>` (multipart file field `file`)

Accepts a VCF that is:

- at most 10 MiB,
- single-chromosome (all records share one CHROM),
- exactly 8 columns (`#CHROM POS ID REF ALT QUAL FILTER INFO`, no FORMAT/sample
  columns),
- biallelic with `REF`/`ALT` containing only `A/C/G/T` (uppercase; no `N`, no
  symbolic alleles, no `*`, no multi-allelic `A,G`).

Records are sorted by position, BGZF-compressed and Tabix-indexed. The upload
carries the version the client last saw (`base_version`, `0` when the store is
empty). If two replacements race — even across the two instances — exactly one
wins; the loser gets `409` and must re-read `/version` and retry.
Returns `201` with `{"version", "chrom", "n_records"}`.

### `GET /variants?chrom=<c>&start=<s>&end=<e>`

1-based closed interval `[start, end]`. Returns every raw VCF record whose
reference span `[POS, POS+len(REF)-1]` intersects the interval — including
deletions that start before the interval. Fetching goes through the Tabix
index (`TabixFile.fetch`), never a full scan. The response carries the
`version` that was actually read; that version's data and index are pinned for
the duration of the query, so version, data and index always match.

### `GET /version`

Current published version, chromosome and record count.

### `POST /read-support?chrom=<c>&start=<s>&end=<e>` (multipart file fields `bam`, `bai`)

Cross-checks the published records against raw aligned reads, for QC of
suspected false-positive calls.

- `chrom` must be the published chromosome; `start`/`end` are a 1-based
  **closed** interval spanning at most **2000 bases**.
- `bam` is a coordinate-sorted BAM, `bai` its index; **at most 20 MiB
  combined** (`413` beyond that).
- The interval must intersect at most **20 published records** (`400`
  otherwise); records are fetched through the same Tabix index as
  `GET /variants`.
- An unreadable BAM/BAI, a non-coordinate-sorted BAM, or a BAM that does not
  contain the chromosome fails the whole request with `400`.

For every intersected record the response carries the original record line
and the pinned `version` it was read from. Single-base substitutions get
molecule counts:

```json
{
  "version": 1, "chrom": "chr1", "start": 100, "end": 110, "n_records": 2,
  "criteria": {"min_mapping_quality": 20, "min_base_quality": 20, "...": "..."},
  "results": [
    {"record": "chr1\t105\t.\tG\tT\t.\t.\t.", "pos": 105, "ref": "G", "alt": "T",
     "applicable": true,
     "support": {"ref": 5, "alt": 3, "other": 1},
     "other_bases": {"C": 1},
     "conflicts": [{"qname": "pairX", "bases": ["A", "T"]}]},
    {"record": "chr1\t100\t.\tACGT\tA\t.\t.\t.", "pos": 100, "ref": "ACGT",
     "alt": "A", "applicable": false, "reason": "not a single-base substitution"}
  ]
}
```

Counting rules: only read ends with mapping quality ≥ 20 and base quality
≥ 20 at the site count; duplicate, secondary and supplementary alignments
and positions deleted/skipped in the read (CIGAR `D`/`N`) are excluded;
non-`ACGT` calls support nothing. Both ends of a read pair covering the
site are one molecule: agreement counts once, qualified disagreement is
listed under `conflicts` and counts for nothing. Non-SNV records are
marked `applicable: false` individually.

The pileup runs in a worker thread while holding only a read pin on the
version it describes, so a slow `/read-support` job never blocks uploads
(correcting the variant list) or `GET /variants` queries; the response
still reports the fixed version the statistics belong to.

## Consistency design

- **Atomic publish**: a new version is fully built (sorted, compressed,
  indexed) in a per-instance staging directory, moved into `versions/`, and
  only then the current-version pointer in SQLite is flipped inside a single
  `BEGIN IMMEDIATE` transaction. Readers never see a half-built file set.
- **Optimistic concurrency**: the pointer flip checks `base_version` inside
  the same transaction, so concurrent replacements from both instances are
  serialized and only one succeeds.
- **Old queries finish**: a query pins its version in a `refs` table before
  opening files. A background reaper deletes a superseded version's data and
  index only when no pin references it. Instances heartbeat; pins of crashed
  instances expire and never block reclamation.
- **Crash recovery**: if the process dies mid-build or mid-publish, the
  staging/orphan directories left behind are swept (with a grace period) and
  the previously published version keeps serving. The pointer flip is the last
  step, so a restart never exposes a partial file set.

## Local development

```bash
pip install -r requirements.txt
python -m pytest tests/          # 56 tests
uvicorn --factory app.main:create_app   # single dev instance (uses /data or $VCF_DATA_DIR)
```

Configuration via environment: `VCF_DATA_DIR`, `VCF_INSTANCE_ID`,
`VCF_MAX_UPLOAD_BYTES`, `VCF_MAX_READS_UPLOAD_BYTES`, `VCF_REAPER_INTERVAL`,
`VCF_HEARTBEAT_INTERVAL`, `VCF_DEAD_AFTER`, `VCF_ORPHAN_GRACE`.
