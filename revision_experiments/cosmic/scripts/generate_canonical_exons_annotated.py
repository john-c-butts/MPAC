#!/usr/bin/env python3
"""
Generate a corrected, annotated BED of protein_coding, Ensembl_canonical exons.

Source of truth for the target gene set:
  gencode.v44.protein.coding.1kb.promoters.autosomes.v2.full.ID.bed
  (protein_coding, ensembl_canonical, has HGNC name, not read_through, not protein_coding_LoF)

GFF source:
  /projects/tewhey-lab/buttsj/genomes/gencode.v44.basic.annotation.gff3.gz

Output:
  gencode.v44.protein.coding.exons.autosomes.canonical.annotated.bed
  (columns: chrom, start, end, gene_name, score, strand,
            gene_id, transcript_id, exon_id, exon_number)
"""

import gzip
import sys
import re
from collections import defaultdict

PROMOTER_BED = (
    "/pod/2/tewhey-lab/buttsj/Variant_Effects/revision_experiments/cosmic/"
    "raw_data/gencode.v44.protein.coding.1kb.promoters.autosomes.v2.full.ID.bed"
)
GFF_FILE = "/projects/tewhey-lab/buttsj/genomes/gencode.v44.basic.annotation.gff3.gz"
EXON_BED = (
    "/pod/2/tewhey-lab/buttsj/Variant_Effects/revision_experiments/cosmic/"
    "raw_data/gencode.v44.protein.coding.exons.autosomes.v2.bed"
)
OUT_BED = (
    "/pod/2/tewhey-lab/buttsj/Variant_Effects/revision_experiments/cosmic/"
    "raw_data/gencode.v44.protein.coding.exons.autosomes.canonical.annotated.bed"
)

# ── 1. Parse the canonical gene set from the promoter BED ─────────────────────
# Format: chrom  start  end  ENSG_ENST_GENENAME  score  strand
print("Step 1: Parsing canonical gene set from promoter BED...", file=sys.stderr)

canonical_enst = {}   # enst_id (with version) -> gene_name
canonical_ensg = {}   # enst_id -> ensg_id
gene_to_enst   = {}   # gene_name -> enst_id

with open(PROMOTER_BED) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        name_field = parts[3]          # e.g. ENSG00000186092.7_ENST00000641515.2_OR4F5
        tokens = name_field.split("_")
        ensg = tokens[0]
        enst = tokens[1]
        gene = "_".join(tokens[2:])    # some gene names have underscores
        canonical_enst[enst] = gene
        canonical_ensg[enst] = ensg
        gene_to_enst[gene]   = enst

print(f"  {len(canonical_enst):,} canonical transcripts loaded", file=sys.stderr)

# ── 2. Parse the GFF for exons belonging to canonical transcripts ─────────────
print("Step 2: Parsing GFF for canonical exons...", file=sys.stderr)

AUTOSOMES = {f"chr{i}" for i in range(1, 23)}

ATTR_RE = re.compile(r'([^=;]+)=([^;]+)')

def parse_attrs(attr_str):
    return dict(ATTR_RE.findall(attr_str))

canonical_exons = []   # list of dicts

with gzip.open(GFF_FILE, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        chrom, source, feat, start, end, score, strand, phase, attrs = parts
        if feat != "exon":
            continue
        if chrom not in AUTOSOMES:
            continue

        a = parse_attrs(attrs)
        transcript_id = a.get("transcript_id", "")
        if transcript_id not in canonical_enst:
            continue

        # BED format: 0-based start, 1-based end  (GFF is 1-based)
        bed_start = int(start) - 1
        bed_end   = int(end)
        gene_name = canonical_enst[transcript_id]
        gene_id   = canonical_ensg[transcript_id]
        exon_id   = a.get("exon_id", ".")
        exon_num  = a.get("exon_number", ".")

        canonical_exons.append({
            "chrom":         chrom,
            "start":         bed_start,
            "end":           bed_end,
            "gene_name":     gene_name,
            "strand":        strand,
            "gene_id":       gene_id,
            "transcript_id": transcript_id,
            "exon_id":       exon_id,
            "exon_number":   exon_num,
        })

print(f"  {len(canonical_exons):,} canonical exons found", file=sys.stderr)

# ── 3. Sort: chrom (natural), then start ──────────────────────────────────────
def chrom_key(e):
    c = e["chrom"].replace("chr", "")
    return (int(c) if c.isdigit() else 99, e["start"])

canonical_exons.sort(key=chrom_key)

# ── 4. Write annotated BED ────────────────────────────────────────────────────
print(f"Step 3: Writing annotated BED to {OUT_BED}...", file=sys.stderr)
with open(OUT_BED, "w") as out:
    for e in canonical_exons:
        out.write(
            f"{e['chrom']}\t{e['start']}\t{e['end']}\t{e['gene_name']}\t"
            f"0\t{e['strand']}\t{e['gene_id']}\t{e['transcript_id']}\t"
            f"{e['exon_id']}\t{e['exon_number']}\n"
        )
print(f"  Wrote {len(canonical_exons):,} lines", file=sys.stderr)

# ── 5. Diagnostic summary ──────────────────────────────────────────────────────
print("\n========== SUMMARY ==========", file=sys.stderr)
print(f"Canonical gene set (promoter BED): {len(canonical_enst):,} transcripts", file=sys.stderr)

# Count unique transcripts found in GFF
found_enst = {e["transcript_id"] for e in canonical_exons}
missing_enst = set(canonical_enst.keys()) - found_enst
print(f"Canonical transcripts with exons in GFF: {len(found_enst):,}", file=sys.stderr)
if missing_enst:
    print(f"WARNING: {len(missing_enst)} canonical transcripts had NO exons in GFF:", file=sys.stderr)
    for t in sorted(missing_enst)[:20]:
        print(f"  {t}  ({canonical_enst[t]})", file=sys.stderr)

# Compare gene counts
canonical_genes_with_exons = {canonical_enst[t] for t in found_enst}
print(f"Unique genes with canonical exons written: {len(canonical_genes_with_exons):,}", file=sys.stderr)

# ── 6. Compare against the original exon BED ─────────────────────────────────
print("\n--- Comparing against original exon BED ---", file=sys.stderr)

orig_entries = 0
orig_genes   = set()
orig_coords  = set()
with open(EXON_BED) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        orig_entries += 1
        orig_genes.add(p[3])
        orig_coords.add((p[0], int(p[1]), int(p[2]), p[3]))

new_coords = {(e["chrom"], e["start"], e["end"], e["gene_name"]) for e in canonical_exons}

print(f"Original exon BED: {orig_entries:,} entries, {len(orig_genes):,} unique genes", file=sys.stderr)
print(f"New canonical BED: {len(canonical_exons):,} entries, {len(canonical_genes_with_exons):,} unique genes", file=sys.stderr)

reduction = orig_entries - len(canonical_exons)
print(f"Reduction: {reduction:,} entries ({reduction/orig_entries*100:.1f}% removed)", file=sys.stderr)

# Genes in original but not new
genes_only_orig = orig_genes - canonical_genes_with_exons
genes_only_new  = canonical_genes_with_exons - orig_genes

print(f"\nGenes in original BED but NOT in canonical BED: {len(genes_only_orig)}", file=sys.stderr)
if genes_only_orig:
    for g in sorted(genes_only_orig)[:30]:
        print(f"  {g}", file=sys.stderr)

print(f"\nGenes in canonical BED but NOT in original BED: {len(genes_only_new)}", file=sys.stderr)
if genes_only_new:
    for g in sorted(genes_only_new)[:30]:
        print(f"  {g}", file=sys.stderr)

# Coordinate-level comparison (accounting for off-by-1 in original)
# Original may have start+1 error; adjust
orig_coords_adj = {(c, s+1, e, g) for c, s, e, g in orig_coords}  # +1 to original start
in_orig_not_new = orig_coords - new_coords
in_orig_not_new_adj = orig_coords_adj - new_coords
in_new_not_orig = new_coords - orig_coords

print(f"\nCoordinate-level comparison:", file=sys.stderr)
print(f"  Intervals in original NOT in canonical (exact match): {len(in_orig_not_new):,}", file=sys.stderr)
print(f"  Intervals in original NOT in canonical (off-by-1 adjusted): {len(in_orig_not_new_adj):,}", file=sys.stderr)
print(f"  Canonical intervals NOT in original: {len(in_new_not_orig):,}", file=sys.stderr)

print("\nDone.", file=sys.stderr)
