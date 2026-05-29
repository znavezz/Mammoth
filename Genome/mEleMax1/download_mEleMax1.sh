#!/bin/bash
# Download mEleMax1 (Asian elephant + mammoth-mito hybrid) reference files,
# matching the alignment context used in Mármol-Sánchez et al., Cell 2025.
#
# Two sources, used together:
#   - Mármol-Sánchez Figshare (article 29590415):
#       Pre-built hybrid genome FASTA = mEleMax1 (Elephas maximus, GCF_024166365.1)
#       nuclear + DQ188829.2 mammoth mitochondrial. Pre-built bowtie2 index +
#       matched GTF. Using these makes our alignments byte-for-byte comparable
#       to the paper.
#   - UCSC GenArk Assembly Hub (GCF_024166365.1):
#       RepeatMasker .out for the mEleMax1 nuclear genome. Figshare does not
#       publish RepeatMasker output for the hybrid; UCSC does for the underlying
#       Elephas maximus assembly. We'll harmonize chromosome naming between
#       Figshare's hybrid FASTA and UCSC's RMSK in the build step.
#
# What this script downloads:
#   Figshare:
#     - AsianElephant_MitoMammoth.fa.gz                   (1.1 GB)  → genome.fa (+ .fai)
#     - AsianElephant_MitoMammoth.gtf.gz                  (20 MB)   → annotation.gtf.gz
#     - AsianElephant_MitoMammoth.gff3.gz                 (24 MB)   → annotation.gff3.gz
#     - AsianElephant_MitoMammoth.{1..4,rev.1,rev.2}.bt2  (~5 GB)   → bowtie2_index/
#     - DQ188829.2_MitoMammoth.bed                        (1.5 KB)  → mito annotation
#     - AsianElephant_MitoMammoth_miRNA_genes.bed         (19 KB)   → miRNA annotation
#   UCSC GenArk:
#     - GCF_024166365.1.repeatMasker.out.gz               (~136 MB) → for AEI SINE BED
#   NCBI:
#     - GCF_024166365.1_assembly_report.txt               (KB)      → chr-name crosswalk
#
# Not downloaded here:
#   - SNP/dbSNP — does not exist for this assembly anywhere; must be derived
#     from matched mammoth WGS (separate script under WGS_for_SNPs/).
#   - Expression BED — built per-run from your own sample counts.
#
# Usage:
#   bash download_mEleMax1.sh
# Progress:
#   tail -f /private8/Projects/Nave/Mammoth/Genome/mEleMax1/download_mEleMax1.log

set -e

GENOME_BASE="/private8/Projects/Nave/Mammoth/Genome/mEleMax1"
LOG_FILE="${GENOME_BASE}/download_mEleMax1.log"
BT2_DIR="${GENOME_BASE}/bowtie2_index"

# ── Sources ───────────────────────────────────────────────────────────────
FIGSHARE_BASE="https://ndownloader.figshare.com/files"
UCSC_BASE="https://hgdownload.soe.ucsc.edu/hubs/GCF/024/166/365/GCF_024166365.1"
NCBI_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/024/166/365/GCF_024166365.1_mEleMax1_primary_haplotype"

# Figshare file IDs (from API listing of article 29590415, 2026)
FIG_FA_GZ_ID=56345834           # AsianElephant_MitoMammoth.fa.gz
FIG_FAI_ID=56345519              # AsianElephant_MitoMammoth.fa.fai
FIG_GTF_ID=56345849              # AsianElephant_MitoMammoth.gtf.gz
FIG_GFF3_ID=56345846             # AsianElephant_MitoMammoth.gff3.gz
FIG_BT2_1_ID=56345597            # .1.bt2
FIG_BT2_2_ID=56345591            # .2.bt2
FIG_BT2_3_ID=56345522            # .3.bt2
FIG_BT2_4_ID=56345588            # .4.bt2
FIG_BT2_R1_ID=56345600           # .rev.1.bt2
FIG_BT2_R2_ID=56345594           # .rev.2.bt2
FIG_MITO_BED_ID=56361398         # DQ188829.2_MitoMammoth.bed
FIG_MIRNA_BED_ID=56362283        # AsianElephant_MitoMammoth_miRNA_genes.bed

# ── Local paths ───────────────────────────────────────────────────────────
GENOME_FA_GZ="${GENOME_BASE}/AsianElephant_MitoMammoth.fa.gz"
GENOME_FAI_SRC="${GENOME_BASE}/AsianElephant_MitoMammoth.fa.fai"
GENOME_FA="${GENOME_BASE}/genome.fa"
GTF_GZ="${GENOME_BASE}/annotation.gtf.gz"
GFF3_GZ="${GENOME_BASE}/annotation.gff3.gz"
MITO_BED="${GENOME_BASE}/DQ188829.2_MitoMammoth.bed"
MIRNA_BED="${GENOME_BASE}/AsianElephant_MitoMammoth_miRNA_genes.bed"
RMSK_OUT_GZ="${GENOME_BASE}/GCF_024166365.1.repeatMasker.out.gz"
ASM_REPORT="${GENOME_BASE}/assembly_report.txt"

THREADS=8

# ── Self-background ───────────────────────────────────────────────────────
if [[ "${1:-}" != "--bg" ]]; then
    mkdir -p "$GENOME_BASE" "$BT2_DIR"
    nohup bash "$0" --bg >> "$LOG_FILE" 2>&1 &
    echo "mEleMax1 download launched in background (PID $!)"
    echo "Follow progress with:"
    echo "  tail -f ${LOG_FILE}"
    exit 0
fi

exec >> "$LOG_FILE" 2>&1
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

log "============================================="
log "  mEleMax1 hybrid download started"
log "  GENOME_BASE   : $GENOME_BASE"
log "  Figshare base : $FIGSHARE_BASE"
log "  UCSC base     : $UCSC_BASE"
log "============================================="

cd "$GENOME_BASE"
mkdir -p "$BT2_DIR"

# ── Helper: idempotent download ───────────────────────────────────────────
fetch() {
    local url="$1" out="$2"
    if [[ -s "$out" ]]; then
        log "  exists ($(du -sh "$out" | cut -f1)): $(basename "$out") — skipping"
    else
        log "  downloading $(basename "$out")"
        log "    from: $url"
        wget --show-progress -O "$out" "$url"
        log "    got:  $(du -sh "$out" | cut -f1)"
    fi
}

# ── Step 1/6: hybrid genome FASTA from Figshare ───────────────────────────
log "[Step 1/6] Hybrid genome FASTA (mEleMax1 nuclear + DQ188829.2 mito)"
fetch "${FIGSHARE_BASE}/${FIG_FA_GZ_ID}" "$GENOME_FA_GZ"
fetch "${FIGSHARE_BASE}/${FIG_FAI_ID}"   "$GENOME_FAI_SRC"

if [[ -s "$GENOME_FA" && -s "${GENOME_FA}.fai" ]]; then
    log "  $GENOME_FA already decompressed and indexed — skipping"
else
    log "  decompressing → genome.fa"
    pigz -dc -p "$THREADS" "$GENOME_FA_GZ" > "$GENOME_FA" 2>/dev/null \
        || gunzip -c "$GENOME_FA_GZ" > "$GENOME_FA"
    log "  samtools faidx genome.fa"
    samtools faidx "$GENOME_FA"
    log "  genome.fa: $(du -sh "$GENOME_FA" | cut -f1), $(wc -l < "${GENOME_FA}.fai") sequences"
fi

# ── Step 2/6: GTF + GFF3 from Figshare ────────────────────────────────────
log "[Step 2/6] Gene annotation (GTF + GFF3)"
fetch "${FIGSHARE_BASE}/${FIG_GTF_ID}"  "$GTF_GZ"
fetch "${FIGSHARE_BASE}/${FIG_GFF3_ID}" "$GFF3_GZ"

# ── Step 3/6: pre-built bowtie2 index from Figshare ───────────────────────
log "[Step 3/6] Pre-built bowtie2 index (~5 GB)"
fetch "${FIGSHARE_BASE}/${FIG_BT2_1_ID}"  "${BT2_DIR}/AsianElephant_MitoMammoth.1.bt2"
fetch "${FIGSHARE_BASE}/${FIG_BT2_2_ID}"  "${BT2_DIR}/AsianElephant_MitoMammoth.2.bt2"
fetch "${FIGSHARE_BASE}/${FIG_BT2_3_ID}"  "${BT2_DIR}/AsianElephant_MitoMammoth.3.bt2"
fetch "${FIGSHARE_BASE}/${FIG_BT2_4_ID}"  "${BT2_DIR}/AsianElephant_MitoMammoth.4.bt2"
fetch "${FIGSHARE_BASE}/${FIG_BT2_R1_ID}" "${BT2_DIR}/AsianElephant_MitoMammoth.rev.1.bt2"
fetch "${FIGSHARE_BASE}/${FIG_BT2_R2_ID}" "${BT2_DIR}/AsianElephant_MitoMammoth.rev.2.bt2"

# ── Step 4/6: small annotation BEDs from Figshare ─────────────────────────
log "[Step 4/6] Mitochondrial + miRNA BEDs"
fetch "${FIGSHARE_BASE}/${FIG_MITO_BED_ID}"  "$MITO_BED"
fetch "${FIGSHARE_BASE}/${FIG_MIRNA_BED_ID}" "$MIRNA_BED"

# ── Step 5/6: RepeatMasker .out from UCSC GenArk ──────────────────────────
log "[Step 5/6] RepeatMasker .out (UCSC GenArk)"
log "  NOTE: UCSC uses chr* naming; Figshare hybrid likely uses NC_* accessions."
log "        Chr-name harmonization happens in build_mEleMax1_AEI_references.sh"
log "        using the assembly_report.txt crosswalk."
fetch "${UCSC_BASE}/GCF_024166365.1.repeatMasker.out.gz" "$RMSK_OUT_GZ"

# ── Step 6/6: NCBI assembly report (chr-name crosswalk) ───────────────────
log "[Step 6/6] Assembly report (from NCBI, for chr-name crosswalk)"
fetch \
    "${NCBI_BASE}/GCF_024166365.1_mEleMax1_primary_haplotype_assembly_report.txt" \
    "$ASM_REPORT"

# ── Sanity check: report chromosome naming on each side ───────────────────
log ""
log "── Chromosome-naming sanity check ──"
log "Figshare hybrid genome.fa first 5 sequences:"
head -5 "${GENOME_FA}.fai" | awk '{print "    " $1 "  (length " $2 ")"}' | tee -a "$LOG_FILE"
log "UCSC RMSK first 5 unique sequences:"
zcat "$RMSK_OUT_GZ" | tail -n +4 | awk '{print $5}' | sort -u | head -5 | sed 's/^/    /' | tee -a "$LOG_FILE"

log "============================================="
log "  ALL DOWNLOADS COMPLETE"
log "  genome.fa            : $GENOME_FA"
log "  annotation.gtf.gz    : $GTF_GZ"
log "  bowtie2 index        : $BT2_DIR/"
log "  RepeatMasker .out.gz : $RMSK_OUT_GZ"
log "  assembly_report.txt  : $ASM_REPORT"
log ""
log "Next step:  bash build_mEleMax1_AEI_references.sh"
log "  (will harmonize chr names between Figshare FASTA and UCSC RMSK,"
log "   then build the 8-col RefSeq BED and per-SINE-family region BEDs)"
log "============================================="
