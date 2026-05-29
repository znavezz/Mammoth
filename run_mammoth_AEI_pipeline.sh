#!/usr/bin/env bash
###############################################################################
# run_mammoth_AEI_pipeline.sh
#
#   Reproducible, reviewable record of the woolly-mammoth per-SINE A-to-I editing
#   index (AEI) pipeline. Re-running this reproduces every step. Each step is one
#   distinct bioinformatic operation, calling a single reusable step-script
#   (alignment and AEI are SEPARATE steps, never merged).
#
#   Organism : woolly mammoth "Reference-assisted 3D" assembly (Sandoval-Velasco
#              et al., Cell 2024). Data + indexes: Mármol-Sánchez Figshare 29590415.
#
#   PROVENANCE / KEY DECISIONS (for review):
#     • Inputs (FASTA, GTF, prebuilt bowtie2 index, Mammoth1.fastq.gz) = Figshare
#       29590415; assembly is NOT in GenBank. The paper publishes NO RMSK and NO BAM.
#     • RMSK generated here (native mammoth coords — the valid choice for AEI;
#       liftover from elephant would corrupt per-base editing signal).
#     • RepeatMasker -species "Elephas maximus" (mammoth's closest living relative;
#       Dfam Mammalia partition 7, installed by setup_dfam_lib.sh).
#     • Aligner = bowtie2 --end-to-end --sensitive (ancient/small-RNA reads); NOT
#       STAR/Salmon (those are for modern RNA-seq). aRNA preprocessing is small-RNA
#       style (cutadapt small-RNA adapter + fastx_collapser + 4+4nt UMI), NOT the
#       lab's standard downloadAndPreprocess.
#     • AEI = Levanon editing_index.nf via runAEI.sh, fed the bowtie2 BAM dir.
#       (Gene quantification is NOT needed for AEI; an expression branch, if wanted,
#       would be added separately using the paper's bowtie2+bedtools-coverage method.)
#
#   USAGE:
#     bash run_mammoth_AEI_pipeline.sh        # self-detaches; Slacks every step
#     tail -f run_mammoth_AEI_pipeline.log
#     _BG_DETACHED=1 bash run_mammoth_AEI_pipeline.sh   # foreground
#
#   Steps run in the FOREGROUND here (via _BG_DETACHED=1) so they execute in strict
#   dependency order; each sub-script's own Slack is suppressed and this script
#   posts one message per step. Needs $NF_WEBHOOK_URL (or $SLACK_WEBHOOK_URL).
###############################################################################
set -euo pipefail

# ── Config (edit here to reproduce for another assembly) ────────────────────
R=/private8/Projects/Nave/reusable
M=/private8/Projects/Nave/Mammoth
G=$M/Genome/Mammoth_SandovalVelasco2024
FAMDB=$R/dfam_famdb
PREFIX=mammoth_SV2024
SPECIES="Elephas maximus"
PA=16
FASTQ=$G/Mammoth1.fastq.gz                 # raw aRNA reads (Figshare file 56345462)
FASTQ_URL="https://ndownloader.figshare.com/files/56345462"
BT2_INDEX=$G/bowtie2_index/Reference-assisted_3D_woolly_mammoth_assembly
PREP_FASTA=$G/preprocess/Mammoth1_Trimmed_Collapsed.fasta
BAMDIR=$G/bowtie2_bam
TX_INDEX=$G/transcriptome_bowtie2_index/$PREFIX     # built in step 7
TX2GENE=$G/${PREFIX}.tx2gene_biotype.tsv            # built in step 7
QUANTDIR=$G/quant                                   # expression counts (step 8)
LOG=$M/run_mammoth_AEI_pipeline.log

# ── Self-detach (survives logout). Foreground: _BG_DETACHED=1 ───────────────
if [[ "${_BG_DETACHED:-}" != "1" ]]; then
  _BG_DETACHED=1 setsid bash "$0" "$@" >"$LOG" 2>&1 </dev/null &
  echo "mammoth AEI pipeline detached (pid $!)."
  echo "  watch: tail -f $LOG"
  echo "  slack: one message per step (needs \$NF_WEBHOOK_URL)"
  exit 0
fi

slack(){ "$R/notify_send.sh" --subject "$1" --text "$2
host: $(hostname) | $(date '+%F %T')" >/dev/null 2>&1 || true; }
banner(){ echo; echo "===================================================================="; \
          echo ">>> $* :: $(date '+%F %T')"; \
          echo "===================================================================="; }
CURRENT="init"
trap 'slack "❌ mammoth AEI pipeline FAILED" "step: $CURRENT — see $LOG"' ERR

slack "▶️ mammoth AEI pipeline START" "prefix=$PREFIX
0+1 download → 2 RepeatMasker → 3 refs → 4 preprocess → 5 align → 6 AEI
→ 7 transcriptome → 8 expression quant"

# ── Step 0+1: downloads (independent → parallel) ────────────────────────────
CURRENT="0+1 downloads (Figshare inputs + reads + Dfam Mammalia partition)"
banner "$CURRENT"
( bash "$G/download_mammoth_SV2024.sh" --bg                       # FASTA + GTF + bowtie2 index
  [[ -s "$FASTQ" ]] || wget -q -O "$FASTQ" "$FASTQ_URL" ) & p0=$! # raw reads (not in the download script)
_BG_DETACHED=1 "$R/setup_dfam_lib.sh" --out_dir "$FAMDB" --partitions 7 --notify false & p1=$!
wait $p0; wait $p1
slack "✅ step 0+1 DONE" "inputs+reads under $G ; Dfam famdb at $FAMDB"

# ── Step 2: RepeatMasker (native mammoth coordinates) ───────────────────────
CURRENT="2 RepeatMasker ($SPECIES, pa=$PA)"
banner "$CURRENT"
_BG_DETACHED=1 "$R/buildRMSK.sh" --fasta "$G/genome.fa" --prefix "$PREFIX" \
  --famdb_dir "$FAMDB" --species "$SPECIES" --pa "$PA" --notify false
slack "✅ step 2 DONE" "$G/RepeatMasker/$PREFIX.repeatMasker.out.gz"

# ── Step 3: SINE-family BEDs + RefSeq BED ───────────────────────────────────
CURRENT="3 SINE-family BEDs + RefSeq BED"
banner "$CURRENT"
_BG_DETACHED=1 "$R/buildSINEbeds.sh" \
  --rmsk "$G/RepeatMasker/$PREFIX.repeatMasker.out.gz" \
  --out_dir "$G/SINE_by_repFamily_bed3_merged" --repeat_class SINE --notify false
_BG_DETACHED=1 "$R/buildAEIrefs.sh" \
  --gtf "$G/annotation.gtf.gz" --out_dir "$G" --prefix "$PREFIX" --notify false
slack "✅ step 3 DONE" "SINE beds + ${PREFIX}.refseq.bed.gz"

# ── Step 4: aRNA preprocessing (cutadapt → collapse → UMI→header) ───────────
CURRENT="4 aRNA preprocess"
banner "$CURRENT"
_BG_DETACHED=1 "$R/preprocess_aRNA.sh" \
  --fastq "$FASTQ" --sample Mammoth1 --out_dir "$G/preprocess" --notify false
slack "✅ step 4 DONE" "$PREP_FASTA"

# ── Step 5: bowtie2 alignment → UMI-dedup BAM ───────────────────────────────
CURRENT="5 bowtie2 align + UMI dedup"
banner "$CURRENT"
_BG_DETACHED=1 "$R/runBowtie2.sh" \
  --reads "$PREP_FASTA" --index "$BT2_INDEX" --sample Mammoth1 \
  --out_dir "$BAMDIR" --notify false
slack "✅ step 5 DONE" "BAM in $BAMDIR"

# ── Step 6: per-SINE AEI (Levanon editing_index.nf via runAEI.sh) ───────────
# NOTE: editing_index.nf is profile-driven (hg38/mm10). For this non-model
# assembly the genome/refseq must be supplied as overrides. VERIFY the exact
# override flag names against the lab's editing_index.nf before the real run.
CURRENT="6 AEI per-SINE"
banner "$CURRENT"
_BG_DETACHED=1 "$R/runAEI.sh" \
  --project_dir "$G" --bams_dir "$BAMDIR" --single_end true \
  --regions "$G/SINE_by_repFamily_bed3_merged/*.bed" --multi_regions true \
  --genome_file "$G/genome.fa" --refseq_file "$G/${PREFIX}.refseq.bed.gz" \
  --notify false
slack "✅ step 6 DONE — AEI complete" "Results under $G/output/.../Editing_index/AEI_perSINE/"

# ── Step 7: transcriptome + bowtie2 index + tx→gene/biotype (expression ref) ─
CURRENT="7 transcriptome annotation"
banner "$CURRENT"
_BG_DETACHED=1 "$R/build_transcriptome_annot.sh" \
  --gtf "$G/annotation.gtf.gz" --genome "$G/genome.fa" --prefix "$PREFIX" \
  --out_dir "$G" --notify false
slack "✅ step 7 DONE" "transcriptome + index + $TX2GENE"

# ── Step 8: aRNA expression quantification (bowtie2 --norc → per-gene counts) ─
CURRENT="8 expression quantification"
banner "$CURRENT"
_BG_DETACHED=1 "$R/runQuantBowtie2.sh" \
  --reads "$PREP_FASTA" --tx_index "$TX_INDEX" --tx2gene "$TX2GENE" \
  --sample Mammoth1 --out_dir "$QUANTDIR" --notify false
slack "✅ step 8 DONE — expression quant complete" "counts in $QUANTDIR"

# ── Summary (audit record) ──────────────────────────────────────────────────
banner "PIPELINE COMPLETE"
cat <<EOF
Outputs:
  $G/genome.fa, annotation.gtf.gz, bowtie2_index/          inputs (Figshare)
  $G/RepeatMasker/$PREFIX.repeatMasker.out.gz              RMSK (this run)
  $G/SINE_by_repFamily_bed3_merged/<family>.bed            per-SINE regions
  $G/${PREFIX}.refseq.bed.gz                               RefSeq BED (AEI)
  $G/preprocess/Mammoth1_Trimmed_Collapsed.fasta           preprocessed reads
  $BAMDIR/Mammoth1_UMIDEDUP.sorted.bam                     alignment
  $G/output/.../Editing_index/AEI_perSINE/<family>/        AEI results
  $G/${PREFIX}.transcriptome.fa, transcriptome_bowtie2_index/   expression ref
  $QUANTDIR/Mammoth1.{transcript,gene,biotype}_counts.tsv  expression counts
  $FAMDB/                                                  Dfam root+Mammalia (reusable)
EOF
slack "🏁 mammoth AEI pipeline FINISHED" "see $LOG"
