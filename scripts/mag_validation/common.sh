#!/bin/bash
# common.sh — shared paths/tools/params for the mag_validation pipeline.
# Sourced by every stage script. Override anything via environment before sourcing.

# ---- HPC work dir (this pipeline) ----
# NB: the hpc2021 job environment pre-sets SCRIPTS (and possibly other generic
# names) even under sbatch --export=NONE, so overridable vars use MV_ prefixes
# and the plain names are always ours.
WORK=${MV_WORK:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation}
SCRIPTS=${MV_SCRIPTS:-$WORK/scripts}
LISTS=${MV_LISTS:-$WORK/lists}
SAMPLE_LIST=${MV_SAMPLE_LIST:-$LISTS/sample_list.tsv}   # dataset \t sample \t reads_path

# ---- Input data ----
DATA=${MV_DATA:-/lustre1/g/aos_shihuang/data}
STRAIN_DIR=${STRAIN_DIR:-$DATA/cami2_strain}
MARINE_DIR=${MARINE_DIR:-$DATA/cami2_marine}
GTDB_DIR=${GTDB_DIR:-$DATA/gtdb-r207/genomes_all}
# strain: $STRAIN_DIR/short_read/sample_N/reads/anonymous_reads.fq.gz ; sources: $STRAIN_DIR/source_genomes/short_read/source_genomes/*.fasta
# marine: $MARINE_DIR/marmgCAMI2_sample_N_reads/2018.08.15_09.49.32_sample_N/reads/anonymous_reads.fq ; sources: $MARINE_DIR/genomes/*.fasta

# ---- Tools (absolute paths; nothing is in PATH by default on hpc2021) ----
S2B=${S2B:-/lustre1/g/aos_shihuang/Syn2bANI-hi95/target/release/syn2bani}   # 0.1.0 @ 068119c
SKANI=${SKANI:-/group/aos_shihuang/conda/envs/gtdbtk310/bin/skani}          # 0.3.1
FASTANI=${FASTANI:-/group/aos_shihuang/conda/envs/fastani/bin/fastANI}
ANVIO_BIN=${ANVIO_BIN:-/group/aos_shihuang/conda/envs/anvio/bin}            # dnadiff, nucmer
MEGAHIT=${MEGAHIT:-/group/aos_shihuang/conda/envs/megahit/bin/megahit}      # 1.2.9
MB2_BIN=${MB2_BIN:-/group/aos_shihuang/conda/envs/metabat2/bin}             # metabat2, jgi_summarize_bam_contig_depths
BT2_BIN=${BT2_BIN:-/group/aos_shihuang/conda/envs/bowtie2/bin}              # bowtie2, bowtie2-build, samtools
MINIMAP2=${MINIMAP2:-/group/aos_shihuang/conda/envs/coverm/bin/minimap2}
CHECKM2=${CHECKM2:-/group/aos_shihuang/conda/envs/checkm2/bin/checkm2}
CHECKM2DB=${CHECKM2DB:-/lustre1/g/aos_shihuang/databases/CheckM2_db/uniref100.KO.1.dmnd}  # NOTE: design doc said db/, actual is databases/
PY=${PY:-/group/aos_shihuang/conda/bin/python3}                             # 3.12, has pandas

# ---- Parameters (design doc section 3/4) ----
MIN_BIN_BP=${MIN_BIN_BP:-100000}       # keep bins >= 100 kbp
MB2_MIN_CONTIG=${MB2_MIN_CONTIG:-2500}
ASSIGN_MIN_COV=${ASSIGN_MIN_COV:-0.80} # contig -> source genome: >=80% of contig length
ASSIGN_MIN_ID=${ASSIGN_MIN_ID:-0.95}   #                    ... at >=95% identity
SPECIES_ANI=${SPECIES_ANI:-95.0}       # source-genome species clustering threshold (skani)
AF_STRICT=${AF_STRICT:-60}             # dnadiff AF tier: >=60% strict acceptance
DNADIFF_CHUNK=${DNADIFF_CHUNK:-15}     # pairs per truth array task
ENZYMES=${ENZYMES:-BcgI,AlfI,AloI,FalI}

# reads path helper: reads_path <dataset> <sample>
# marine layout differs by sample: sample_0 has the dated subdir with plain .fq;
# samples 1-9 are extracted flat as reads/anonymous_reads.fq.gz. Probe both.
reads_path() {
    local ds=$1 sm=$2
    if [ "$ds" = "strain" ]; then
        echo "$STRAIN_DIR/short_read/${sm}/reads/anonymous_reads.fq.gz"
    else
        local base="$MARINE_DIR/marmgCAMI2_${sm}_reads"
        if [ -s "$base/2018.08.15_09.49.32_${sm}/reads/anonymous_reads.fq" ]; then
            echo "$base/2018.08.15_09.49.32_${sm}/reads/anonymous_reads.fq"
        else
            echo "$base/reads/anonymous_reads.fq.gz"
        fi
    fi
}
