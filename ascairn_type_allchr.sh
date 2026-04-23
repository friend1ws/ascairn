#! /bin/bash

set -e

if [ $# -lt 4 ] || [ $# -gt 5 ]; then
    echo "Usage: $0 <BAM_FILE> <OUTPUT_PREFIX> <DATA_DIR> <REFERENCE: hg38|chm13> [THREAD_NUM]"
    echo ""
    echo "Arguments:"
    echo "  BAM_FILE        Input BAM file"
    echo "  OUTPUT_PREFIX   Prefix for output files"
    echo "  DATA_DIR        Directory containing resource files"
    echo "  REFERENCE       Reference genome build (must be 'hg38' or 'chm13')"
    echo "  THREAD_NUM      Number of threads (optional, default=8)"
    exit 1
fi

BAM_FILE="$1"
OUTPUT_PREFIX="$2"
DATA_DIR="$3"
REFERENCE="$4"
THREAD_NUM="${5:-8}"

export POLARS_MAX_THREADS="${THREAD_NUM}"


if [ "$REFERENCE" != "hg38" ] && [ "$REFERENCE" != "chm13" ]; then
    echo "Error: REFERENCE must be either 'hg38' or 'chm13'."
    exit 1
fi

OUTPUT_DIR="$(dirname "${OUTPUT_PREFIX}")"
if [ ! -d "${OUTPUT_DIR}" ]
then
    mkdir -p "${OUTPUT_DIR}"
fi



##########
ascairn check_depth \
    "${BAM_FILE}" \
    -o "${OUTPUT_PREFIX}.depth.txt" \
    --baseline_region "${DATA_DIR}/chr22_long_arm_${REFERENCE}.bed" \
    --x_region "${DATA_DIR}/chrX_short_arm_${REFERENCE}.bed" \
    -t "${THREAD_NUM}"


ascairn kmer_count \
    "${BAM_FILE}" \
    -o "${OUTPUT_PREFIX}.kmer_count.txt" \
    --kmer_file "${DATA_DIR}/rare_kmer_list.fa" \
    --cen_region "${DATA_DIR}/cen_region_curated_margin_${REFERENCE}.bed" \
    -t "${THREAD_NUM}"


SEX="$(grep Sex ${OUTPUT_PREFIX}.depth.txt | cut -d ' ' -f 2)"

FIRST_CHR=1

for CHR_IND in `seq 1 22` X
do

    commands=(
        "ascairn" "cen_type" \
        "${OUTPUT_PREFIX}.kmer_count.txt" \
        "-o" "${OUTPUT_PREFIX}.chr${CHR_IND}" \
        "--kmer_info" "${DATA_DIR}/kmer_info/chr${CHR_IND}.kmer_info.txt.gz" \
        "--hap_info" "${DATA_DIR}/hap_info/chr${CHR_IND}.hap_info.txt" \
        "--depth_file" "${OUTPUT_PREFIX}.depth.txt"
    )

    if [ $CHR_IND = "X" -a $SEX = "male" ]
    then
        commands+=("--single_hap")
    fi

    ${commands[@]}

    # Aggregate per-chromosome cen_type.txt into a single file
    if [ $FIRST_CHR -eq 1 ]
    then
        echo -ne "Chr\t" > ${OUTPUT_PREFIX}.cen_type_all.txt
        head -n1 ${OUTPUT_PREFIX}.chr${CHR_IND}.cen_type.txt >> ${OUTPUT_PREFIX}.cen_type_all.txt
        FIRST_CHR=0
    fi
    echo -ne "chr${CHR_IND}\t" >> ${OUTPUT_PREFIX}.cen_type_all.txt
    tail -n1 ${OUTPUT_PREFIX}.chr${CHR_IND}.cen_type.txt >> ${OUTPUT_PREFIX}.cen_type_all.txt

done

