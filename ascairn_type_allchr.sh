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
    "${DATA_DIR}/chr22_long_arm_${REFERENCE}.bed" \
    "${OUTPUT_PREFIX}.depth.txt" \
    --x_region_file "${DATA_DIR}/chrX_short_arm_${REFERENCE}.bed" \
    --threads "${THREAD_NUM}"


ascairn kmer_count \
    "${BAM_FILE}" \
    "${DATA_DIR}/rare_kmer_list.fa" \
    "${DATA_DIR}/cen_region_curated_margin_${REFERENCE}.bed" \
    "${OUTPUT_PREFIX}.kmer_count.txt" \
    --threads "${THREAD_NUM}"


DEPTH="$(grep Coverage ${OUTPUT_PREFIX}.depth.txt | cut -d ' ' -f 2)"
SEX="$(grep Sex ${OUTPUT_PREFIX}.depth.txt | cut -d ' ' -f 2)"

echo -e "Chr\tCluster_1\tCluster_2\tHaplotype_1\tHaplotpe_2" > ${OUTPUT_PREFIX}.cen_type.result.txt

for CHR_IND in `seq 1 22` X
do

    commands=(
        "ascairn" "cen_type" \
        "${OUTPUT_PREFIX}.kmer_count.txt" \
        "${OUTPUT_PREFIX}.chr${CHR_IND}" \
        "${DATA_DIR}/kmer_info/chr${CHR_IND}.kmer_info.txt.gz" \
        "${DATA_DIR}/cluster_m3/chr${CHR_IND}.cluster_marker_count.txt.gz" \
        "${DEPTH}" \
        "--cluster_haplotype_file" \
        "${DATA_DIR}/cluster_m3/chr${CHR_IND}.hap_cluster.txt"
    )

    if [ $CHR_IND = "X" -a $SEX = "male" ]
    then
        commands+=("--is_single_hap")
    fi

    ${commands[@]}

    cl1="$(head -n2 ${OUTPUT_PREFIX}.chr${CHR_IND}.cluster.hap_pair.txt | tail -n1 | cut -f 1)"
    hap1="$(head -n2 ${OUTPUT_PREFIX}.chr${CHR_IND}.haplotype.hap_pair.txt | tail -n1 | cut -f 1)"
    if [ $CHR_IND = "X" -a $SEX = "male" ]
    then
        cl2="NA"
        hap2="NA"
    else
        cl2="$(head -n2 ${OUTPUT_PREFIX}.chr${CHR_IND}.cluster.hap_pair.txt | tail -n1 | cut -f 2)"
        hap2="$(head -n2 ${OUTPUT_PREFIX}.chr${CHR_IND}.haplotype.hap_pair.txt | tail -n1 | cut -f 2)"
    fi

    echo -e "chr${CHR_IND}\t${cl1}\t${cl2}\t${hap1}\t${hap2}" >> ${OUTPUT_PREFIX}.cen_type.result.txt

done

