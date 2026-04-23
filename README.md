<div align="center">
  <img src="image/karamatsu_cairn.png" alt="Cairn at Mt. Karamatsu" width="250">
  <p>Maruyama cairn at Mt. Karamatsu</p>
</div>

# ascairn

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![CI](https://github.com/friend1ws/ascairn/actions/workflows/python-test.yml/badge.svg)

`ascairn` (alpha-satellite cairn) is software for estimating centromere variation from short-read sequencing data using rare k-mers within centromere sequences.
For each chromosome, ascairn identifies the most likely pair of centromere clusters and the nearest haplotypes from a reference panel ([Shiraishi et al., bioRxiv, 2025](https://doi.org/10.1101/2025.07.26.666712)).

Both GRCh38 (hg38) and T2T-CHM13 (chm13) reference genomes are supported. Centromere region BED files for each reference are included in the resource repository.

<div align="center">
  <img src="image/categorization_overview.png" alt="Overview of centromere classification using rare k-mer" width="750">
</div>

## Prerequisites

### Software
- [samtools](https://github.com/samtools/samtools)
- [Jellyfish](https://github.com/gmarcais/Jellyfish)
- [mosdepth](https://github.com/brentp/mosdepth)

### Python packages
- click
- scipy
- polars
- boto3 (required only for accessing CRAM files in Amazon S3)

## Installation

1. Install prerequisite software and ensure they are accessible via your `PATH`.

2. Install `ascairn`.
```bash
git clone https://github.com/friend1ws/ascairn.git
cd ascairn
pip install .   # or pip install -e . for development
```

3. Download [ascairn resource files](https://github.com/friend1ws/ascairn_resource).
This repository contains the reference data required by ascairn, including:
   - Rare k-mer list for alpha satellite sequences (`rare_kmer_list.fa`)
   - Centromere region BED files for hg38 and chm13
   - Per-chromosome k-mer information (`kmer_info/`)
   - Per-chromosome haplotype-to-cluster mapping (`hap_info/`)

```bash
git clone https://github.com/friend1ws/ascairn_resource.git
```

After installation, your directory structure should look like this:
```
ascairn/
ascairn_resource/
└── resource/
    ├── common/
    │   ├── cen_region_curated_margin_hg38.bed
    │   ├── cen_region_curated_margin_chm13.bed
    │   ├── chr22_long_arm_hg38.bed
    │   ├── chr22_long_arm_chm13.bed
    │   ├── chrX_short_arm_hg38.bed
    │   └── chrX_short_arm_chm13.bed
    └── panel/
        └── ascairn_paper_2025/
            ├── rare_kmer_list.fa
            ├── kmer_info/
            │   ├── chr1.kmer_info.txt.gz
            │   ├── ...
            │   └── chrX.kmer_info.txt.gz
            └── hap_info/
                ├── chr1.hap_info.txt
                ├── ...
                └── chrX.hap_info.txt
```

## Quick Start

### 1. Prepare the sequence data

We use NA12877 from the 1000 Genomes Project, whose GRCh38-aligned CRAM is publicly available on AWS S3 and FTP.

**Option 1: Direct S3 path (no download required)**

If samtools is properly installed with S3 support, ascairn can read CRAM files directly from S3 without downloading. The file is:
```
s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/ERR3989340/NA12877.final.cram
```

**Option 2: Download locally**

Either via AWS CLI (public access, no AWS credentials required):
```
aws s3 cp --no-sign-request s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/ERR3989340/NA12877.final.cram seq_data/
aws s3 cp --no-sign-request s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/ERR3989340/NA12877.final.cram.crai seq_data/
```

Or via FTP:
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989340/NA12877.final.cram -P seq_data/
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989340/NA12877.final.cram.crai -P seq_data/
```

### 2. Run the workflow

The simplest way to run ascairn is via `type_all`, which executes all steps (`check_depth`, `kmer_count`, `cen_type`) for chr1-22, chrX, and chrY (males only) in a single command.

**Using Option 1 (direct S3 path):**
```bash
ascairn type_all \
    s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/ERR3989340/NA12877.final.cram \
    -o output/NA12877 \
    --resource_dir ascairn_resource/resource/panel/ascairn_paper_2025 \
    --reference hg38 \
    -t 8
```

**Using Option 2 (downloaded CRAM):**
```bash
ascairn type_all \
    seq_data/NA12877.final.cram \
    -o output/NA12877 \
    --resource_dir ascairn_resource/resource/panel/ascairn_paper_2025 \
    --reference hg38 \
    -t 8
```

### 3. Check the results

After successful execution, the main output file is:
```
output/NA12877.cen_type_all.txt
```

See [Output Format](#output-format) for how to interpret the results.

### Running with CHM13-aligned data

For users with CHM13-aligned data, here is an equivalent example using NA12878, whose CHM13-aligned CRAM is publicly available at the DDBJ mirror.

**Option 1: Direct HTTPS URL (no download required)**

If samtools is properly installed with libcurl support, ascairn can read CRAM files directly from HTTPS URLs:

```bash
ascairn type_all \
    https://ddbj.nig.ac.jp/public/public-human-genomes/CHM13/1000Genomes/CRAM/NA12878/NA12878.cram \
    -o output/NA12878 \
    --resource_dir ascairn_resource/resource/panel/ascairn_paper_2025 \
    --reference chm13 \
    -t 8
```

**Option 2: Download locally**

```
wget https://ddbj.nig.ac.jp/public/public-human-genomes/CHM13/1000Genomes/CRAM/NA12878/NA12878.cram -P seq_data/
wget https://ddbj.nig.ac.jp/public/public-human-genomes/CHM13/1000Genomes/CRAM/NA12878/NA12878.cram.crai -P seq_data/

ascairn type_all \
    seq_data/NA12878.cram \
    -o output/NA12878 \
    --resource_dir ascairn_resource/resource/panel/ascairn_paper_2025 \
    --reference chm13 \
    -t 8
```

The only differences from the GRCh38 workflow are (i) the input CRAM aligned to CHM13 and (ii) `--reference chm13` in the command.

## Commands

### `type_all`

Runs the full workflow (`check_depth` → `kmer_count` → `cen_type` for all chromosomes) in a single command. Automatically determines biological sex and applies single-haplotype mode for chrX in males.

```bash
ascairn type_all \
    seq_data/NA12877.final.cram \
    -o output/NA12877 \
    --resource_dir ascairn_resource/resource/panel/ascairn_paper_2025 \
    --reference hg38 \
    -t 8
```

| Option | Description | Required |
|--------|-------------|----------|
| `BAM_FILE` | Path to BAM or CRAM file (positional argument) | Yes |
| `-o` | Output path prefix | Yes |
| `--resource_dir` | Path to panel resource directory | Yes |
| `--reference` | Reference genome build (`hg38` or `chm13`) | Yes |
| `-t` | Number of threads | No (default: 8) |

### Individual commands

The workflow consists of three commands that are normally run in order.
You can run them individually for more control.

### `check_depth`

Checks sequence coverage in a reference region (chr22 long arm) and determines biological sex by comparing chrX coverage.

```bash
ascairn check_depth \
    seq_data/NA12877.final.cram \
    -o output/NA12877.depth.txt \
    --baseline_region ascairn_resource/resource/common/chr22_long_arm_hg38.bed \
    --x_region ascairn_resource/resource/common/chrX_short_arm_hg38.bed \
    -t 8
```

**Output** (`output/NA12877.depth.txt`):
```
Coverage: 37.24
Baseline region file: ascairn_resource/resource/common/chr22_long_arm_hg38.bed
Sex: male
ChrX region file: ascairn_resource/resource/common/chrX_short_arm_hg38.bed
ChrX coverage: 18.0
```

### `kmer_count`

Extracts reads aligned to alpha satellite regions and counts occurrences of predefined rare k-mers.

```bash
ascairn kmer_count \
    seq_data/NA12877.final.cram \
    -o output/NA12877.kmer_count.txt \
    --kmer_file ascairn_resource/resource/panel/ascairn_paper_2025/rare_kmer_list.fa \
    --cen_region ascairn_resource/resource/common/cen_region_hg38.bed \
    -t 8
```

**Output** (`output/NA12877.kmer_count.txt`): a two-column TSV with k-mer sequences and their counts.

### `cen_type`

Identifies the most likely centromere cluster pair and nearest haplotype pair for a given chromosome. Requires the depth file from `check_depth` output.

```bash
ascairn cen_type \
    output/NA12877.kmer_count.txt \
    -o output/NA12877.chr22 \
    --kmer_info ascairn_resource/resource/panel/ascairn_paper_2025/kmer_info/chr22.kmer_info.txt.gz \
    --hap_info ascairn_resource/resource/panel/ascairn_paper_2025/hap_info/chr22.hap_info.txt \
    --depth_file output/NA12877.depth.txt
```

For male samples on chrX, add the `--single_hap` option (`type_all` handles this automatically based on `check_depth` output).

**Output files** (using `output/NA12877.chr22` as the output prefix):

| File | Description |
|------|-------------|
| `*.cen_type.txt` | Best cluster pair, haplotype pair, and hap_info annotations (e.g., Contig_len) |
| `*.cluster.hap_pair.txt` | Ranked cluster pairs with log-likelihoods. The first data row is the best pair. |
| `*.haplotype.hap_pair.txt` | Ranked haplotype pairs with log-likelihoods. The first data row is the best pair. |
| `*.cluster.marker_prob.txt` | Per-marker copy number probabilities used for cluster assignment |
| `*.haplotype.marker_prob.txt` | Per-marker copy number probabilities used for haplotype assignment |

## Output Format

### `cen_type_all.txt`

The main result file (generated by `type_all`) is a TSV with one row per chromosome:

| Column | Description |
|--------|-------------|
| Chr | Chromosome (chr1-chr22, chrX, chrY) |
| Cluster_1 | Best-matching cluster ID for the first allele |
| Cluster_2 | Best-matching cluster ID for the second allele |
| Haplotype_1 | Nearest haplotype ID for the first allele (e.g., `HG00268.mat`) |
| Haplotype_2 | Nearest haplotype ID for the second allele |
| *additional columns* | Any extra columns from `hap_info.txt` are automatically appended with `_1`/`_2` suffixes (e.g., `Contig_len_1`, `Contig_len_2`) |

For male chrX, `_2` columns are `NA`.

### Per-chromosome `cen_type.txt`

Each `cen_type` run produces a `*.cen_type.txt` file with the same columns as above (without the `Chr` column). The `cen_type_all.txt` file is the concatenation of these per-chromosome files.

### Cluster and haplotype pair files

These files list all candidate pairs ranked by log-likelihood (higher = better fit):

| Column | Description |
|--------|-------------|
| Cluster1 / Haplotype1 | First member of the pair |
| Cluster2 / Haplotype2 | Second member of the pair |
| Loglikelihood | Log-likelihood of the pair given the observed k-mer counts |

## Performance

We evaluated ascairn using leave-one-individual-out cross-validation on the reference panel: for each individual, both parental aHOR-haps were removed from the panel, and the resulting model was used to infer the haplogroup pair from the individual's short-read WGS data. To assess robustness to sequencing depth, each sample was downsampled to 1–30x coverage.

The figures below show haplogroup pair assignment accuracy (precision) for each chromosome at various downsampled coverages, for short-read WGS data aligned to GRCh38 (upper) and CHM13 v2.0 (lower):

**GRCh38-aligned data:**

<div align="center">
  <img src="image/cluster_match_hg38.png" alt="Accuracy on GRCh38-aligned data" width="750">
</div>

**CHM13-aligned data:**

<div align="center">
  <img src="image/cluster_match_chm13.png" alt="Accuracy on CHM13-aligned data" width="750">
</div>

Accuracy is generally high (>90%) for most chromosomes at coverage ≥ 5x, with comparable performance between GRCh38- and CHM13-aligned data. See [Shiraishi et al., bioRxiv, 2025](https://doi.org/10.1101/2025.07.26.666712) for details.

## Notes

- **CHM13 support**: ascairn supports both GRCh38 (hg38) and T2T-CHM13 (chm13) aligned BAM/CRAM files. Specify the reference genome via `--reference hg38` or `--reference chm13`. The centromere region BED files for each reference are provided in `resource/common/`.
- **chrY**: chrY is processed only for male samples (determined automatically by `check_depth` from the chrX coverage ratio).
- **Test data**: The Quick Start example uses NA12877 from the 1000 Genomes Project high-coverage dataset. This sample was chosen because it is publicly accessible via both AWS S3 and FTP, enabling reproducible testing without restricted data access.

## Citation

Shiraishi Y, Ochi Y, Sugawa M, Sakamoto Y, Kimura K, Tsujimura T, Okada A, Okuda R, Namba S, Miyauchi T, Mateos RN, Suzuki H, Chiba K, Ito Y, Nakamura W, Ohka F, Motomura K, Yamamoto T, Kawai Y, Okada Y, Suzuki H, Kato M, Saito R, Garrison E, Logsdon GA, Ogawa S. Rare k-mers reveal centromere haplogroups underlying human diversity and cancer translocations. *bioRxiv*. 2025. doi: [10.1101/2025.07.26.666712](https://doi.org/10.1101/2025.07.26.666712)
