import click
import os
import subprocess

from ascairn.logger import get_logger
logger = get_logger(__name__)


def read_sex_from_depth_file(depth_file):
    """Parse Sex value from depth.txt file."""
    with open(depth_file) as f:
        for line in f:
            if line.startswith("Sex:"):
                return line.split()[1]
    return None


@click.command()
@click.argument("bam_file", type=click.Path(exists=False))
@click.option("-o", "--output_prefix", required=True, type=click.Path())
@click.option("--resource_dir", required=True, type=click.Path(exists=True),
              help="Path to panel resource directory (e.g. ascairn_resource/resource/panel/ascairn_paper_2025).")
@click.option("--reference", required=True, type=click.Choice(["hg38", "chm13"]),
              help="Reference genome build for the BAM file.")
@click.option("-t", "--threads", default=8, help="Number of threads to use.")
@click.option("--debug", is_flag=True, default=False,
              help="Keep per-chromosome intermediate files (cen_type.txt, cluster/haplotype tables).")
def type_all_command(bam_file, output_prefix, resource_dir, reference, threads, debug):
    """Run the full ascairn workflow (check_depth, kmer_count, cen_type for all chromosomes)."""

    os.environ["POLARS_MAX_THREADS"] = str(threads)

    # make directory for the output prefix
    output_dir = os.path.dirname(output_prefix)
    if output_dir != '' and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    common_dir = os.path.join(resource_dir, "..", "..", "common")

    depth_file = f"{output_prefix}.depth.txt"
    kmer_count_file = f"{output_prefix}.kmer_count.txt"

    # Step 1: check_depth
    logger.info("Step 1: Checking sequence depth")
    subprocess.run([
        "ascairn", "check_depth",
        bam_file,
        "-o", depth_file,
        "--baseline_region", os.path.join(common_dir, f"chr22_long_arm_{reference}.bed"),
        "--x_region", os.path.join(common_dir, f"chrX_short_arm_{reference}.bed"),
        "-t", str(threads),
    ], check=True)

    # Step 2: kmer_count
    logger.info("Step 2: Counting rare k-mers")
    subprocess.run([
        "ascairn", "kmer_count",
        bam_file,
        "-o", kmer_count_file,
        "--kmer_file", os.path.join(resource_dir, "rare_kmer_list.fa"),
        "--cen_region", os.path.join(common_dir, f"cen_region_curated_margin_{reference}.bed"),
        "-t", str(threads),
    ], check=True)

    # Read sex from depth file
    sex = read_sex_from_depth_file(depth_file)

    # Step 3: cen_type for each chromosome
    chromosomes = [str(i) for i in range(1, 23)] + ["X", "Y"]
    first_chr = True

    for chrom in chromosomes:
        # chrY only exists in males
        if chrom == "Y" and sex != "male":
            continue

        chr_prefix = f"{output_prefix}.chr{chrom}"
        logger.info(f"Step 3: Running cen_type for chr{chrom}")

        cmd = [
            "ascairn", "cen_type",
            kmer_count_file,
            "-o", chr_prefix,
            "--kmer_info", os.path.join(resource_dir, "kmer_info", f"chr{chrom}.kmer_info.txt.gz"),
            "--hap_info", os.path.join(resource_dir, "hap_info", f"chr{chrom}.hap_info.txt"),
            "--depth_file", depth_file,
        ]
        if chrom in ("X", "Y") and sex == "male":
            cmd.append("--single_hap")

        subprocess.run(cmd, check=True)

        # Aggregate per-chromosome cen_type.txt
        result_file = f"{output_prefix}.cen_type_all.txt"
        chr_result = f"{chr_prefix}.cen_type.txt"

        if first_chr:
            with open(chr_result) as f_in, open(result_file, 'w') as f_out:
                header = f_in.readline().rstrip('\n')
                f_out.write(f"Chr\t{header}\n")
                data = f_in.readline().rstrip('\n')
                f_out.write(f"chr{chrom}\t{data}\n")
            first_chr = False
        else:
            with open(chr_result) as f_in, open(result_file, 'a') as f_out:
                f_in.readline()  # skip header
                data = f_in.readline().rstrip('\n')
                f_out.write(f"chr{chrom}\t{data}\n")

        # Remove per-chromosome intermediate files unless --debug is set
        if not debug:
            for suffix in ("cen_type.txt", "cluster.hap_pair.txt", "cluster.marker_prob.txt",
                           "haplotype.hap_pair.txt", "haplotype.marker_prob.txt"):
                path = f"{chr_prefix}.{suffix}"
                if os.path.exists(path):
                    os.remove(path)

    logger.info(f"Completed. Results written to {output_prefix}.cen_type_all.txt")
