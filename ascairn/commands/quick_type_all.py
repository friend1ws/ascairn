import click
import os
import importlib.resources
from ascairn.utils import *
from ascairn.match import *
# from ascairn.utils import bam_processing, dummy_scripts

from ascairn.logger import get_logger
logger = get_logger(__name__)

@click.command()
@click.argument("bam_file", type=click.Path(exists=False))
@click.argument("output_prefix", type=click.Path())
@click.option("--reference", type=click.Choice(["hg38", "chm13"]), default="hg38", help="Reference genome (hg38 or chm13).")
@click.option("--sex", type=click.Choice(["male", "female", "unknown"]), default="unknown", help="Sex (male, femail or unknown).")
@click.option("--threads", default=4, help="Number of threads to use.")
# @click.option("--kmer_file", type=click.Path(), default="ascairn/data/chr22.rare_kmer.pruned.annot.long.txt")
# @click.option("--cluster_file", type=click.Path(), default="ascairn/data/chr22.cluster_marker_count.txt")
# @click.option("--hap_file", type=click.Path(), default="ascairn/data/chr22.hap_cluster.txt")
@click.option("--baseline_region_file", type=click.Path(), default=None)
@click.option("--cen_region_file", type=click.Path(), default=None)
def quick_type_all_command(bam_file, output_prefix, reference, sex, threads, baseline_region_file, cen_region_file): # kmer_file, cluster_file, hap_file):


    # check if the executables exist
    is_tool("samtools")
    is_tool("jellyfish")
    is_tool("mosdepth")

    # check input file existences
    is_exists_bam(bam_file)

    # make directory for the output prefix
    output_dir = os.path.dirname(output_prefix)
    if output_dir != '' and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if not os.path.exists(output_prefix + ".match"):
        os.makedirs(output_prefix + ".match")

    ##########
    # preparing several preset files
    if baseline_region_file is None:
        if reference == "hg38":
            baseline_region_file = importlib.resources.files("ascairn.data").joinpath("chr22_long_arm_hg38.bed")
        else:
            baseline_region_file = importlib.resources.files("ascairn.data").joinpath("chr22_long_arm_chm13.bed")

    if cen_region_file is None:
        if reference == "hg38":
            cen_region_file = importlib.resources.files("ascairn.data").joinpath("cen_region_curated_margin_hg38.bed")
        else:
            cen_region_file = importlib.resources.files("ascairn.data").joinpath("cen_region_curated_margin_chm13.bed")

    rare_kmer_file = importlib.resources.files("ascairn.data").joinpath("rare_kmer_list.fa")
    ##########

    logger.info("Checking the sequence coverage in the baseline region") 
    depth = check_depth(bam_file, output_prefix + ".baseline.depth.txt", baseline_region_file, num_threads = threads)
    
    if sex == "unknown":
        if reference == "hg38":
            chrX_short_region_file = importlib.resources.files("ascairn.data").joinpath("chrX_short_arm_hg38.bed")
        else:
            chrX_short_region_file = importlib.resources.files("ascairn.data").joinpath("chrX_short_arm_chm13.bed")

        logger.info("Checking the sequence coverage in the short arm of chromosome X short-arm region for sex determination")
        depth_chrX = check_depth(bam_file, output_prefix + ".chrX_short.depth.txt", chrX_short_region_file, num_threads = threads)

        logger.info("The depth ratio: %s" % str(depth_chrX / depth))
        if depth_chrX / depth < 0.8: 
            logger.info("Sex has been set to male")
            sex = "male"
        else:
            logger.info("Sex has been set to female")
            sex = "female"

    # logger.info("Parsing rare kmer from the BAM file")
    # gather_rare_kmer(bam_file, output_prefix, cen_region_file, rare_kmer_file, kmer_size = 27, num_threads = threads)
    count_rare_kmer(bam_file, output_prefix + ".kmer_count.txt", cen_region_file, kmer_file_fasta, kmer_size, threads)

    # depth = 41.0
    for cen_id in [str(x) for x in range(1, 23)] + ["X"]:  

        kmer_info_file = importlib.resources.files("ascairn.data").joinpath(f"kmer_info/chr{cen_id}.kmer_info.txt.gz")
        cluster_kmer_count_file = importlib.resources.files("ascairn.data").joinpath(f"cluster/chr{cen_id}.cluster_marker_count.txt.gz")
        cluster_haplotype_file = importlib.resources.files("ascairn.data").joinpath(f"cluster/chr{cen_id}.hap_cluster.txt")

        
        if cen_id == "X" and sex == "male":

            logger.info("Matching haplotypes of chr%s" % cen_id)
            match_cluster_haplotype_single(output_prefix + ".kmer_count.txt", output_prefix + ".match/chr" + cen_id,
                kmer_info_file, cluster_kmer_count_file, depth,
                cluster_haplotype_file, cluster_ratio = 0.1, pseudo_count = 0.1, nbinom_size_0 = 0.5, nbinom_size = 8, nbinom_mu_0 = 0.8, nbinom_mu_unit = 0.4)
        else:

            logger.info("Matching haplotypes of chr%s" % cen_id)
            match_cluster_haplotype(output_prefix + ".kmer_count.txt", output_prefix + ".match/chr" + cen_id, 
                kmer_info_file, cluster_kmer_count_file, depth,
                cluster_haplotype_file, cluster_ratio = 0.1, pseudo_count = 0.1, nbinom_size_0 = 0.5, nbinom_size = 8, nbinom_mu_0 = 0.8, nbinom_mu_unit = 0.4)

        
