import click
import os
import importlib.resources as pkg_resources
from ascairn.utils import *
from ascairn.match import *
# from ascairn.utils import bam_processing, dummy_scripts

@click.command()
@click.argument("bam_file", type=click.Path(exists=True))
@click.argument("output_prefix", type=click.Path())
@click.option("--reference", type=click.Choice(["hg38", "chm13"]), default="hg38", help="Reference genome (hg38 or chm13).")
@click.option("--threads", default=4, help="Number of threads to use.")
# @click.option("--kmer_file", type=click.Path(), default="ascairn/data/chr22.rare_kmer.pruned.annot.long.txt")
# @click.option("--cluster_file", type=click.Path(), default="ascairn/data/chr22.cluster_marker_count.txt")
# @click.option("--hap_file", type=click.Path(), default="ascairn/data/chr22.hap_cluster.txt")
@click.option("--baseline_region_file", type=click.Path(), default=None)
@click.option("--cen_region_file", type=click.Path(), default=None)
def type_command(bam_file, output_prefix, reference, threads, baseline_region_file, cen_region_file): # kmer_file, cluster_file, hap_file):

    """
    if kmer_file is None:
        with pkg_resources.path("ascairn.data", "chr22.rare_kmer.pruned.annot.long.txt") as default_kmer_file:
            kmer_file = str(default_kmer_file)
    if cluster_file is None:
        with pkg_resources.path("ascairn.data", "chr22.cluster_marker_count.txt") as default_cluster_file:
            cluster_file = str(default_cluster_file)

    if hap_file is None:
        with pkg_resources.path("ascairn.data", "chr22.hap_cluster.txt") as default_hap_file:
            hap_file = str(default_hap_file)
    """

    # check if the executables exist
    is_tool("samtools")
    is_tool("jellyfish")
    # is_tool("mosdepth")

    # check input file existences
    is_exists_bam(bam_file)

    # make directory for the output prefix
    output_dir = os.path.dirname(output_prefix)
    if output_dir != '' and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    ##########
    # preparing several preset files
    if baseline_region_file is None:
        if reference == "hg38":
            with pkg_resources.path("ascairn.data", "chr22_long_arm_hg38.bed") as default_baseline_region_file:
                baseline_region_file = str(default_baseline_region_file)
        else:
            with pkg_resources.path("ascairn.data", "chr22_long_arm_chm13.bed") as default_baseline_region_file:
                baseline_region_file = str(default_baseline_region_file)

    if cen_region_file is None:
        if reference == "hg38":
            with pkg_resources.path("ascairn.data", "cen_region_curated_margin_hg38.bed") as default_cen_region_file:
                cen_region_file = str(default_cen_region_file)
        else:
            with pkg_resources.path("ascairn.data", "cen_region_curated_margin_chm13.bed") as default_cen_region_file:
                cen_region_file = str(default_cen_region_file)

    with pkg_resources.path("ascairn.data", "rare_kmer_list.fa") as default_rare_kmer_file:
        rare_kmer_file = str(default_rare_kmer_file)

    with pkg_resources.path("ascairn.data", "rare_kmer_list.fa") as default_rare_kmer_file:
        rare_kmer_file = str(default_rare_kmer_file)

    ##########


    depth = check_depth(bam_file, output_prefix, baseline_region_file)

    gather_rare_kmer(bam_file, output_prefix, cen_region_file, rare_kmer_file, kmer_size = 27, num_threads = threads)


    match_cluster_haplotype(output_prefix + ".kmer_count.txt", output_prefix + ".match/chr" + cen_id, 
        kmer_info_file, cluster_kmer_count_file, depth,
        cluster_haplotype_file, cluster_ratio = 0.1, pseudo_count = 0.1, nbinom_size_0 = 0.5, nbinom_size = 8, nbinom_mu_0 = 0.8, nbinom_mu_unit = 0.4):

