import click
import os
import importlib.resources
from ascairn.utils import *
from ascairn.match import *
# from ascairn.utils import bam_processing, dummy_scripts

from ascairn.logger import get_logger
logger = get_logger(__name__)

@click.command()
@click.argument("kmer_count_file", type=click.Path(exists=True))
@click.argument("output_prefix", type=click.Path())
@click.argument("kmer_info_file", type=click.Path(exists=True))
@click.argument("cluster_kmer_count_file", type=click.Path(exists=True))
@click.argument("depth", type=float)
@click.option("--cluster_haplotype_file", type=click.Path(exists=True))
@click.option("--is_single_hap", is_flag=True, show_default=True, default=False)
@click.option("--cluster_ratio", default=0.1)
@click.option("--pseudo_count", default=0.1)
@click.option("--nbinom_size_0", default=0.5)
@click.option("--nbinom_size", default=8.0)
@click.option("--nbinom_mu_0", default=0.8)
@click.option("--nbinom_mu_unit", default=0.4)
def cen_type_command(kmer_count_file, output_prefix, kmer_info_file, cluster_kmer_count_file, depth, cluster_haplotype_file, is_single_hap,
                 cluster_ratio = 0.1, pseudo_count = 0.1, nbinom_size_0 = 0.5, nbinom_size = 8, nbinom_mu_0 = 0.8, nbinom_mu_unit = 0.4): 


    # make directory for the output prefix
    output_dir = os.path.dirname(output_prefix)
    if output_dir != '' and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    logger.info("Matching haplotypes")
    if is_single_hap == True:

        match_cluster_haplotype_single(kmer_count_file, output_prefix, kmer_info_file, cluster_kmer_count_file, depth,
            cluster_haplotype_file, cluster_ratio = 0.1, pseudo_count = 0.1, nbinom_size_0 = 0.5, nbinom_size = 8, nbinom_mu_0 = 0.8, nbinom_mu_unit = 0.4)
    else:

        match_cluster_haplotype(kmer_count_file, output_prefix, kmer_info_file, cluster_kmer_count_file, depth,
            cluster_haplotype_file, cluster_ratio = 0.1, pseudo_count = 0.1, nbinom_size_0 = 0.5, nbinom_size = 8, nbinom_mu_0 = 0.8, nbinom_mu_unit = 0.4)

    logger.info("Completed.")
 
