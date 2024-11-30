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
@click.argument("output_file", type=click.Path())
@click.option("--reference", type=click.Choice(["hg38", "chm13"]), default="hg38", help="Reference genome (hg38 or chm13).")
@click.option("--kmer_file", type=click.Path())
@click.option("--threads", default=4, help="Number of threads to use.")
# @click.option("--kmer_file", type=click.Path(), default="ascairn/data/chr22.rare_kmer.pruned.annot.long.txt")
@click.option("--cen_region_file", type=click.Path(), default=None)
def kmer_count_command(bam_file, kmer_file, output_file, reference, threads, cen_region_file): # kmer_file, cluster_file, hap_file):


    # check if the executables exist
    is_tool("samtools")
    is_tool("jellyfish")

    # check input file existences
    is_exists_bam(bam_file)

    # make directory for the output prefix
    output_dir = os.path.dirname(output_file)
    if output_dir != '' and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if cen_region_file is None:
        if reference == "hg38":
            cen_region_file = importlib.resources.files("ascairn.data").joinpath("cen_region_curated_margin_hg38.bed")
        else:
            cen_region_file = importlib.resources.files("ascairn.data").joinpath("cen_region_curated_margin_chm13.bed")

    is_make_kmer_file_fasta = False
    if kmer_file is None:
        kmer_file_fasta = importlib.resources.files("ascairn.data").joinpath("rare_kmer_list.fa")
        logger.info("No kmer_file argument provided. Using the default rare k-mer list file.")
    elif kmer_file.endswith(".fa") or kmer_file.endswith(".fasta"):
        logger.info("Since the kmer_file has a suffix of '.fa' or '.fasta', consider it as FASTA format file.")
        kmer_file_fasta = kmer_file
    else:
        logger.info("Convert the TSV format kmer_file to FASTA format")
        convert_tsv_to_fasta(kmer_file, output_file + ".tmp.kmer_list.fa")
        kmer_file_fasta = output_file + ".tmp.kmer_list.fa"
        is_make_kmer_file_fasta = True

    kmer_size = check_kmer_size_from_kmer_fasta(kmer_file_fasta)

    logger.info("Counting rare kmer from the BAM file")
    count_rare_kmer(bam_file, output_file, cen_region_file, kmer_file_fasta, kmer_size, threads)

    if is_make_kmer_file_fasta:
        os.remove(output_file + ".tmp.kmer_list.fa")



