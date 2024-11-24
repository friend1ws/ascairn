import subprocess, sys, os, csv, shutil
# import pysam

from ascairn.logger import get_logger
logger = get_logger(__name__)

def is_exists_bam(input_file):

    if input_file.startswith("s3://"):
        is_exists_s3(input_file)
    else:
        is_exists(input_file)

def is_exists(input_file):
    
    if not os.path.exists(input_file):
        logger.error("Input not exists: %s" % input_file)
        sys.exit(1)


def is_tool(executable):

    from shutil import which
    if which(executable) is None:
        logger.error("Executable does not exist: " + executable)
        sys.exit(1) 

    return True


def is_exists_s3(bam_object):

    from urllib.parse import urlparse
    import boto3

    obj_p = urlparse(bam_object)

    tbucket = obj_p.netloc
    tkey = obj_p.path
    if tkey.startswith("/"): tkey = tkey[1:]

    client = boto3.client("s3")
    try:
        response = client.head_object(Bucket = tbucket, Key = tkey)
    except:
        logger.error("Input not exists: %s: " % bam_object)
        sys.exit(1)




def check_depth(bam_file, output_file, baseline_region_file, num_threads = 4):

    # make directory for the output file 
    tmp_dir = output_file + ".tmp_dir.check_depth"
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
 
    # 出力ファイルパスの設定
    baseline_bam = tmp_dir + "/baseline.bam"
    mosdepth_prefix = tmp_dir + "/baseline"
    mosdepth_summary = f'{mosdepth_prefix}.mosdepth.summary.txt'

    subprocess.run(["samtools", "view", "-bh", bam_file, "-L", baseline_region_file, "-M", "-@", str(num_threads), "-o", baseline_bam], check=True)
    
    subprocess.run(["samtools", "index", baseline_bam], check=True)
    
    subprocess.run(["mosdepth", mosdepth_prefix, baseline_bam, "-b", baseline_region_file, "-t", str(num_threads)], check=True)

    depth = None
    with open(mosdepth_summary, 'r') as hin:
        for F in csv.DictReader(hin, delimiter = '\t'):
            if F["chrom"] == "total_region": depth = float(F["mean"])

    with open(output_file, 'w') as hout:
        print(depth, file = hout)
    
    shutil.rmtree(tmp_dir)

    return depth 

 


def gather_rare_kmer(bam_file, output_prefix, cen_region_file, rare_kmer_file, kmer_size = 27, num_threads = 4):

    tmp_dir = output_prefix + ".tmp_dir.gather_rare_kmer"
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    tmp_bam = tmp_dir + "/centromere.bam"
    tmp_fasta = tmp_dir + "/centromere.fasta"
    tmp_kmer_jf = tmp_dir + "/centromere.rare_kmer.jf"
    tmp_kmer_fa = tmp_dir + "/centromere.rare_kmer.fa"
    kmer_count_file = output_prefix + ".kmer_count.txt"

    subprocess.run(["samtools", "view", "-bh", bam_file, "-L", cen_region_file, "-M", "-@", str(num_threads), "-o", tmp_bam], check=True)

    with open(tmp_fasta, 'w') as hout:
        subprocess.run(["samtools", "fasta", "-@", str(num_threads), tmp_bam], stdout=hout, stderr=subprocess.DEVNULL, check=True)

    subprocess.run(["jellyfish", "count", "-s", "100M", "-C", "-m", str(kmer_size), "-t", str(num_threads), "--if", rare_kmer_file, "-o", tmp_kmer_jf, tmp_fasta], check = True) 
        
    with open(tmp_kmer_fa, 'w') as hout:
        subprocess.run(["jellyfish", "dump", tmp_kmer_jf], stdout=hout, check=True)

    
    kmer = None
    count = None
    with open(tmp_kmer_fa, 'r') as hin, open(kmer_count_file, 'w') as hout:
        for line in hin:
            line = line.rstrip('\n')
            if line.startswith('>'): 
                if kmer is not None: 
                    print(f'{kmer}\t{count}', file = hout)
                count = line.lstrip('>').rstrip('\n')
            else:
                kmer = line.rstrip('\n')

        if kmer is not None:
            print(f'{kmer}\t{count}', file = hout)
    
    shutil.rmtree(tmp_dir)


