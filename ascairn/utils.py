import subprocess, sys, os, csv, shutil
# import pysam

from ascairn.logger import get_logger
logger = get_logger(__name__)


"""
def is_exists_bam(input_file):

    if input_file.startswith("s3://"):
        is_exists_s3(input_file)
    else:
        is_exists(input_file)

def is_exists(input_file):
    
    if not os.path.exists(input_file):
        logger.error("Input not exists: %s" % input_file)
        sys.exit(1)
"""

def is_exists_bam(input_file):

    logger.info("Checking accessibility of the input BAM file.")
    try:
        subprocess.run(["samtools", "view", "-H", input_file], check = True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    except Exception as e:
        logger.error(
            "Failed to access the file '%s' using samtools. Please verify the following: "
            "1) The path of the BAM file is correct, "
            "2) The BAM file and its index (.bai) are accessible, and "
            "3) Samtools is properly installed and configured." % input_file
        )
        sys.exit(1)
        
def is_tool(executable):

    from shutil import which
    if which(executable) is None:
        logger.error("Executable does not exist: " + executable)
        sys.exit(1) 

    return True

"""
def is_configured_aws():

    import boto3
    from botocore.exceptions import NoCredentialsError, PartialCredentialsError

    # try:
    #     session = boto3.Session()
    #     credentials = session.get_credentials()

    #     if credentials is None:
    #         return False
    #     if credentials.access_key and credentials.secret_key:
    #         return True
    #     else:
    #         return False
    # except (NoCredentialsError, PartialCredentialsError):
    #     return False

    try:
        session = boto3.Session()
        credentials = session.get_credentials()

        if credentials is None or not (credentials.access_key and credentials.secret_key):
            return False

        sts = session.client('sts')
        sts.get_caller_identity()
        return True

    except (NoCredentialsError, PartialCredentialsError, ClientError):
        return False


def is_exists_s3(bam_object):

    from urllib.parse import urlparse
    import boto3

    obj_p = urlparse(bam_object)

    tbucket = obj_p.netloc
    tkey = obj_p.path
    if tkey.startswith("/"): tkey = tkey[1:]

    if is_configured_aws():
        logger.info("Using configured AWS credentials.")
        client = boto3.client("s3")
    else:
        logger.info("No AWS credentials found. Using unsigned client.")
        config = Config(signature_version='UNSIGNED')
        client = boto3.client("s3", config=config)

    try:
        response = client.head_object(Bucket = tbucket, Key = tkey)
    except:
        logger.error("Input not exists: %s: " % bam_object)
        sys.exit(1)
"""



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

 


def gather_rare_kmer(bam_file, output_file, cen_region_file, rare_kmer_file, kmer_size = 27, num_threads = 4):

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



def count_rare_kmer(bam_file, output_file, cen_region_file, rare_kmer_file, kmer_size = 27, num_threads = 4):

    output_dir = os.path.dirname(output_file)
    if output_dir != '' and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    tmp_bam = output_file + ".tmp.centromere.bam"
    tmp_fasta = output_file + ".centromere.fasta"
    tmp_kmer_jf = output_file + ".centromere.rare_kmer.jf"
    tmp_kmer_fa = output_file + ".centromere.rare_kmer.fa"

    logger.info("Extracting reads aligned to the target centromere region from the BAM file.")
    subprocess.run(["samtools", "view", "-bh", bam_file, "-L", cen_region_file, "-M", "-@", str(num_threads), "-o", tmp_bam], check=True)

    logger.info("Generating a FASTA format file of extracted reads for Jellyfish execution.")
    with open(tmp_fasta, 'w') as hout:
        subprocess.run(["samtools", "fasta", "-@", str(num_threads), tmp_bam], stdout=hout, stderr=subprocess.DEVNULL, check=True)

    logger.info("Counting rare k-mers using Jellyfish.")
    subprocess.run(["jellyfish", "count", "-s", "100M", "-C", "-m", str(kmer_size), "-t", str(num_threads), "--if", rare_kmer_file, "-o", tmp_kmer_jf, tmp_fasta], check = True) 
        
    logger.info("Organizing the result of Jellyfish.")
    with open(tmp_kmer_fa, 'w') as hout:
        subprocess.run(["jellyfish", "dump", tmp_kmer_jf], stdout=hout, check=True)

    
    kmer = None
    count = None
    with open(tmp_kmer_fa, 'r') as hin, open(output_file, 'w') as hout:
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
    
    os.remove(tmp_bam)
    os.remove(tmp_fasta)
    os.remove(tmp_kmer_jf)
    os.remove(tmp_kmer_fa)


def convert_tsv_to_fasta(kmer_file_tsv, kmer_file_fasta):

    with open(kmer_file_tsv, 'r') as hin, open(kmer_file_fasta, 'w') as hout:
        tind = 0
        for line in hin:
            F = line.rstrip('\n').split('\t')
            print(f'>kmer_{tind}\n{F[0]}', file = hout)
            tind = tind + 1


def check_kmer_size_from_kmer_fasta(kmer_file_fasta):

    tmp_kmer_size = None
    tind = 0
    with open(kmer_file_fasta, 'r') as hin:
        for line in hin:
            if line.startswith('>'): continue
            if tmp_kmer_size is not None and tmp_kmer_size != len(line.rstrip('\n')):
                logger.error("Sizes of rare kmers are inconsistent!")
                sys.exit(1)
            tmp_kmer_size = len(line.rstrip('\n'))
            tind = tind + 1
            if tind >= 1000: break

    if tmp_kmer_size is None:
        logger.info("Could not read the rare kmers")
        sys.exit(1)

    return(tmp_kmer_size)

            
