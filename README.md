<div align="center">
  <img src="image/karamatsu_cairn.png" alt="Cairn at Mt. Karamatsu" width="250">
  <p>Maruyama cairn at Mt. Karamatsu</p>
</div>

# ascairn
Centromere sequence analysis using rare k-mer markers

# Dependency
## Software
- [samtools](https://github.com/samtools/samtools)
- [Jellyfish](https://github.com/gmarcais/Jellyfish)
- [mosdepth](https://github.com/brentp/mosdepth)

## Python
- click
- boto3 (required only for accessing CRAM files in Amazon S3)
  
# Install

1. Install all the prerequisite software and install `ascairn`.
```
pip install ascairn (--user)
```

2. Download the resource files for ascairn (in preparation).

# Mini-tutorial

1. Prepare the sequence data

2. Calculate the sequence depth and determine the sex for the sample.

3. Parse the rare-kmer counts of the target BAM file in the alpha satellite region.

4. 
5. Perform centromere typing on a local BAM/CRAM file:
```
ascairn type /path/to/local/file.bam /path/to/output
```
Example command using a CRAM file located in Amazon S3 (it will take 1~2 hours):
```
ascairn type s3://1000genomes/1000G_2504_high_coverage/data/ERR3240114/HG00096.final.cram output/HG00096 --threads 8
```
