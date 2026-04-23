# Test for ascairn cen_type command
import os
import subprocess
import filecmp
import pytest

# Test parameters: (chromosome, is_single_hap)
# NA12877 is male, so chrX uses --single_hap
# Note: chr1 excluded due to high memory requirement (~6GB)
CHROMOSOMES = [
    ("22", False),
    ("X", True),
]


@pytest.mark.parametrize("chrom,is_single_hap", CHROMOSOMES)
def test_cen_type(chrom, is_single_hap, panel_dir, output_dir, expected_dir):
    # Use expected files as input (independent of other tests)
    kmer_count_file = os.path.join(expected_dir, "NA12877.kmer_count.txt")
    depth_file = os.path.join(expected_dir, "NA12877.depth.txt")

    output_prefix = os.path.join(output_dir, f"NA12877.chr{chrom}")

    cmd = [
        "ascairn", "cen_type",
        kmer_count_file,
        "-o", output_prefix,
        "--kmer_info", os.path.join(panel_dir, "kmer_info", f"chr{chrom}.kmer_info.txt.gz"),
        "--hap_info", os.path.join(panel_dir, "hap_info", f"chr{chrom}.hap_info.txt"),
        "--depth_file", depth_file,
    ]
    if is_single_hap:
        cmd.append("--single_hap")

    subprocess.run(cmd, check=True, env=os.environ)

    output_file = f"{output_prefix}.cen_type.txt"
    expected_file = os.path.join(expected_dir, f"NA12877.chr{chrom}.cen_type.txt")
    assert os.path.exists(output_file), f"Output file not found: {output_file}"
    assert filecmp.cmp(output_file, expected_file, shallow=False), "Mismatch: cen_type.txt"
