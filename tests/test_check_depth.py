# Test for ascairn check_depth command
import os
import subprocess


def compare_depth_values(output_file, expected_file):
    """Compare depth files ignoring file path lines (which vary by environment)."""
    path_prefixes = ("Baseline region file:", "ChrX region file:")
    with open(output_file) as f1, open(expected_file) as f2:
        out_lines = [l.rstrip("\n") for l in f1 if not l.startswith(path_prefixes)]
        exp_lines = [l.rstrip("\n") for l in f2 if not l.startswith(path_prefixes)]
    return out_lines == exp_lines


def test_check_depth(s3_cram, common_dir, output_dir, expected_dir):
    output_file = os.path.join(output_dir, "NA12877.depth.txt")
    expected_file = os.path.join(expected_dir, "NA12877.depth.txt")

    subprocess.run(
        [
            "ascairn", "check_depth",
            s3_cram,
            "-o", output_file,
            "--baseline_region", os.path.join(common_dir, "chr22_long_arm_hg38.bed"),
            "--x_region", os.path.join(common_dir, "chrX_short_arm_hg38.bed"),
            "-t", "4",
        ],
        check=True,
    )

    assert os.path.exists(output_file)
    assert compare_depth_values(output_file, expected_file)
