# Test for ascairn kmer_count command
import os
import subprocess
import filecmp


def test_kmer_count(s3_cram, common_dir, panel_dir, output_dir, expected_dir):
    output_file = os.path.join(output_dir, "NA12877.kmer_count.txt")
    expected_file = os.path.join(expected_dir, "NA12877.kmer_count.txt")

    subprocess.run(
        [
            "ascairn", "kmer_count",
            s3_cram,
            "-o", output_file,
            "--kmer_file", os.path.join(panel_dir, "rare_kmer_list.fa"),
            "--cen_region", os.path.join(common_dir, "cen_region_curated_margin_hg38.bed"),
            "-t", "4",
        ],
        check=True,
    )

    assert os.path.exists(output_file)
    assert filecmp.cmp(output_file, expected_file, shallow=False)
