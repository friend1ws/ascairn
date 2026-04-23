# pytest fixtures for ascairn tests
import os
import subprocess
import pytest

# Limit Polars threads to avoid memory issues during tests
os.environ["POLARS_MAX_THREADS"] = "1"

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))

PANEL_NAME = "ascairn_paper_2025"


@pytest.fixture(scope="session")
def _ascairn_resource_path():
    """Clone ascairn_resource (main branch) if not exists and return base path."""
    resource_path = os.path.join(TESTS_DIR, "ascairn_resource")
    if not os.path.exists(resource_path):
        subprocess.run(
            ["git", "clone", "-b", "main", "https://github.com/friend1ws/ascairn_resource.git"],
            cwd=TESTS_DIR,
            check=True,
        )
    return resource_path


@pytest.fixture(scope="session")
def common_dir(_ascairn_resource_path):
    """Shared BED files (baseline/x_region/cen_region)."""
    return os.path.join(_ascairn_resource_path, "resource", "common")


@pytest.fixture(scope="session")
def panel_dir(_ascairn_resource_path):
    """Panel resource directory (kmer_info/ + hap_info/ + rare_kmer_list.fa)."""
    return os.path.join(_ascairn_resource_path, "resource", "panel", PANEL_NAME)


@pytest.fixture(scope="session")
def output_dir():
    path = os.path.join(TESTS_DIR, "output")
    os.makedirs(path, exist_ok=True)
    return path


@pytest.fixture
def expected_dir():
    return os.path.join(TESTS_DIR, "expected", PANEL_NAME)


@pytest.fixture
def s3_cram():
    return "s3://1000genomes/1000G_2504_high_coverage/additional_698_related/data/ERR3989340/NA12877.final.cram"
