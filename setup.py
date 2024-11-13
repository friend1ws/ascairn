from setuptools import setup, find_packages

setup(
    name="ascairn",
    version="0.1",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "ascairn": ["data/*"],
    },
    install_requires=[
        "click",
        # "pysam",
    ],
    extras_require = {
        "s3": ["boto3"],
    },
    entry_points={
        "console_scripts": [
            "ascairn=ascairn.cli:main",
        ],
    },
)
