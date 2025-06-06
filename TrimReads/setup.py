#!/usr/bin/env python
from setuptools import setup, find_packages
import os
import re

# Get package version
VERSION = "1.0.0"
#def get_version():
#    version_file = os.path.join(os.path.dirname(__file__), "src", "trimreads", "__init__.py")
#    with open(version_file, "r") as f:
#        version_content = f.read()
#    version_match = re.search(r"^__version__\s*=\s*['\"]([^'\"]*)['\"]", version_content, re.M)
#    if version_match:
#        return version_match.group(1)
#    raise RuntimeError("Unable to find version string.")

# Get long description from README
def get_long_description():
    with open("README.md", "r", encoding="utf-8") as f:
        return f.read()

# Package requirements
requirements = [
    'biopython>=1.79',
    'numpy>=1.21',
    'tqdm>=4.62',  # For progress bars
]

# Test requirements
test_requirements = [
    'pytest>=6.2',
    'pytest-cov>=3.0',
]

# Setup configuration
setup(
    name="TrimReads",
    version=VERSION,
    #version=get_version(),
    author="Bioinformatics Team",
    author_email="bioinfo@example.com",
    description="A high-throughput sequencing data quality trimming tool",
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/TrimReads",
    project_urls={
        "Bug Tracker": "https://github.com/yourusername/TrimReads/issues",
        "Documentation": "https://trimreads.readthedocs.io",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "Development Status :: 4 - Beta",
    ],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    include_package_data=True,
    package_data={
        "trimreads": [
            "tests/test_data/*.fastq",
            "data/*.fastq",
            "docs/*.md",
        ],
    },
    install_requires=requirements,
    tests_require=test_requirements,
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "trimreads = trimreads.trimreads:main",
            "fastq_utils = trimreads.fastq_parser:main",
        ],
    },
    extras_require={
        "full": ["matplotlib", "seaborn"],  # For advanced visualization in demos
        "dev": test_requirements + [
            "flake8",
            "black",
            "twine",
            "sphinx",
            "sphinx_rtd_theme",
        ],
    },
    options={
        "bdist_wheel": {
            "universal": True,
        },
    },
    license="MIT",
    keywords="bioinformatics sequencing fastq quality-trimming",
)