#! /usr/bin/env python3

import os
try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages

PACKAGE = 'lqc'
VMAJOR, VMINOR, VMICRO = 0, 0, 5
VERSION = '{}.{}.{}'.format(VMAJOR, VMINOR, VMICRO)


def readme():
    readme_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'README.md'
    )
    with open(readme_file) as f:
        return f.read()


long_description = readme()


setup(
    name = PACKAGE,
    version = VERSION,
    description = 'The Long-read RNA-seq quality control software.',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    classifiers = [
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    url = 'https://github.com/xiaolab/LQC',
    author = 'Zhiheng Liu',
    author_email = 'wolfsonliu@live.com',
    license = 'GPL',
    packages = find_packages(),
    install_requires = [
        'numpy',
        'pandas',
        'matplotlib',
        'pysam'
    ],
    python_requires = '>=3.7',
    scripts = [
        'bin/lqc'
    ],
    package_dir={'lqc': 'lqc'},
    package_data = {
        # If any package contains *.txt or *.rst files, include them:
        'lqc.test': ['*.test_data'],
        'lqc.template': ['*.html', '*.svg'],
    },
    zip_safe=False
)
