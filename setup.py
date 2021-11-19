#! /usr/bin/env python3

import os
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

PACKAGE = 'mes'
VMAJOR, VMINOR, VMICRO = 0, 0, 2
VERSION = '{}.{}.{}'.format(VMAJOR, VMINOR, VMICRO)

def readme():
    readme_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'README.md'
    )
    with open(readme_file) as f:
        return f.read()

long_description=readme()

setup(
    name = PACKAGE,
    version = VERSION,
    description = "The mes software.",
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    classifiers = [
        'Development Status :: 1 - Planning',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)'
    ],
    url = 'https://github.com/xiaolab/mes',
    author = 'Zhiheng Liu',
    author_email = 'wolfsonliu@live.com',
    license = 'GPL',
    packages = find_packages(),
    install_requires = [
        'numpy>=1.10',
        'pandas>=0.16',
        'matplotlib>=2.0.0',
        'pysam'
    ],
    python_requires = '>=3.7',
    scripts = [
        'bin/mes'
    ],
    package_dir={'mes': 'mes'},
    package_data = {
        # If any package contains *.txt or *.rst files, include them:
        'mes.test': ['*.test_data'],
        'mes.template': ['*.html', '*.svg'],
    },
    zip_safe=False
)
