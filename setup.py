#!/usr/bin/env python

from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name = 'savnet',
    version = '0.3.0b1',
    description='Python tools for detecting mutations causing splicing changes',
    url = 'https://github.com/friend1ws/SAVNet',
    author = 'Yuichi Shiraishi',
    author_email = 'friend1ws@gamil.com',
    license = 'GPLv3',

    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: Unix',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],

    packages = find_packages(exclude = ['tests']),
    package_data={'sv_utils': ['data/*']},

    install_requires = ["annot_utils", "pysam", "junc_utils", "intron_retention_utils", "chimera_utils"]
    install_requires = [],
    entry_points = {'console_scripts': ['savnet = savnet:main']}

)


