#!/usr/bin/env python

from distutils.core import setup

setup(name='genomon_splicing_mutation',
      version='0.2.0',
      description='Python tools for detecting mutations causing splicing changes',
      author='Yuichi Shiraishi',
      author_email='friend1ws@gamil.com',
      url='https://github.com/Genomon-Project/GenomonSplicingMutation.git',
      package_dir = {'': 'lib'},
      packages=['genomon_splicing_mutation'],
      scripts=['genomon_splicing_mutation'],
      license='GPL-3'
     )

