#!/usr/bin/env python

from distutils.core import setup

setup(name='savnet',
      version='0.2.1',
      description='Python tools for detecting mutations causing splicing changes',
      author='Yuichi Shiraishi',
      author_email='friend1ws@gamil.com',
      url='https://github.com/Genomon-Project/GenomonSplicingMutation.git',
      package_dir = {'': 'lib'},
      packages=['savnet'],
      scripts=['savnet'],
      license='GPL-3'
     )

