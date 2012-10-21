#!/usr/bin/env python

import os
from cpg_islands import metadata

# auto-install and download distribute
import distribute_setup
distribute_setup.use_setuptools()

from setuptools import setup, find_packages


# credit: <http://packages.python.org/an_example_pypi_project/setuptools.html>
# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

install_requirements = []

# see here for more options:
# <http://packages.python.org/distribute/setuptools.html>
setup(name=metadata.title,
      version=metadata.version,
      author=metadata.authors[0],
      author_email=metadata.emails[0],
      maintainer=metadata.authors[0],
      maintainer_email=metadata.emails[0],
      url=metadata.url,
      description=metadata.description,
      long_description=read('README.rst'),
      download_url=metadata.url,
      # find a list of classifiers here:
      # <http://pypi.python.org/pypi?%3Aaction=list_classifiers>
      classifiers=[
          'Development Status :: 2 - Pre-Alpha',
          'Environment :: Console',
          'Intended Audience :: Developers',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Healthcare Industry',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Information Analysis',
          'Topic :: Scientific/Engineering :: Medical Science Apps.',
          'Topic :: Scientific/Engineering :: Visualization',
      ],
      packages=find_packages(),
      install_requires=install_requirements,
      zip_safe=False,  # don't use eggs
      entry_points={
          'console_scripts': [
              'cpg_islands = cpg_islands.cli.main:main'
          ],
          # if you have a gui, use this
          # 'gui_scripts': [
          #     'my_module_gui = my_module.gui:main'
          # ]
      }
      )
