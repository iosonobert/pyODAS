#!/usr/bin/env python

from setuptools import setup, find_packages
from setuptools.dist import Distribution

class BinaryDistribution(Distribution):
    def is_pure(self):
        return False

		
setup(name='pyODAS',
      version='1.0.0',
      url='https://github.com/iosonobert/pyODAS/',
      description='partial translation of RSI ODAS library',
      author='Andrew Zulberti',
      author_email='andrew.zulberti@gmail.com',
      packages=find_packages(),
      install_requires=['numpy','matplotlib', 'xarray'],
      license='unlicensed to all but author',
      include_package_data=True,
      distclass=BinaryDistribution,
    )
