# -*- coding: utf-8 -*-
#                                                     #
#  __author__ = Adarsh Kalikadien                     #
#  __institution__ = TU Delft                         #
#  __contact__ = a.v.kalikadien@tudelft.nl            #

from setuptools import setup, find_packages

setup(
    name='chemspax',
    version='0.1.0',
    packages=find_packages(include=['chemspax', 'chemspax.*']),
    url='github.com/epics-group/chemspax',
    license='MIT',
    author='Adarsh Kalikadien',
    author_email='a.v.kalikadien@tudelft.nl',
    description='A Python tool for local chemical space exploration of any structure based on their 3D geometry',

    # use entry point to create a command that will run the main function (if needed)
    # entry_points={
    #     'console_scripts': ['main = chemspax.main:main']
    # },
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],

    # ToDo: add package data and change paths to pkg_data_dirs in attach_substitutions.py
    # source: https://kiwidamien.github.io/making-a-python-package-vi-including-data-files.html
    # include_package_data=True,
    package_data={'chemspax': ['substituents_xyz', 'skeletons']},

    long_description=open('docs/README.md').read(),


)

