#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.md') as history_file:
    history = history_file.read()

requirements = [
'numpy==1.19.2',
'pycairo==1.20.0',
'scipy==1.5.2',
'pysam==0.16.0.1',
'pandas==1.1.3',
'matplotlib==3.0.0',
'gtfparse==1.2.1',
'pyBigWig==0.3.17'
]

setup_requirements = ['pytest-runner', 'setuptools>=30.3.0', 'wheel', 'setuptools_scm', "Cython"]

test_requirements = ['pytest>=3','pytest-cov','coverage']

setup(
    author="Tomas Drozd, Andreas Tuerk",
    author_email='tomas.drozd@lexogen.com, andreas.tuerk@lexogen.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
	'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="SIRVsuite - a tool for spike-in analysis in NGS data",
    entry_points={
        'console_scripts': [
            'SIRVsuite=SIRVsuite.__main__:main',
        ],
    },
    install_requires=requirements,
    long_description=readme + '\n\n' + history,
    long_description_content_type='text/markdown',
    include_package_data=True,
    keywords='SIRVsuite',
    name='SIRVsuite',
    packages=find_packages(include=['SIRVsuite', 'SIRVsuite.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/TomasDrozd/SIRVsuite',
    version="0.0.dev4",
    zip_safe=False,
)
