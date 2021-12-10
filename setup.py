#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('./docs/README_summary.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.md') as history_file:
    history = history_file.read()

requirements = [
'cycler==0.10.0',
'gtfparse==1.2.1',
'kiwisolver==1.2.0',
'matplotlib==3.3.1',
'numpy==1.19.2',
'pandas==1.1.3',
'pyBigWig==0.3.18',
'pycairo==1.20.0',
'pyparsing==2.4.7',
'pysam==0.16.0.1',
'python-dateutil==2.8.1',
'pytz==2020.1',
'scipy==1.5.2',
'six==1.15.0'
]

setup_requirements = ['pytest-runner', 'setuptools>=30.3.0', 'wheel', 'setuptools_scm', "Cython"]

test_requirements = ['pytest>=3','pytest-cov','coverage']

setup(
    author="Tomas Drozd, Andreas Tuerk",
    author_email='tomas.drozd@lexogen.com, andreas.tuerk@lexogen.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
	'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'License :: Other/Proprietary License',
    ],
    description="SIRVsuite - a command line tool to QC an RNA-Seq workflow using Lexogen's SIRV spike-in controls",
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
    project_urls={
    	'Documentation': 'https://github.com/Lexogen-Tools/SIRVsuite/blob/master/README.md'
    },
    use_scm_version={
    	"write_to": "SIRVsuite/_version.py",
    	"version_scheme": "post-release",
    },
    zip_safe=False,
)
