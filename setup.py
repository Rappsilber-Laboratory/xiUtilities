#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    "pandas>=2.2.2",
    "numpy>=1.26.4",
    "setuptools>=71.1.0",
    "numpy-indexed>=0.3.7",
    "tqdm>=4.66.0",
    "pyteomics>=4.6",
    "seaborn>=0.12.2",
    "pytest>=8.3.1",
    "pydocstyle>=6.3.0",
    "pytest-cov>=4.0.0",
    "pytest-flake8>=1.0.6",
    "pytest-pydocstyle>=2.3.2",
    "py>=1.11.0",
    "flake8==4.0.1",
    "networkx>=3.1",
    "multiprocess>=0.70.16",
    "deepmerge~=1.1.0",
    "fastparquet>=2024.5.0",
    "pyarrow>=17.0.0",
    "scipy>=1.12.0",
    "matplotlib",
]

test_requirements = ['pytest>=3', ]

setup(
    author="Falk Boudewijn Schimweg",
    author_email='git@falk.schimweg.de',
    python_requires='>=3.7',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    description="A collection of corss-link mass spectrometry tools.",
    install_requires=requirements,
    license="GNU General Public License v3",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='xi_utils',
    name='xi_utils',
    packages=find_packages(include=['xi_utils', 'xi_utils.*'], exclude=["tests"]),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/Rappsilber-Laboratory/xiUtils',
    zip_safe=False,
)
