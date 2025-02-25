#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    "pandas>=1.3.5",
    "numpy>=1.7.2",
    "setuptools",
    "numpy-indexed>=0.3.7",
    "tqdm>=4.66.0",
    "pyteomics>=4.6",
    "seaborn>=0.11.2",
    "pytest",
    "pydocstyle",
    "pytest-cov",
    "pytest-flake8",
    "pytest-pydocstyle",
    "flake8",
    "networkx>=2.6.3",
    "multiprocess>=0.70.15",
    "deepmerge>=1.1.0",
    "scipy>=1.0.1",
    "matplotlib>=3.5.3",
    "setuptools_scm>=7.1",
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
    keywords='xiutilities',
    name='xiutilities',
    packages=find_packages(include=['xiutilities', 'xiutilities.*'], exclude=["tests"]),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/Rappsilber-Laboratory/xiUtils',
    zip_safe=False,
)
