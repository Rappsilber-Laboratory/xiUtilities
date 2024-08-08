from importlib.metadata import version

from pkg_resources import DistributionNotFound

"""Top-level package for xiRescore."""

__author__ = """Falk Boudewijn Schimweg"""
__email__ = 'git@falk.schimweg.de'
try:
    __version__ = version('xiutilities')
except DistributionNotFound:
    # package is not installed
    __version__ = '0.0.0'
