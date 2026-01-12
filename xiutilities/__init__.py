import setuptools_scm

from pkg_resources import DistributionNotFound

"""Top-level package for xiRescore."""

def get_version():
    # Try to get version from version file
    try:
        from ._version import __version__
        return __version__
    except:
        pass
    # Try to get version from workdir
    try:
        return setuptools_scm.get_version()
    except:
        pass
    # No version found
    return "unkown"

__author__ = """Falk Boudewijn Schimweg"""
__email__ = 'git@falk.schimweg.de'
__version__ = get_version()
