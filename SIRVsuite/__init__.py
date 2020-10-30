try:
    from SIRVsuite._version import version
    __version__ = version
except ImportError:
    __version__ = "not-installed"
