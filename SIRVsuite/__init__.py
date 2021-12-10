try: # pragma: no cover
    from SIRVsuite._version import version # pragma: no cover
    __version__ = version # pragma: no cover
except ImportError: # pragma: no cover 
    __version__ = "not-installed" # pragma: no cover
