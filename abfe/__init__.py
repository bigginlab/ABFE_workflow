import os

root_path = os.path.dirname(__file__)

from abfe.calculate_abfe import calculate_abfe
from abfe._version import __version__, __version_tuple__
__version__ = __version__
__version_tuple__ = __version_tuple__