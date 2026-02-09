"""
Utility components for the nexus package.

Provides geometry functions (Numba-accelerated distance, angle, and wrapping
calculations), display helpers, and performance tracking utilities.
"""

from .geometry import *
from .aesthetics import *
from .performance import *

__all__ = [
    'aesthetics',
    'geometry',
    'performance'
]