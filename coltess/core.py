#!/usr/bin/env python3
"""
Core data structures for coltess.
"""

from dataclasses import dataclass, field
from typing import Optional
import numpy as np


@dataclass
class StarData:
    """
    This class stores the star's identifying information (name, coordinates,
    Gaia ID) and, once analysis is performed, its light-curve data (times
    and fluxes).

    The object is intentionally mutable: analysis functions are expected to
    populate photometry fields after creation.
    """

    # Identity / catalog info
    name: str
    ra: float
    dec: float
    gaia_id: Optional[str] = None

    # Photometry results (optional)
    times: Optional[np.ndarray] = field(default=None, repr=False)
    fluxes: Optional[np.ndarray] = field(default=None, repr=False)

    # Metadata
    sector: Optional[int] = None

    def has_photometry(self) -> bool:
        return self.times is not None and self.fluxes is not None
