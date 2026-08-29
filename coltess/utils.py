"""
Utility helpers for coltess.

Note
----
The checkpoint helpers below are currently unused. Do not rely on
pickle for persistence; switch to NPZ format before using them.
"""

import os
import pickle
from typing import Optional


def save_checkpoint(value: int, filename: str = "checkpoint.pkl") -> None:
    """Save checkpoint"""
    with open(filename, "wb") as f:
        pickle.dump(value, f)


def load_checkpoint(filename: str = "checkpoint.pkl") -> Optional[int]:
    """Load checkpoint"""
    if os.path.exists(filename):
        with open(filename, "rb") as f:
            return pickle.load(f)
    return None
