
import os
import pickle #DONT USE PICKLE, CHANGE IT TO NPZ WHEN YOU IMPLEMENT THIS!!!!!
from typing import Optional

"""
DONT USE PICKLE, CHANGE IT TO NPZ WHEN YOU IMPLEMENT THIS!!!!!
"""

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
