from pathlib import Path
import os


def ensure_dir_exists(dir: Path) -> None:
    '''Ensures directory exists by creating the directory recursively if it doesn't exist yet'''
    
    if not dir.exists():
        os.makedirs(dir)