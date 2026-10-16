#!/usr/bin/env python
"""Entry point: see `python run_clone_geno.py --help`."""
import sys

from crispr_clone.cli import main

if __name__ == "__main__":
    sys.exit(main())
