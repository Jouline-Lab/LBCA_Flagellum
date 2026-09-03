# -*- coding: utf-8 -*-
"""
Drop-in replacement for the `pastml` CLI command that works around a bug
in pastml.acr: it calls sys.setrecursionlimit(recursion_limit) without
ever importing sys in that module, so any run using --recursion_limit
crashes with NameError: name 'sys' is not defined (observed with
pastml installed from PyPI, Python 3.9).

Rather than editing the installed package file, this injects `sys` into
pastml.acr's own namespace before calling its normal CLI entry point, so
the rest of the module works exactly as installed. All command-line
arguments are forwarded through untouched -- this is a straight
substitute for the `pastml` command, not a different interface.

Usage: python pastml_wrapper.py <same arguments you would pass to `pastml`>
"""

import sys
import pastml.acr as acr

acr.sys = sys  # patches the missing "import sys" without touching the installed file

if __name__ == "__main__":
    sys.exit(acr.main())
