#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simple functions used by the rest of the module.

Functions
---------
output_json
    Return a string rendition of a dictionary
run
    Run a system command and return stdout, stderr, exit_code
run_cmnd
    Run a system command with run() and return result split by newline
open_zipped
    Open a regular, gzipped, bz2zipped, or open filehandle agnostically
get_length
    Convert kb and mb lengths into integers
chrom_sort_key
    For use with sorted(), allows correct ordering of chromosomes.
"""
import re as _re
import sys as _sys
import bz2 as _bz2
import gzip as _gzip
import subprocess as _sub

import json

__all__ = ["output_json", "run", "run_cmnd", "open_zipped", "get_length",
           "chrom_sort_key"]


def output_json(output):
    """Return JSON formated string."""
    return json.dumps(output, sort_keys=True, indent=2)


def run(command, raise_on_error=False):
    """Run a command with subprocess the way it should be.

    Parameters
    ----------
    command : str
        A command to execute, piping is fine.
    raise_on_error : bool
        Raise a subprocess.CalledProcessError on exit_code != 0

    Returns
    -------
    stdout : str
    stderr : str
    exit_code : int
    """
    pp = _sub.Popen(command, shell=True, universal_newlines=True,
                    stdout=_sub.PIPE, stderr=_sub.PIPE)
    out, err = pp.communicate()
    code = pp.returncode
    if raise_on_error and code != 0:
        _sys.stderr.write(
            '{}\nfailed with:\nCODE: {}\nSTDOUT:\n{}\nSTDERR:\n{}\n'
            .format(command, code, out, err)
        )
        raise _sub.CalledProcessError(
            returncode=code, cmd=command, output=out, stderr=err
        )
    return out, err, code


def run_cmnd(cmnd):
    """Run a command and return the output split into a list by newline."""
    output, _, _ = _run.run(cmnd, raise_on_error=True)
    return output.strip().split('\n')


def open_zipped(infile, mode='r'):
    """Return file handle of file regardless of compressed or not.

    Also returns already opened files unchanged, text mode automatic for
    compatibility with python2.
    """
    # Return already open files
    if hasattr(infile, 'write'):
        return infile
    # Make text mode automatic
    if len(mode) == 1:
        mode = mode + 't'
    if not isinstance(infile, str):
        raise ValueError("I cannot open a filename that isn't a string.")
    if infile.endswith('.gz'):
        return _gzip.open(infile, mode)
    if infile.endswith('.bz2'):
        if hasattr(_bz2, 'open'):
            return _bz2.open(infile, mode)
        return _bz2.BZ2File(infile, mode)
    return open(infile, mode)


def get_length(length_string):
    """Convert a string length like 50kb to an integer length."""
    # Calculate length_string
    if isinstance(length_string, str):
        if length_string.isdigit():
            return int(length_string)
        dist, mod = _re.match(r'(\d+)(\D+)', length_string).groups()
        dist = int(dist)
        mod = mod.lower()
        if mod:
            if not mod in ['kb', 'mb']:
                raise ValueError('Cannot parse {}, must be in kb or mb only'
                                 .format(length_string))
            if mod == 'kb':
                dist = dist * 1000
            elif mod == 'mb':
                dist = dist * 1000000
    elif isinstance(length_string, int):
        return length_string
    else:
        raise ValueError('length_string must be either a string or an integer, is '
                         '{}'.format(type(length_string)))
    return dist


def chrom_sort_key(x):
    """Return an integer for sorting from a chromosome."""
    if x.startswith('chr'):
        x = x[3:]
    if x.upper() == 'X':
        return 100
    elif x.upper() == 'Y':
        return 101
    elif x.upper().startswith('M'):
        return 150
    elif x.isdigit():
        return int(x)
    return x
