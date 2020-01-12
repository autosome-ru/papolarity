import gzip
import sys
from .nullcontext import nullcontext

def choose_open_function(filename, force_gzip=None):
    '''
    If `force_gzip` is True or False, use corresponding open function.
    If `force_gzip` is None, it's guessed based on filename extension
    '''
    if (force_gzip is True) or ((force_gzip is None) and filename.endswith('.gz')):
        return gzip.open
    elif (force_gzip is False) or ((force_gzip is None) and not filename.endswith('.gz')):
        return open
    else:
        raise ValueError("`force_gzip` should be one of True/False/None")

def open_for_write(filename, force_gzip=None, mode='wt', **kwargs):
    if filename and (filename != '-'):
        open_func = choose_open_function(filename=filename, force_gzip=force_gzip)
        return open_func(filename, mode, **kwargs)
    else:
        return nullcontext(sys.stdout)

def open_for_read(filename, force_gzip=None, mode='rt', **kwargs):
    if filename and (filename != '-'):
        open_func = choose_open_function(filename=filename, force_gzip=force_gzip)
        return open_func(filename, mode, **kwargs)
    else:
        return nullcontext(sys.stdin)
