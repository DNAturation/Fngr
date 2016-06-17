from functools import partial, wraps
import json
import os

pretty_json = partial(json.dump, indent=4, separators=(', ', ': '))

def cache(func):
    '''Caches the results of func in an external JSON file'''

    @wraps(func)
    def wrapper(*args):

        # args[0] is 'self'
        cache_obj = args[0].cache_obj
        sequence = args[1]

        # if caching has been turned off
        if cache_obj is None:
            out = func(*args)

        # the sequence is in the cache
        elif sequence in cache_obj:

            out = cache_obj[sequence]

        # the sequence is not in the cache
        else:

            out = func(*args)
            cache_obj[sequence] = out

        return out, cache_obj
    return wrapper

def load_cache(path):

    if path is not None:

        if not os.access(path, os.F_OK):
            cache_obj = {}

        else:

            with open(path) as f:
                cache_obj = json.load(f)

        return cache_obj

def write_cache(obj, path):

    if path is not None:
        with open(path, 'w') as f:
            pretty_json(obj, f)
