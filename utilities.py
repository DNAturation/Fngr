from functools import partial, wraps
import json
import os
import sys

pretty_json = partial(json.dump, indent=4, separators=(', ', ': '))

def cache(func):
    '''Caches the results of func in an external JSON file'''

    @wraps(func)
    def wrapper(*args):

        cache_file = args[0].cache_file
        sequence = args[1]

        # if caching has been turned off
        if not cache_file:
            out = func(*args)

        # if the cache hasn't been created yet
        elif not os.access(cache_file, os.F_OK):

            out = func(*args)
            data = {sequence: out}

            with open(cache_file, 'w') as f:
                pretty_json(data, f)

        # if the cache exists
        else:

            with open(cache_file, 'r+') as f:

                data = json.load(f)

                if sequence in data:
                    out = data[sequence]

                else:

                    f.seek(0)
                    out = func(*args)
                    data[sequence] = out
                    pretty_json(data, f)

        return out
    return wrapper
