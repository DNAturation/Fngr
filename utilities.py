from functools import partial, wraps
import json
import os

pretty_json = partial(json.dump, indent=4, separators=(', ', ': '))

def cache(func, cache_file):
    '''Caches the results of func in an external JSON file'''

    @wraps(func)
    def wrapper(sequence):

        # if caching has been turned off
        if not cache_file:
            out = func(sequence)

        # if the cache hasn't been created yet
        elif not os.access(cache_file, os.F_OK):

            out = func(sequence)
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

                    out = func(sequence)
                    data[sequence] = out

                    pretty_json(data, f)

        return out
