#!/usr/bin/python

class RangedDict(dict):
    """
    A dictionary that supports setting items en masse by ranges, but also supports normal keys.

    The core difference between this and any other dict is that by passing a tuple of 2 to 3 numeric values, an
    inclusive range of keys will be set to that value. An example usage is:

    >>> d = RangedDict({
    ...   (1, 5): "foo"
    ... })
    >>> print d[1]  # prints "foo"
    >>> print d[4]  # also prints "foo"
    >>> print d[5]  # still prints "foo" because ranges are inclusive
    >>> d['bar'] = 'baz'
    >>> print d['bar']  # prints 'baz' because this also works like a normal dict

    Do note, ranges are inclusive, so 5 is also set.

    The third, optional, parameter that can be given to a range tuple is a step parameter (analogous to the step
    parameter in xrange), like so: `(1, 5, 2)`, which would set keys 1, 3, and 5 only. For example:

    >>> d[(11, 15, 2)] = "bar"
    >>> print d[13]  # prints "bar"
    >>> print d[14]  # raises KeyError because of step parameter

    NOTE: ALL tuples are strictly interpreted as attempts to set a range tuple. This means that any tuple that does NOT
    conform to the range tuple definition above (e.g., `("foo",)`) will raise a ValueError.


    code was obtained from http://codereview.stackexchange.com/questions/31907/multi-key-dictionary
    written by http://codereview.stackexchange.com/users/67875/hjc1710
    thank your for the code


    """
    def __init__(self, data=None):
        # we add data as a param so you can just wrap a dict literal in the class constructor and it works, instead of
        # having to use kwargs
        if data is None:
            data = {}

        for k, v in data.items():
            if isinstance(k, tuple):
                self._store_tuple(k, v)
            else:
                self[k] = v

    def __setitem__(self, key, value):
        if isinstance(key, tuple):
            self._store_tuple(key, value)
        else:
            # let's go ahead and prevent that infinite recursion, mmmmmkay
            dict.__setitem__(self, key, value)

    def _store_tuple(self, _tuple, value):
        if len(_tuple) > 3 or len(_tuple) < 2:
            # eventually, it would be nice to allow people to store tuples as keys too. Make a new class like: RangeKey
            # to do this
            raise ValueError("Key: {} is invalid! Ranges are described like: (start, stop[, step])")

        step = _tuple[2] if len(_tuple) == 3 else 1
        start = _tuple[0]
        stop = _tuple[1]

        # +1 for inclusivity
        for idx in xrange(start, stop + 1, step):
            dict.__setitem__(self, idx, value)
