""":mod:`cpg_islands.utils` --- Miscellaneous utilities
"""


# credit: <http://stackoverflow.com/a/2022629>
class Event(list):
    """Event subscription.

    A list of callable objects. Calling an instance of this will cause a
    call to each item in the list in ascending order by index.

    Example Usage:
    >>> def f(x):
    ...     print 'f({0})'.format(x)
    >>> def g(x):
    ...     print 'g({0})'.format(x)
    >>> e = Event()
    >>> e()
    >>> e.append(f)
    >>> e(123)
    f(123)
    >>> e.remove(f)
    >>> e()
    >>> e += (f, g)
    >>> e(10)
    f(10)
    g(10)
    >>> del e[0]
    >>> e(2)
    g(2)
    """
    def __call__(self, *args, **kwargs):
        for f in self:
            f(*args, **kwargs)

    def __repr__(self):
        return 'Event({0})'.format(list.__repr__(self))
