""":mod:`cpg_islands.views` --- View interfaces
"""


class BaseApplicationView(object):
    def start(self):
        """Start the view."""
        raise NotImplementedError()
