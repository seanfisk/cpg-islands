""":mod:`cpg_islands.views` --- View interfaces
"""


class BaseApplicationView(object):
    def start(self):
        """Start the view."""
        raise NotImplementedError()

    def set_locations(self, locations):
        """Set the CpG island locations.

        :param locations: CpG island locations
        :type locations: :class:`list` of :class:`tuple`
        """
        raise NotImplementedError()
