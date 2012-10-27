""":mod:`cpg_islands.views` --- View interfaces
"""

from cpg_islands.utils import Event


class BaseApplicationView(object):
    submitted = Event()
    """Called when the form is submitted, i.e., submit is clicked by
    the user."""

    sequence_changed = Event()
    """Called when the sequence text is changed."""

    def start(self):
        """Start the view."""
        raise NotImplementedError()

    def set_sequence(self, sequence_str):
        """Set the sequence text.

        :param sequence_str: the sequence
        :type sequence_str: :class:`str`
        """
        raise NotImplementedError()

    def set_locations(self, locations):
        """Set the CpG island locations.

        :param locations: CpG island locations
        :type locations: :class:`list` of :class:`tuple`
        """
        raise NotImplementedError()

    def _submit_clicked(self):
        """Called when user clicks submit."""
        raise NotImplementedError()
