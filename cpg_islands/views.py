""":mod:`cpg_islands.views` --- View interfaces
"""

from cpg_islands.utils import Event


class BaseApplicationView(object):
    file_load_requested = Event()
    """Called when the user requests to load a file."""

    def __init__(self, sequence_input_view, results_view):
        """Initialize the main application view with docked
        SeqenceInputView and ResultsView.

        :param sequence_input_view: the input view
        :type sequence_input_view: :class:`BaseSequenceInputView`
        :param results_view: the results view
        :type results_view: :class:`BaseResultsView`
        """
        raise NotImplementedError()

    def start(self):
        """Start the view."""
        raise NotImplementedError()


class BaseSequenceInputView(object):
    submitted = Event()
    """Called when the form is submitted, i.e., submit is clicked by
    the user.
    """

    def set_sequence(self, sequence_str):
        """Set the sequence text.

        :param sequence_str: the sequence in string form
        :type sequence_str: :class:`str`
        """
        raise NotImplementedError()


class BaseResultsView(object):
    def set_locations(self, locations):
        """Set the CpG island locations.

        :param locations: CpG island locations
        :type locations: :class:`list` of :class:`tuple`
        """
        raise NotImplementedError()

    def show_error(self, message):
        """Show the user an error dialog.

        :param message: error message
        :type message: :class:`str`
        """
        raise NotImplementedError()
