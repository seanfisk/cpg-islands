""":mod:`cpg_islands.views` --- View interfaces
"""

from cpg_islands.utils import Event


class BaseAppView(object):
    file_load_requested = Event()
    """Called when the user requests to load a file. Callbacks should
    look like:

    .. function:: callback(file_path)

        :param file_path: path the the file to load
        :type file_path: :class:`str`
    """

    def start(self):
        """Start the view."""
        raise NotImplementedError()


class BaseSeqInputView(object):
    submitted = Event()
    """Called when the form is submitted, i.e., submit is clicked by
    the user. Callbacks should look like:

    .. function:: callback(seq_str, island_size_str, minimum_gc_ratio_str)

        :param seq_str: the sequence as a string
        :type seq_str: :class:`str`
        :param island_size_str: number of bases which an island may contain
        :type island_size_str: :class:`str`
        :param minimum_gc_ratio_str: the ratio of GC to other bases
        :type minimum_gc_ratio_str: :class:`str`
    """

    def set_seq(self, seq_str):
        """Set the sequence text.

        :param seq_str: the sequence in string form
        :type seq_str: :class:`str`
        """
        raise NotImplementedError()

    def show_error(self, message):
        """Show the user an error dialog.

        :param message: error message
            :type message: :class:`str`
            """
        raise NotImplementedError()


class BaseResultsView(object):
    def set_locations(self, locations):
        """Set the CpG island locations.

        :param locations: CpG island locations
        :type locations: :class:`list` of :class:`tuple`
        """
        raise NotImplementedError()
