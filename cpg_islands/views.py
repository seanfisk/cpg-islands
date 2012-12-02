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

    def show_seq_input(self):
        """Show the sequence input view."""
        raise NotImplementedError()

    def show_results(self):
        """Show the results view."""
        raise NotImplementedError()


class BaseSeqInputView(object):
    submitted = Event()
    """Called when the form is submitted, i.e., submit is clicked by
    the user. Callbacks should look like:

    .. function:: callback(seq_str, island_size_str, min_gc_ratio_str)

        :param seq_str: the sequence as a string
        :type seq_str: :class:`str`
        :param island_size_str: number of bases which an island may contain
        :type island_size_str: :class:`str`
        :param min_gc_ratio_str: the ratio of GC to other bases
        :type min_gc_ratio_str: :class:`str`
    """

    def set_seq(self, seq_str):
        """Set the sequence text.

        :param seq_str: the sequence in string form
        :type seq_str: :class:`str`
        """
        raise NotImplementedError()

    def set_island_size(self, island_size_str):
        """Set the size of the CpG island.

        :param island_size: number of bases in the island as a string
        :type island_size: :class:`str`
        """
        raise NotImplementedError()

    def set_algorithms(self, algorithm_names):
        """Set the list of algorithm names.

        :param algorithm_names: list of algorithm names
        :type algorithm_names: :class:`list` of :class:`str`
        """
        raise NotImplementedError()

    def set_min_gc_ratio(self, min_gc_ratio):
        """Set minimum GC ratio.

        :param min_gc_ratio: the ratio of Guanine/Cytosine as a string
        :type min_gc_ratio: :class:`str`
        """
        raise NotImplementedError()

    def set_min_obs_exp_cpg_ratio(self, min_obs_exp_cpg_ratio):
        """Set the minimum observed/expected CpG ratio.

        :param min_gc_ratio: the ratio of Guanine/Cytosine as a string
        :type min_gc_ratio: :class:`str`
        """
        raise NotImplementedError()

    def show_error(self, message):
        """Show the user an error dialog.

        :param message: error message
        :type message: :class:`str`
        """
        raise NotImplementedError()


class BaseResultsView(object):
    island_selected = Event()

    def set_islands(self, islands):
        """Set the CpG island locations.

        :param islands: CpG island locations
        :type islands: :class:`list` of :class:`tuple`
        """
        raise NotImplementedError()

    def set_algo_name(self, algo_name):
        """Set the name of the algorithm used.

        :param algo_name: the algorithm name
        :type algo_name: :class:`str`
        """
        raise NotImplementedError()

    def set_exec_time(self, exec_time_str):
        """Set the execution time of the algorithm.

        :param exec_time_str: execution time as a string
        :type exec_time: :class:`str`
        """
        raise NotImplementedError()

    def set_global_seq(self, seq_str):
        """Set the global sequence string.

        :param seq_str: DNA sequence of feature
        :type seq_str: :class:`str`
        """
        raise NotImplementedError()

    def highlight_global_seq(self, start, end):
        """Highlight the subsequence within the global sequence.

        :param start: start index of the currently selected island
        :type start: :class:`str`
        :param end: end index of the currently selected island
        :type end: :class:`str`
        """
        raise NotImplementedError()

    def set_start(self, start_str):
        """Set start index of the island.

        :param start_str: start index of island
        :type start_str: :class:`str`
        """
        raise NotImplementedError()

    def set_end(self, end_str):
        """Set end index of the island

        :param end_str: end index of island
        :type end_str: :class:`str`
        """
        raise NotImplementedError()

    def set_length(self, length_str):
        """Set length of the island.

        :param length: length of island
        :type length: :class:`str`
        """
        raise NotImplementedError()

    def set_subseq(self, seq_str):
        """Set the subsequence.

        :param seq_str: DNA sequence of island
        :type seq_str: :class:`str`
        """
        raise NotImplementedError()

    def clear_subseq(self):
        """Clear subsequence field."""
        raise NotImplementedError()

    def set_gc_ratio(self, gc_ratio_str):
        """Set GC ratio of the island.

        :param gc_ratio_str: GC ratio of island
        :type gc_ratio: :class:`str`
        """
        raise NotImplementedError()

    def set_obs_exp_cpg_ratio(self, obs_exp_cpg_ratio_str):
        """Set observed/expected CpG ratio of the island.

        :param obs_exp_cpg_ratio_str: observed/expected ratio
        :type obs_exp_cpg_ratio_str: :class:`str`
        """
        raise NotImplementedError()


class BaseEntrezView(object):
    search_requested = Event()
    query_changed = Event()
    result_selected = Event()
    load_requested = Event()

    def set_suggestion(self, suggestion):
        """Set the suggestions based on spelling.

        :param result: the encoded text
        :type result: :class:`str`
        """
        raise NotImplementedError()

    def set_query_translation(self, query_translation):
        """Set the query text which has been translated to Entrez
        search terms.

        :param query_translation: the translated query
        :type query_translation: :class:`str`
        """
        raise NotImplementedError()

    def set_seq_locus(self, locus, ncbi_url):
        """Set the selected sequence's locus and URL to access on
        NCBI. Basically, identify the sequence on NCBI.

        :param locus: sequence's locus
        :type locus: :class:`str`
        :param ncbi_url: URL to sequence on NCBI
        :type ncbi_url: :class:`str`
        """
        raise NotImplementedError()

    def set_seq_desc(self, desc):
        """Set the selected sequence's description.

        :param desc: sequence's description
        :type desc: :class:`str`
        """
        raise NotImplementedError()

    def set_selected_seq(self, seq_str):
        """Set the selected sequence.

        :param seq_str: the sequence data
        :type seq_str: :class:`str`
        """
        raise NotImplementedError()

    def set_result(self, results):
        """Set the list of sequence ids.

        :param result: list of sequence ids
        :type result: :class:`list`
        """
        raise NotImplementedError()
