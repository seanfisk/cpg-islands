""":mod:`cpg_islands.models` --- Application models
"""

import argparse
from abc import ABCMeta, abstractmethod

from Bio import SeqIO

from cpg_islands import metadata, algorithms
from cpg_islands.utils import Event


class MetaAppModel(object):
    """Overlying application model interface."""
    __metaclass__ = ABCMeta

    started = Event()
    """Fired when the application starts. Callbacks should look like:

    .. function:: callback()
    """

    islands_computed = Event()
    """Fired when island locations have been computed. Callbacks
    should look like:

    .. function:: callback()

    """

    @abstractmethod
    def register_for_events(self):
        """Register for events fired by other models."""
        raise NotImplementedError()

    @abstractmethod
    def run(self, argv=None):
        """Called when the application starts."""
        raise NotImplementedError()

    @abstractmethod
    def load_file(self, file_path):
        """Direct pass-through to :func:`MetaSeqInputModel.load_file()`.

        :param file_path: the path to the sequence file
        :type file_path: :class:`str`
        """
        raise NotImplementedError()


class MetaSeqInputModel(object):
    __metaclass__ = ABCMeta

    file_loaded = Event()
    """Fired after a file has been loaded into memory. Callbacks
    should look like:

    .. function:: callback(file_contents)

        :param file_contents: contents of the loaded file
        :type file_contents: :class:`str`
    """

    error_raised = Event()
    """Fired when an error occurs. Callbacks should look like:

    .. function:: callback(error_message)

        :param error_message: the error message
        :type error_message: :class:`str`
    """

    island_definition_defaults_set = Event()
    """Fired when default island definitions are set. Callbacks should
    look like:

    .. function:: callback(island_size, min_gc_ratio)

        :param island_size: number of bases in the island
        :type island_size: :class:`int`
        :param min_gc_ratio: minimum ratio of Guanine/Cytosine
        :type min_gc_ratio: :class:`float`
    """

    islands_computed = Event()
    """Fired when island locations have been computed. Callbacks
    should look like:

    .. function:: callback()

    """

    @abstractmethod
    def set_island_definition_defaults():
        """Set the default island definitions of an island size of 200
        and a GC ratio of 60%.
        """
        raise NotImplementedError()

    @abstractmethod
    def load_file(self, file_path):
        """Load a sequence file into memory.

        :param file_path: the path to the sequence file
        :type file_path: :class:`str`
        """
        raise NotImplementedError()

    @abstractmethod
    def compute_islands(self, seq, island_size, min_gc_ratio):
        """Create a list of CpG island features in a sequence.

        :param seq: the sequence to analyze
        :type seq: :class:`Bio.Seq.Seq`
        :param island_size: the number of bases which an island may contain
        :type island_size: :class:`int`
        :param min_gc_ratio: the ratio of GC to other bases
        :type min_gc_ratio: :class:`float`
        :raise: :exc:`ValueError` when parameters are invalid
        """
        raise NotImplementedError()


class MetaResultsModel(object):
    islands_computed = Event()
    """Fired when island locations have been computed. Callbacks
    should look like:

    .. function:: callback(islands)

    :param islands: list of island locations
    :type islands: :class:`list` of :class:`SeqFeature`
    """

    @abstractmethod
    def set_results(self, seq, islands):
        """Set the results of island computation.

        :param seq: the sequence analyzed
        :type seq: :class:`Seq`
        :param islands: list of computed island locations
        :type islands: :class:`list` or :class:`SeqFeature`
        """
        raise NotImplementedError()

    @abstractmethod
    def get_global_seq(self):
        """Returns the entire sequence.

        :return: the full sequence
        :rtype: :class:`str`
        """
        raise NotImplementedError()

    @abstractmethod
    def get_local_seq(self, index):
        """Returns a string representation of a feature based on its index.

        :param index: the index of the feature to be returned
        :type index: :class:`int`
        """
        raise NotImplementedError()


class AppModel(MetaAppModel):
    def __init__(self, seq_input_model):
        """Constructor.

        :param type: :class:`MetaSeqInputModel`
        """
        self.seq_input_model = seq_input_model

    def register_for_events(self):
        self.seq_input_model.islands_computed.append(self.islands_computed)

    def run(self, argv):
        author_strings = []
        for name, email in zip(metadata.authors, metadata.emails):
            author_strings.append('Author: {0} <{1}>'.format(name, email))
        version_str = '{0} {1}'.format(metadata.nice_title, metadata.version)
        epilog = '''{version_str}

{authors}
URL: <{url}>
'''.format(
            title=metadata.nice_title,
            version_str=version_str,
            authors='\n'.join(author_strings),
            url=metadata.url)

        arg_parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=metadata.description,
            epilog=epilog)
        arg_parser.add_argument('--version', '-V',
                                action='version', version=version_str)
        arg_parser.parse_args(args=argv[1:])

        self.seq_input_model.set_island_definition_defaults()
        self.started()

    def load_file(self, file_path):
        self.seq_input_model.load_file(file_path)


class SeqInputModel(MetaSeqInputModel):
    def __init__(self, results_model):
        """Constructor.

        :param type: :class:`MetaResultsModel`
        """
        self.results_model = results_model

    def set_island_definition_defaults(self):
        self.island_definition_defaults_set(200, 0.6)

    def load_file(self, file_path):
        try:
            seq_record = SeqIO.read(file_path, 'genbank')
        except ValueError as error:
            self.error_raised(str(error))
            return
        self.file_loaded(str(seq_record.seq))

    def compute_islands(self, seq, island_size, min_gc_ratio):
        islands = \
            algorithms.registry[0].algorithm(seq, island_size, min_gc_ratio)
        self.results_model.set_results(seq, islands)
        self.islands_computed()


class ResultsModel(MetaResultsModel):
    def __init__(self):
        self._islands = []
        self._seq = ''

    def set_results(self, seq, islands):
        self._seq = seq
        self._islands = islands
        self.islands_computed(islands)

    def get_local_seq(self, index):
        return str(self._islands[index].extract(self._seq))

    def get_global_seq(self):
        return str(self._seq)
