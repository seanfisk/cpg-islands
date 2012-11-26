""":mod:`cpg_islands.models` --- Application models
"""

import argparse
from abc import ABCMeta, abstractmethod

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

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
        :param min_obs_exp_cpg_ratio: minimum of observed CpG's to expected
        :type min_obs_exp_cpg_ratio: :class:`float`
    """

    algorithms_loaded = Event()
    """Fired when all algorithms have been loaded.

    .. function:: callback(algorithm_names)

        :param algorithm_names: list of algorithm names
        :types algorithm_names: :class:`list` of :class:`str`
    """

    algorithms_loaded = Event()
    """Fired when all available algorithms have been loaded. Callbacks
    should look like:

    .. function:: callback(algorithm_names)

        :param algorithm_names: list of algorithm names
        :type algorithm_names: :class:`list` of :class:`str`
    """

    islands_computed = Event()
    """Fired when island locations have been computed. Callbacks
    should look like:

    .. function:: callback()

    """

    @abstractmethod
    def set_island_definition_defaults(self):
        """Set the default island definitions.
        """
        raise NotImplementedError()

    @abstractmethod
    def load_algorithms(self):
        """Load all available algorithms."""
        raise NotImplementedError()

    @abstractmethod
    def load_file(self, file_path):
        """Load a sequence file into memory.

        :param file_path: the path to the sequence file
        :type file_path: :class:`str`
        """
        raise NotImplementedError()

    @abstractmethod
    def compute_islands(
            self, seq, island_size, min_gc_ratio,
            min_obs_exp_cpg_ratio, algo_index):
        """Create a list of CpG island features in a sequence.

        :param seq: the sequence to analyze
        :type seq: :class:`Bio.Seq.Seq`
        :param island_size: the number of bases which an island may contain
        :type island_size: :class:`int`
        :param min_gc_ratio: the ratio of GC to other bases
        :type min_gc_ratio: :class:`float`
        :param min_obs_exp_cpg_ratio: minimum observed to expected CpG's
        :type min_obs_exp_cpg_ratio: :class:`float`
        :param algo_index: the index of the algorithm to use
        :type algo_index: :class:`int`
        :raise: :exc:`ValueError` when parameters are invalid
        """
        raise NotImplementedError()


class MetaResultsModel(object):
    islands_computed = Event()
    """Fired when island locations have been computed. Callbacks
    should look like:

    .. function:: callback(seq_record)

        :param seq_record: seq record with features
        :type seq_record: :class:`Bio.SeqRecord.SeqRecord`
    """

    @abstractmethod
    def set_results(self, seq_record):
        """Set the results of island computation.

        :param seq_record: the seq record with features
        :type seq_record: :class:`Bio.SeqRecord.SeqRecord`
        """
        raise NotImplementedError()

    def get_results(self):
        """Return the results of computation.

        :return: the seq record with features
        :rtype: :class:`Bio.SeqRecord.SeqRecord`
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
        self.seq_input_model.load_algorithms()
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
        self.island_definition_defaults_set(200, 0.5, 0.65)

    def load_algorithms(self):
        algorithm_names = [algo.name for algo in algorithms.registry]
        self.algorithms_loaded(algorithm_names)

    def load_file(self, file_path):
        try:
            seq_record = SeqIO.read(file_path, 'genbank')
        except ValueError as error:
            self.error_raised(str(error))
            return
        self.file_loaded(str(seq_record.seq))

    def compute_islands(
            self, seq_record, island_size, min_gc_ratio,
            min_obs_exp_cpg_ratio, algo_index):
        seq_record = \
            algorithms.registry[algo_index].algorithm(
                seq_record, island_size, min_gc_ratio, min_obs_exp_cpg_ratio)
        self.results_model.set_results(seq_record)
        self.islands_computed()


class ResultsModel(MetaResultsModel):
    def __init__(self):
        self._seq_record = SeqRecord(Seq('', IUPAC.unambiguous_dna))

    def set_results(self, seq_record):
        self._seq_record = seq_record
        self.islands_computed(seq_record)

    def get_results(self):
        return self._seq_record
