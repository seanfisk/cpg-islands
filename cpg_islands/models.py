""":mod:`cpg_islands.models` --- Application models
"""

import argparse
from abc import ABCMeta, abstractmethod
import timeit

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from cpg_islands import metadata, algorithms
from cpg_islands.utils import Event


class IslandInfo(object):
    """Container class for island information."""
    def __init__(self, start, end, length, subseq, gc_ratio,
                 obs_exp_cpg_ratio):
        """Constructor.

        :param start: island start index
        :type start: :class:`int`
        :param end: island end index
        :type end: :class:`int`
        :param length: length of the island
        :type length: :class:`int`
        :param subseq: the bases that make up the island
        :type subseq: :class:`str`
        :param gc_ratio: island GC ratio
        :type gc_ratio: :class:`float`
        :param obs_exp_cpg_ratio: island observed/expected CpG ratio
        :type obs_exp_cpg_ratio: :class:`float`
        """
        self.start = start
        self.end = end
        self.length = length
        self.subseq = subseq
        self.gc_ratio = gc_ratio
        self.obs_exp_cpg_ratio = obs_exp_cpg_ratio

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not self == other


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

    .. function:: callback(global_seq, feature_tuples, algo_name, exec_time)

        :param global_seq: the full sequence
        :type global_seq: :class:`str`
        :param feature_tuples: tuples containing locations of the features
        :type feature_tuples: :class:`list` of :class:`tuple`
        :param algo_name: the name of the algorithm used
        :type algo_name: :class:`str`
        :param exec_time: algorithm's execution duration
        :type exec_time: :class:`float`
    """

    @abstractmethod
    def set_results(self, results, algo_name, exec_time):
        """Set the results of island computation.

        :param results: algorithm results
        :type results: :class:`AlgoResults`
        :param algo_name: the name of the algorithm used
        :type algo_name: :class:`str`
        :param exec_time: the algorithm's execution time
        :type exec_time: :class:`float`
        """
        raise NotImplementedError()

    @abstractmethod
    def get_island_info(self, island_index):
        """Return the island information corresponding to a certain index.

        :return: island information
        :rtype: :class:`IslandInfo`
        """
        raise NotImplementedError()


class MetaEntrezModel(object):
    started = Event()
    """Responsible for initializing the Entrez email"""
    @abstractmethod
    def search(self, text):
        """Search Entrez database.

        :param text: the text to search
        :type text: :class:`str`
        :return: the results
        :rtype: :class:`dict`
        """
        raise NotImplementedError()

    @abstractmethod
    def suggest(self, text):
        """Suggested Entrez spelling.

        :param text: unchecked text from input
        :type text: :class:`str`
        :return: a suggested query
        :rtype: :class:`list` of `dict`
        """
        raise NotImplementedError()

    @abstractmethod
    def get_seq(self, id):
        """Pull sequence based on id.

        :param id: the id of the sequence
        :type id: :class:`int`
        :return: the sequence
        :rtype: :class:`Seq`
        """
        raise NotImplementedError()

    @abstractmethod
    def load_seq(self, seq):
        """Load sequence into SeqInputView.

        :param seq: the id of the sequence
        :type seq: :class:`str`
        """
        raise NotImplementedError()


class AppModel(MetaAppModel):
    def __init__(self, seq_input_model, entrez_model):
        """Constructor.

        :param type: :class:`MetaSeqInputModel`
        """
        self.entrez_model = entrez_model
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
        self.island_definition_defaults_set(200, 0.5, 0.6)

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
        start = timeit.default_timer()
        algo = algorithms.registry[algo_index]

        seq_record = algo.algorithm(
            seq_record, island_size, min_gc_ratio, min_obs_exp_cpg_ratio)

        end = timeit.default_timer()

        self.results_model.set_results(seq_record, algo.name, end - start)
        self.islands_computed()


class ResultsModel(MetaResultsModel):
    def __init__(self):
        self.island_information = algorithms.AlgoResults(
            SeqRecord(Seq('', IUPAC.unambiguous_dna)),
            [])

    def set_results(self, results, algo_name, exec_time):
        self.results = results
        island_tuples = [(feature.location.start.position,
                          feature.location.end.position) for feature in
                         self.results.seq_record.features]
        self.islands_computed(str(results.seq_record.seq),
                              island_tuples, algo_name, exec_time)

    def get_island_info(self, island_index):
        feature = self.results.seq_record.features[island_index]
        metadata = self.results.island_metadata_list[island_index]

        start = feature.location.start.position
        end = feature.location.end.position
        length = end - start
        subseq = str(feature.extract(self.results.seq_record.seq))
        gc_ratio = metadata.gc_ratio
        obs_exp_cpg_ratio = metadata.obs_exp_cpg_ratio

        return IslandInfo(start, end, length, subseq, gc_ratio,
                          obs_exp_cpg_ratio)


class EntrezModel(MetaEntrezModel):
    def __init__(self, seq_input_model):
        Entrez.email = 'gray.gwizdz@gmail.com'
        self.seq_input_model = seq_input_model
        self.results = {'IdList': [], 'QueryTranslation': ''}

    def search(self, text):
        handle = Entrez.esearch(db='nucleotide', term=text)
        results = Entrez.read(handle)
        self.results = {'IdList': results['IdList'],
                        'QueryTranslation': results['QueryTranslation']}
        return self.results

    def suggest(self, text):
        handle = Entrez.espell(db='pubmed', term=text)
        return Entrez.read(handle)

    def get_seq(self, id):
        handle = Entrez.efetch(db='nucleotide', id=id,
                               rettype='gb', retmode='text')
        self.record = SeqIO.read(handle, 'genbank')
        handle.close()
        return self.record

    def load_seq(self, seq):
        self.seq_input_model.file_loaded(seq)
