""":mod:`cpg_islands.models` --- Application models
"""

from __future__ import print_function
import sys
import argparse
from abc import ABCMeta, abstractmethod

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqFeature import SeqFeature, FeatureLocation

from cpg_islands import metadata
from cpg_islands.utils import Event


class MetaApplicationModel(object):
    __metaclass__ = ABCMeta

    started = Event()
    """Fired when the application starts."""

    def __init__(self, sequence_input_model):
        self.sequence_input_model = sequence_input_model

    @abstractmethod
    def run(self, argv=None):
        raise NotImplementedError()

    @abstractmethod
    def load_file(self, file_path):
        self.sequence_input_model.load_file(file_path)


class MetaSequenceInputModel(object):
    __metaclass__ = ABCMeta

    file_loaded = Event()
    error_raised = Event()

    def __init__(self, results_model):
        self.results_model = results_model

    @abstractmethod
    def load_file(self, file_path):
        raise NotImplementedError()

    @abstractmethod
    def annotate_cpg_islands(self, seq, island_size, minimum_gc_ratio):
        raise NotImplementedError()


class MetaResultsModel(object):
    locations_computed = Event()

    @abstractmethod
    def annotate_cpg_islands(self, seq, island_size, minimum_gc_ratio):
        raise NotImplementedError()


class ApplicationModel(MetaApplicationModel):
    def run(self, argv):
        self.started()

        author_strings = []
        for name, email in zip(metadata.authors, metadata.emails):
            author_strings.append('Author: {0} <{1}>'.format(name, email))

        epilog = '''{title} {version}

{authors}
URL: <{url}>
'''.format(
            title=metadata.nice_title,
            version=metadata.version,
            authors='\n'.join(author_strings),
            url=metadata.url)

        arg_parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=metadata.description,
            epilog=epilog)

        args = arg_parser.parse_args(args=argv[1:])

    def load_file(self, file_path):
        self.sequence_input_model.load_file(file_path)


class SequenceInputModel(MetaSequenceInputModel):
    def load_file(self, file_path):
        """Load a sequence file and return the contents.

        :param file_path: the path to the sequence file
        :type file_path: :class:`str`
        :return: the sequence text
        :rtype: :class:`str`
        """
        try:
            seq_record = SeqIO.read(file_path, 'genbank')
        except ValueError as error:
            self.error_raised(str(error))
            return
        self.file_loaded(str(seq_record.seq))

    def annotate_cpg_islands(self, seq, island_size, minimum_gc_ratio):
        self.results_model.annotate_cpg_islands(
            seq, island_size, minimum_gc_ratio)


class ResultsModel(MetaResultsModel):
    def annotate_cpg_islands(self, seq, island_size, minimum_gc_ratio):
        """Create a list of CpG island features in a sequence.

        :param seq: the sequence to annotate
        :type seq: :class:`Bio.Seq.Seq`
        :param island_size: the number of bases which an island may contain
        :type island_size: :class:`int`
        :param minimum_gc_ratio: the ratio of GC to other bases
        :type minimum_gc_ratio: :class:`float`
        :return: list of features
        :rtype: :class:`list` of :class:`Bio.SeqFeature.SeqFeature`
        """
        if island_size <= 0:
            raise ValueError(
                'Invalid island size: {0}'.format(island_size))
        seq_len = len(seq)
        if island_size > seq_len:
            raise ValueError(
                'Island size ({0}) must be less than or '
                'equal to sequence length ({1})'.format(island_size, seq_len))
        if not (0 <= minimum_gc_ratio <= 1):
            raise ValueError('Invalid GC ratio for ratio between '
                             'zero and one: {0}'.format(minimum_gc_ratio))
        minimum_gc_percentage = minimum_gc_ratio * 100
        features = []
        for start_index in xrange(len(seq) - island_size + 1):
            end_index = start_index + island_size
            if GC(seq[start_index:end_index]) >= minimum_gc_percentage:
                feature = SeqFeature(
                    FeatureLocation(start_index, end_index))
                features.append(feature)
        self.locations_computed(features)
