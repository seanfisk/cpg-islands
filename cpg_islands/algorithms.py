""":mod:`cpg_islands.algorithms` --- Algorithms for finding CpG Islands
"""

from __future__ import division
from abc import ABCMeta, abstractmethod, abstractproperty

from Bio.SeqUtils import GC
from Bio.SeqFeature import SeqFeature, FeatureLocation


def _is_gc(base):
    """Test whether a specified base is Guanine or Cytosine.

    :param base: the base to test
    :type base: :class:`str`
    :return: whether the base in G or C
    :rtype: :class:`bool`
    """
    return base == 'G' or base == 'C'


def _count_gc(partial_seq_str):
    """Count the number of Guanine or Cytosine bases in a sequence.

    :param partial_seq_str: the partial sequence
    :type partial_seq_str: :class:`str`
    :return: the GC count
    :rtype: :class:`int`
    """
    gc_count = 0
    for base in partial_seq_str:
        if _is_gc(base):
            gc_count += 1
    return gc_count


def _make_feature(start, end):
    """Create a feature given the start and end indices.

    :param start: inclusive start index
    :type start: :class:`int`
    :param end: exclusive end index
    :type end: :class:`int`
    :return: the created feature
    :rtype: :class:`SeqFeature`
    """
    return SeqFeature(FeatureLocation(start, end))


class MetaAlgorithm(object):
    __metaclass__ = ABCMeta

    @abstractproperty
    def name(self):
        """Return the nice name of this algorithm."""
        raise NotImplementedError()

    @property
    def id(self):
        """Return the identifier for this algorithm."""
        return self.name.lower().replace(' ', '_')

    @abstractmethod
    def algorithm(self, seq_record, island_size, min_gc_ratio):
        """Create a list of CpG island features in a sequence.

        :param seq_record: the sequence record to annotate
        :type seq_record: :class:`SeqRecord`
        :param island_size: the number of bases which an island may contain
        :type island_size: :class:`int`
        :param min_gc_ratio: the ratio of GC to other bases
        :type min_gc_ratio: :class:`float`
        :return: a list of features found in the sequence
        :rtype: :class:`list` of :class:`Bio.SeqFeature.SeqFeature`
        :raise: :exc:`ValueError` when parameters are invalid
        """
        # This is helper code to validate parameters. Call with super.
        if island_size <= 0:
            raise ValueError(
                'Invalid island size: {0}'.format(island_size))
        seq_len = len(seq_record)
        if island_size > seq_len:
            raise ValueError(
                'Island size ({0}) must be less than or '
                'equal to sequence length ({1})'.format(island_size, seq_len))
        if not (0 <= min_gc_ratio <= 1):
            raise ValueError('Invalid GC ratio for ratio between '
                             'zero and one: {0}'.format(min_gc_ratio))


class SlidingWindowBiopythonGCAlgorithm(MetaAlgorithm):
    @property
    def name(self):
        return 'Sliding Window with Biopython GC function'

    def algorithm(self, seq_record, island_size, min_gc_ratio):
        super(SlidingWindowBiopythonGCAlgorithm,
              self).algorithm(seq_record, island_size, min_gc_ratio)
        seq = seq_record.seq
        min_gc_pct = min_gc_ratio * 100
        islands = []
        for start_index in xrange(len(seq) - island_size + 1):
            end_index = start_index + island_size
            if GC(seq[start_index:end_index]) >= min_gc_pct:
                islands.append(_make_feature(start_index, end_index))
        seq_record.features = islands
        return seq_record


class SlidingWindowAccumulator(MetaAlgorithm):
    @property
    def name(self):
        return 'Sliding Window with Count Accumulator'

    def algorithm(self, seq_record, island_size, min_gc_ratio):
        super(SlidingWindowAccumulator,
              self).algorithm(seq_record, island_size, min_gc_ratio)
        seq = seq_record.seq
        seq_len = len(seq)
        start_index = 0
        gc_count = _count_gc(seq[0:island_size])
        islands = []
        while True:
            # check
            # note: end_index is exclusive, like Python's slice
            end_index = start_index + island_size
            if gc_count / island_size >= min_gc_ratio:
                islands.append(_make_feature(start_index, end_index))
            if end_index >= seq_len:
                break
            # update
            if _is_gc(seq[start_index]):
                gc_count -= 1
            if _is_gc(seq[end_index]):
                gc_count += 1
            start_index += 1
        seq_record.features = islands
        return seq_record


# Create instances of each implementation of a MetaAlgorithm, so that
# they are ready to use.
registry = [class_() for class_ in MetaAlgorithm.__subclasses__()]
