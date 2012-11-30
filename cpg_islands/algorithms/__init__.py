""":mod:`cpg_islands.algorithms` --- Algorithms for finding CpG Islands
"""

from __future__ import division
from abc import ABCMeta, abstractmethod, abstractproperty

from Bio.SeqFeature import SeqFeature, FeatureLocation

# from cpg_islands.algorithms import sliding_window_cython


class IslandMetadata(object):
    """Container class for island metadata."""
    def __init__(self, gc_ratio, obs_exp_cpg_ratio):
        """Constructor.

        :param gc_ratio: this island's GC ratio
        :type gc_ratio: :class:`float`
        :param obs_exp_cpg_ratio: this islands's observed/expected CpG ratio
        :type obs_exp_cpg_ratio: :class:`float`
        """
        self.gc_ratio = gc_ratio
        self.obs_exp_cpg_ratio = obs_exp_cpg_ratio

    def __eq__(self, other):
        return (self.gc_ratio == other.gc_ratio and
                self.obs_exp_cpg_ratio == other.obs_exp_cpg_ratio)

    def __ne__(self, other):
        return not self == other


class AlgoResults(object):
    """Container class for algorithm results."""
    def __init__(self, seq_record, island_metadata_list):
        """Constructor.

        :param seq_record: seq record with feature added
        :type seq_record: :class:`SeqRecord`
        :param island_metadata_list: list of metadata for each island
        :type island_metadata_list: :class:``
        """
        self.seq_record = seq_record
        self.island_metadata_list = island_metadata_list

    def extract_features(self):
        return [str(feature.extract(self.seq_record.seq)) for feature
                in self.seq_record.features]

    def __eq__(self, other):
        if str(self.seq_record) != str(other.seq_record):
            return False
        if self.extract_features() != other.extract_features():
            return False
        return self.island_metadata_list == other.island_metadata_list

    def __ne__(self, other):
        return not self == other


def _compute_counts(subseq):
    """Count the number of Guanine bases, Cytosine bases, and CpG in a
    subsequence.

    :param subseq: the partial sequence
    :type subseq: :class:`str`
    :return: a tuple of ``(c_count, g_count, cpg_count)``
    :rtype: :class:`tuple`
    """
    c_count = 0
    g_count = 0
    cpg_count = 0
    seq_len = len(subseq)
    i = 0
    while True:
        base = subseq[i]
        if base == 'C':
            c_count += 1
        elif base == 'G':
            g_count += 1
        if i >= seq_len - 1:
            break
        if subseq[i:i + 2] == 'CG':
            cpg_count += 1
        i += 1
    return (c_count, g_count, cpg_count)


def _compute_ratios(c_count, g_count, cpg_count, subseq_len):
    """Compute ratios given counts.

    :param c_count: number of C's in the subsequence
    :type c_count: :class:`int`
    :param g_count: number of G's in the subsequence
    :type g_count: :class:`int`
    :param cpg_count: number of CpG's in the subsequence
    :type cpg_count: :class:`int`
    :param subseq_len: length of the subsequence
    :type subseq_len: :class:`int`
    :return: GC ratio and observed/expected CpG ratio
    :rtype: :class:`tuple` of (gc_ratio, obs_exp_cpg_ratio)
    """
    gc_ratio = (c_count + g_count) / subseq_len
    obs_exp_cpg_ratio = \
        cpg_count / ((c_count * g_count) / subseq_len)
    return (gc_ratio, obs_exp_cpg_ratio)


def _is_island(gc_ratio, obs_exp_cpg_ratio, min_gc_ratio,
               min_obs_exp_gc_ratio):
    """Determine whether these ratios constitute an island.

    :param gc_ratio: the ratio of GC to other bases
    :type gc_ratio: :class:`float`
    :param obs_exp_cpg_ratio: observed-to-expected CpG ratio
    :type obs_exp_cpg_ratio: :class:`float`
    :param min_gc_ratio: the minimum ratio of GC to other bases
    :type min_gc_ratio: :class:`float`
    :param min_obs_exp_cpg_ratio: minimum observed-to-expected CpG ratio
    :type min_obs_exp_cpg_ratio: :class:`float`
    :return: whether this subsequence is an island
    :rtype: :class:`bool`
    """
    return (gc_ratio >= min_gc_ratio and
            obs_exp_cpg_ratio >= min_obs_exp_gc_ratio)


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
        return (self.name.lower().
                replace(' ', '_').
                replace('(', '').
                replace(')', ''))

    @abstractmethod
    def algorithm(
            self, seq_record, island_size, min_gc_ratio,
            min_obs_exp_cpg_ratio):
        """Create a list of CpG island features in a sequence.

        :param seq_record: the sequence record to annotate
        :type seq_record: :class:`SeqRecord`
        :param island_size: the number of bases which an island may contain
        :type island_size: :class:`int`
        :param min_gc_ratio: the ratio of GC to other bases
        :type min_gc_ratio: :class:`float`
        :param min_obs_exp_cpg_ratio: minimum observed-to-expected CpG ratio
        :type min_obs_exp_cpg_ratio: :class:`float`
        :return: container class of algorithm results
        :rtype: :class:`AlgoResults`
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
            raise ValueError(
                'Invalid GC ratio for ratio between '
                'zero and one: {0}'.format(min_gc_ratio))
        if not (0 <= min_obs_exp_cpg_ratio):
            raise ValueError(
                'Invalid observed-to-expected CpG ratio for ratio greater '
                'than or equal to zero: {0}'.format(min_obs_exp_cpg_ratio))


class SlidingWindowPythonAlgorithm(MetaAlgorithm):
    @property
    def name(self):
        return 'Sliding Window (Python)'

    def algorithm(
            self, seq_record, island_size, min_gc_ratio,
            min_obs_exp_cpg_ratio):
        super(SlidingWindowPythonAlgorithm, self).algorithm(
            seq_record, island_size, min_gc_ratio, min_obs_exp_cpg_ratio)

        seq_str = str(seq_record.seq)
        seq_len = len(seq_str)
        island_features = []
        island_metadata_list = []
        start_index = 0
        end_index = island_size
        gc_ratio = 0
        obs_exp_cpg_ratio = 0

        while end_index <= seq_len:
            # Keep adding bases to the end of the subsequence until
            # the subsequence no longer meets the criteria for being
            # an island.
            while end_index <= seq_len:
                c_count, g_count, cpg_count = _compute_counts(
                    seq_str[start_index:end_index])
                if c_count == 0 or g_count == 0:
                    break
                last_gc_ratio = gc_ratio
                last_obs_exp_cpg_ratio = obs_exp_cpg_ratio
                gc_ratio, obs_exp_cpg_ratio = _compute_ratios(
                    c_count, g_count, cpg_count, end_index - start_index)
                if not _is_island(gc_ratio, obs_exp_cpg_ratio,
                                  min_gc_ratio, min_obs_exp_cpg_ratio):
                    break
                end_index += 1
            # If `end_index' was incremented at least once, we've
            # found an island.
            if end_index > start_index + island_size:
                # This means that the true exclusive end index of the
                # island is one less than `end_index'.
                island_end_index = end_index - 1
                island_features.append(
                    _make_feature(start_index, island_end_index))
                # If we ended by reaching the end of the sequence, the
                # current ratios are the accurate ones. Otherwise, the
                # last ratios computed will be accurate.
                if end_index > seq_len:
                    island_metadata = IslandMetadata(gc_ratio,
                                                     obs_exp_cpg_ratio)
                else:
                    island_metadata = IslandMetadata(last_gc_ratio,
                                                     last_obs_exp_cpg_ratio)
                island_metadata_list.append(island_metadata)
                # Reset the pointer to the start of the subsequence to
                # the exclusive end of the island we just found.
                start_index = island_end_index
            else:
                # If we didn't find an island, increment the pointer
                # to the start of the subsequence.
                start_index += 1
            # Whether we found an island or not, `end_index' needs to
            # be reset to the start index plus the window size.
            end_index = start_index + island_size

        seq_record.features = island_features
        return AlgoResults(seq_record, island_metadata_list)


class AccumulatingSlidingWindowPythonAlgorithm(MetaAlgorithm):
    @property
    def name(self):
        return 'Accumulating Sliding Window (Python)'

    # NOQA is here right now to stop flake8 from whining about the
    # cyclomatic complexity of this function, which it reports as
    # 13. This should be fixed.
    def algorithm(  # NOQA
            self, seq_record, island_size, min_gc_ratio,
            min_obs_exp_cpg_ratio):
        super(AccumulatingSlidingWindowPythonAlgorithm, self).algorithm(
            seq_record, island_size, min_gc_ratio, min_obs_exp_cpg_ratio)

        island_features = []
        island_metadata_list = []
        seq_str = str(seq_record.seq)
        seq_len = len(seq_str)
        start_index = 0
        end_index = island_size
        # Calculate initial counts.
        c_count, g_count, cpg_count = _compute_counts(
            seq_str[start_index:end_index])
        is_island = False
        gc_ratio = 0
        obs_exp_cpg_ratio = 0

        while True:
            was_island = is_island
            # Calculate whether we have an island.
            if c_count == 0 or g_count == 0:
                is_island = False
            else:
                last_gc_ratio = gc_ratio
                last_obs_exp_cpg_ratio = obs_exp_cpg_ratio
                gc_ratio, obs_exp_cpg_ratio = _compute_ratios(
                    c_count, g_count, cpg_count, end_index - start_index)
                is_island = _is_island(gc_ratio, obs_exp_cpg_ratio,
                                       min_gc_ratio, min_obs_exp_cpg_ratio)

            if not is_island:
                if was_island:
                    # We were in an island and we just exited the
                    # island. That means we've found the largest
                    # island we can.

                    # The true exclusive end index of the island is
                    # one less than `end_index'.
                    island_end_index = end_index - 1
                    island_features.append(
                        _make_feature(start_index, island_end_index))
                    island_metadata_list.append(
                        IslandMetadata(last_gc_ratio, last_obs_exp_cpg_ratio))
                    # Reset the pointer to the start of the subsequence to
                    # the exclusive end of the island we just found.
                    start_index = island_end_index
                    # The pointer to the end of the subsequence should
                    # now be the minimal window.
                    end_index = start_index + island_size
                    # Recalculate initial counts.
                    c_count, g_count, cpg_count = \
                        _compute_counts(seq_str[start_index:end_index])
                    is_island = False
                    # Start again looking again.
                    continue
                else:
                    # We are not in and island and we weren't in an
                    # island. First calculate what we are going to
                    # lose from the subsequence and update the counts
                    # appropriately.
                    if seq_str[start_index] == 'C':
                        c_count -= 1
                        if seq_str[start_index + 1] == 'G':
                            cpg_count -= 1
                    elif seq_str[start_index] == 'G':
                        g_count -= 1
                    # Increment `start_index'.
                    start_index += 1
            # Increment `end_index'.
            end_index += 1
            # If `end_index' is greater than length of the sequence,
            # we have reached the end. Exit.
            if end_index > seq_len:
                break
            # `end_index - 1` now refers to the index of the last base
            # in the subsequence. Check to see what we have added to
            # the subsequence and update counts appropriately.
            if seq_str[end_index - 1] == 'C':
                c_count += 1
            elif seq_str[end_index - 1] == 'G':
                g_count += 1
                if seq_str[end_index - 2] == 'C':
                    cpg_count += 1

        # If we ended in an island, it won't be recorded by the
        # loop. Record it now. We also use the current ratios, since
        # the loop exited due to bounds checking.
        if is_island:
            island_features.append(_make_feature(start_index, seq_len))
            island_metadata_list.append(
                IslandMetadata(gc_ratio, obs_exp_cpg_ratio))

        seq_record.features = island_features
        return AlgoResults(seq_record, island_metadata_list)

# class AccumulatingSlidingWindowCythonAlgorithm(MetaAlgorithm):
#     @property
#     def name(self):
#         return 'Accumulating Sliding Window (Cython)'

#     def algorithm(self, seq_record, island_size, min_gc_ratio):
#         super(AccumulatingSlidingWindowCythonAlgorithm,
#               self).algorithm(seq_record, island_size, min_gc_ratio)
#         island_tuples = sliding_window_cython.sliding_window(
#             str(seq_record.seq), island_size, min_gc_ratio)
#         seq_record.features = \
#             [_make_feature(start, end) for start, end in island_tuples]
#         return seq_record


# Create instances of each implementation of a MetaAlgorithm, so that
# they are ready to use.
registry = [class_() for class_ in MetaAlgorithm.__subclasses__()]
