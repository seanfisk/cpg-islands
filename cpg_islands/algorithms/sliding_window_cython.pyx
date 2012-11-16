# -*- mode: python -*-

from __future__ import division

cdef int _is_gc(char base):
    """Test whether a specified base is Guanine or Cytosine.

    :param base: the base to test
    :return: whether the base is G or C
    """
    return base == 'G' or base == 'C'


cdef int _count_gc(char *partial_seq_str, int n):
    """Count the number of Guanine or Cytosine bases in a sequence.

    :param partial_seq_str: the partial sequence (this is not null-terminated)
    :param n: the number of bases to count
    :return: the GC count
    """
    cdef int gc_count = 0
    cdef int i
    for i in xrange(n):
        if _is_gc(partial_seq_str[i]):
            gc_count += 1
    return gc_count


def sliding_window(
        char *seq_str, int island_size, double min_gc_ratio):
    """Cython sliding window algorithm.

    :param seq_str: the sequence as a string
    :param island_size: the number of bases which an island may contain
    :param min_gc_ratio: the ratio of GC to other bases
    """
    # Might want to try passing in length of sequence string
    cdef int seq_len = len(seq_str)
    cdef int start_index = 0
    cdef int gc_count = _count_gc(seq_str, island_size)
    cdef int end_index
    islands = []
    while True:
        # check
        # note: end_index is exclusive, like Python's slice
        end_index = start_index + island_size
        if gc_count / island_size >= min_gc_ratio:
            islands.append((start_index, end_index))
        if end_index >= seq_len:
            break
        # update
        if _is_gc(seq_str[start_index]):
            gc_count -= 1
        if _is_gc(seq_str[end_index]):
            gc_count += 1
        start_index += 1
    return islands
