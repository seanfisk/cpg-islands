from __future__ import division

import pytest

from cpg_islands import algorithms
from tests.helpers import make_features, make_seq_record


def extract_features(seq_record):
    return [str(feature.extract(seq_record.seq)) for feature in
            seq_record.features]


def assert_seq_records_equal(computed, expected):
    assert extract_features(computed) == extract_features(expected)

# Use the old style parametrization code to set explicit ids. This
# makes the output look far better. To see why, uncomment the 2.3
# fixture code and comment the pytest_generate_tests function.

# @pytest.fixture(params=algorithms.registry)
# def algorithm(request):
#     return request.param


def pytest_generate_tests(metafunc):
    if 'algorithm' in metafunc.fixturenames:
        metafunc.parametrize(
            'algorithm',
            [instance.algorithm for instance in algorithms.registry],
            ids=[instance.id for instance in algorithms.registry])


class TestAlgorithms:
    def test_empty_sequence(self, algorithm):
        with pytest.raises(ValueError) as exc_info:
            algorithm(make_seq_record(''), 1, 0)
        assert (str(exc_info.value) ==
                'Island size (1) must be less than or '
                'equal to sequence length (0)')

    def test_zero_island_size(self, algorithm):
        with pytest.raises(ValueError) as exc_info:
            algorithm(make_seq_record(''), 0, 0)
        # exc_info.value returns the actual exception
        assert str(exc_info.value) == 'Invalid island size: 0'

    def test_negative_island_size(self, algorithm):
        with pytest.raises(ValueError) as exc_info:
            algorithm(make_seq_record(''), -1, 0)
        assert str(exc_info.value) == 'Invalid island size: -1'

    def test_island_size_less_than_sequence_size(self, algorithm):
        """When the user submits an island size greater than the
        sequence size, they are shown an error.
        """
        with pytest.raises(ValueError) as exc_info:
            algorithm(make_seq_record('ATATGCGC'), 9, 0.5)
        assert ((str(exc_info.value)) ==
                'Island size (9) must be less than or '
                'equal to sequence length (8)')

    def test_negative_gc_ratio(self, algorithm):
        """When the user submits a negative GC ratio, they are
        shown an error.
        """
        with pytest.raises(ValueError) as exc_info:
            algorithm(make_seq_record('ATGC'), 2, -1.5)
        assert (str(exc_info.value) ==
                'Invalid GC ratio for ratio between zero and one: -1.5')

    def test_greater_than_one_gc_ratio(self, algorithm):
        """When the user submits a GC ratio greater than one, they
        are shown an error.
        """
        with pytest.raises(ValueError) as exc_info:
            algorithm(make_seq_record('ATGC'), 2, 20)
        assert (str(exc_info.value) ==
                'Invalid GC ratio for ratio between zero and one: 20')

    def test_single_base(self, algorithm):
        seq_record = make_seq_record('C')
        computed = algorithm(seq_record, 1, 1)
        seq_record.features = make_features([(0, 1)])
        assert_seq_records_equal(computed, seq_record)

    def test_ratio_one(self, algorithm):
        seq_record = make_seq_record('ATGCCGATTTTA')
        computed = algorithm(seq_record, 4, 1)
        seq_record.features = make_features([(2, 6)])
        assert_seq_records_equal(computed, seq_record)

    def test_ratio_half(self, algorithm):
        seq_record = make_seq_record('ATATGCTAAT')
        computed = algorithm(seq_record, 4, 0.5)
        seq_record.features = make_features([(2, 6), (3, 7), (4, 8)])
        assert_seq_records_equal(computed, seq_record)

    def test_ratio_third(self, algorithm):
        seq_record = make_seq_record('TTATATGCTAATAT')
        size = 6
        computed = algorithm(seq_record, size, 1 / 3)
        seq_record.features = \
            make_features([(i, i + size) for i in xrange(2, 7)])
        assert_seq_records_equal(computed, seq_record)

    def test_island_at_end(self, algorithm):
        seq_record = make_seq_record('ATATATGCGC')
        computed = algorithm(seq_record, 4, 1)
        seq_record.features = make_features([(6, 10)])
        assert_seq_records_equal(computed, seq_record)
