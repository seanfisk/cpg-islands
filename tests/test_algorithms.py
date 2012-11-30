from __future__ import division

import pytest

from cpg_islands import algorithms
from tests.helpers import make_algo_results, make_seq_record


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
            algorithm(make_seq_record(''), 1, 0, 0)
        assert (str(exc_info.value) ==
                'Island size (1) must be less than or '
                'equal to sequence length (0)')

    def test_zero_island_size(self, algorithm):
        with pytest.raises(ValueError) as exc_info:
            algorithm(make_seq_record(''), 0, 0, 0)
        # exc_info.value returns the actual exception
        assert str(exc_info.value) == 'Invalid island size: 0'

    def test_negative_island_size(self, algorithm):
        with pytest.raises(ValueError) as exc_info:
            algorithm(make_seq_record(''), -1, 0, 0.6)
        assert str(exc_info.value) == 'Invalid island size: -1'

    def test_island_size_less_than_sequence_size(self, algorithm):
        """When the user submits an island size greater than the
        sequence size, they are shown an error.
        """
        with pytest.raises(ValueError) as exc_info:
            algorithm(make_seq_record('ATATGCGC'), 9, 0.5, 0.6)
        assert ((str(exc_info.value)) ==
                'Island size (9) must be less than or '
                'equal to sequence length (8)')

    def test_negative_gc_ratio(self, algorithm):
        """When the user submits a negative GC ratio, they are
        shown an error.
        """
        with pytest.raises(ValueError) as exc_info:
            algorithm(make_seq_record('ATGC'), 2, -1.5, 0.6)
        assert (str(exc_info.value) ==
                'Invalid GC ratio for ratio between zero and one: -1.5')

    def test_zero_gc_ratio(self, algorithm):
        """The user can submit a GC ratio of zero."""
        algorithm(make_seq_record('ATCG'), 2, 0, 0.6)

    def test_one_gc_ratio(self, algorithm):
        """The user can submit a GC ratio of one."""
        algorithm(make_seq_record('ATCG'), 2, 1, 0.6)

    def test_greater_than_one_gc_ratio(self, algorithm):
        """When the user submits a GC ratio greater than one, they
        are shown an error.
        """
        with pytest.raises(ValueError) as exc_info:
            algorithm(make_seq_record('ATGC'), 2, 20, 0.6)
        assert (str(exc_info.value) ==
                'Invalid GC ratio for ratio between zero and one: 20')

    def test_negative_obs_to_exp_ratio(self, algorithm):
        """When the user submits a negative observed-to-expected CpG
        ratio, they are shown an error.
        """
        with pytest.raises(ValueError) as exc_info:
            algorithm(make_seq_record('ATGC'), 2, 0.5, -1.5)
        assert (str(exc_info.value) ==
                'Invalid observed-to-expected CpG ratio '
                'for ratio greater than or equal to zero: -1.5')

    def test_zero_obs_to_exp_ratio(self, algorithm):
        """The user can submit an observed-to-expected CpG ratio of zero."""
        algorithm(make_seq_record('ATGC'), 2, 0.5, 0)

    def test_greater_than_one_obs_to_exp_ratio(self, algorithm):
        """The user can submit an observed-to-expected CpG ratio
        greater than one.
        """
        algorithm(make_seq_record('ATGC'), 2, 0.5, 1.5)

    def test_single_cpg(self, algorithm):
        seq_str = 'CG'
        computed = algorithm(make_seq_record(seq_str), 2, 1, 2)
        expected = make_algo_results(seq_str, [(0, 2, 1, 2)])
        print computed.seq_record.features[0].extract(computed.seq_record.seq)
        for a in [computed, expected]:
            print a.island_metadata_list[0].gc_ratio
        assert computed == expected

    def test_island_at_beginning(self, algorithm):
        seq_str = 'CGGATATATA'
        computed = algorithm(make_seq_record(seq_str), 3, 0.5, 0.6)
        expected = make_algo_results(seq_str, [(0, 6, 0.5, 3)])
        assert computed == expected

    def test_island_in_middle(self, algorithm):
        seq_str = 'ATATACACGGAATATT'
        computed = algorithm(make_seq_record(seq_str), 4, 0.5, 0.6)
        expected = make_algo_results(seq_str, [(5, 13, 0.5, 2)])
        assert computed == expected

    def test_island_at_end(self, algorithm):
        seq_str = 'ATATATTATTCAACGAGG'
        computed = algorithm(make_seq_record(seq_str), 5, 0.5, 0.6)
        expected = make_algo_results(seq_str, [(10, 18, 0.625,  4 / 3)])
        assert computed == expected
