from __future__ import division

import pytest

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from cpg_islands import algorithms
from tests.helpers import make_features


def extract_features(feature_list, sequence):
    return [str(feature.extract(sequence)) for feature in feature_list]


def assert_features_equal(computed_features, expected_features, sequence):
    computed_features_extracted = extract_features(computed_features, sequence)
    expected_features_extracted = extract_features(expected_features, sequence)
    assert computed_features_extracted == expected_features_extracted


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
            algorithm(Seq('', IUPAC.unambiguous_dna),
                      1, 0)
        assert (str(exc_info.value) ==
                'Island size (1) must be less than or '
                'equal to sequence length (0)')

    def test_zero_island_size(self, algorithm):
        with pytest.raises(ValueError) as exc_info:
            algorithm(Seq('', IUPAC.unambiguous_dna),
                      0, 0)
        # exc_info.value returns the actual exception
        assert str(exc_info.value) == 'Invalid island size: 0'

    def test_negative_island_size(self, algorithm):
        with pytest.raises(ValueError) as exc_info:
            algorithm(Seq('', IUPAC.unambiguous_dna),
                      -1, 0)
        assert str(exc_info.value) == 'Invalid island size: -1'

    def test_island_size_less_than_sequence_size(self, algorithm):
        """When the user submits an island size greater than the
        sequence size, they are shown an error.
        """
        with pytest.raises(ValueError) as exc_info:
            algorithm(Seq('ATATGCGC',
                          IUPAC.unambiguous_dna),
                      9, 0.5)
        assert ((str(exc_info.value)) ==
                'Island size (9) must be less than or '
                'equal to sequence length (8)')

    def test_negative_gc_ratio(self, algorithm):
        """When the user submits a negative GC ratio, they are
        shown an error.
        """
        with pytest.raises(ValueError) as exc_info:
            algorithm(Seq('ATGC', IUPAC.unambiguous_dna),
                      2, -1.5)
        assert (str(exc_info.value) ==
                'Invalid GC ratio for ratio between zero and one: -1.5')

    def test_greater_than_one_gc_ratio(self, algorithm):
        """When the user submits a GC ratio greater than one, they
        are shown an error.
        """
        with pytest.raises(ValueError) as exc_info:
            algorithm(Seq('ATGC', IUPAC.unambiguous_dna),
                      2, 20)
        assert (str(exc_info.value) ==
                'Invalid GC ratio for ratio between zero and one: 20')

    def test_single_base(self, algorithm):
        seq = Seq('C', IUPAC.unambiguous_dna)
        computed = algorithm(seq, 1, 1)
        assert_features_equal(
            computed,
            make_features([(0, 1)]),
            seq)

    def test_ratio_one(self, algorithm):
        seq = Seq('ATGCCGATTTTA', IUPAC.unambiguous_dna)
        computed = algorithm(seq, 4, 1)
        assert_features_equal(
            computed,
            make_features([(2, 6)]),
            seq)

    def test_ratio_half(self, algorithm):
        seq = Seq('ATATGCTAAT', IUPAC.unambiguous_dna)
        computed = algorithm(seq, 4, 0.5)
        assert_features_equal(
            computed,
            make_features([(2, 6), (3, 7), (4, 8)]),
            seq)

    def test_ratio_third(self, algorithm):
        seq = Seq('TTATATGCTAATAT', IUPAC.unambiguous_dna)
        size = 6
        computed = algorithm(seq, size, 1 / 3)
        assert_features_equal(
            computed,
            make_features([(i, i + size) for i in xrange(2, 7)]),
            seq)

    def test_island_at_end(self, algorithm):
        seq = Seq('ATATATGCGC', IUPAC.unambiguous_dna)
        computed = algorithm(seq, 4, 1)
        assert_features_equal(
            computed,
            make_features([(6, 10)]),
            seq)
