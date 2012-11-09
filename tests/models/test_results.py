from __future__ import division

import pytest
from mock import MagicMock
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation

from cpg_islands.models import ResultsModel
from tests.helpers import make_features


def extract_features(feature_list, sequence):
    return [str(feature.extract(sequence)) for feature in feature_list]


def assert_features_equal(computed_features, expected_features, sequence):
    computed_features_extracted = extract_features(computed_features, sequence)
    expected_features_extracted = extract_features(expected_features, sequence)
    assert computed_features_extracted == expected_features_extracted


@pytest.fixture
def model():
    return ResultsModel()


class TestResultsModel:
    class TestAnnotate:
        def test_empty_sequence(self, model):
            with pytest.raises(ValueError) as exc_info:
                model.annotate_cpg_islands(Seq('', IUPAC.unambiguous_dna),
                                           1, 0)
                assert (str(exc_info.value) ==
                        'Island size (1) must be less than or '
                        'equal to sequence length (0)')

        def test_zero_island_size(self, model):
            with pytest.raises(ValueError) as exc_info:
                model.annotate_cpg_islands(Seq('', IUPAC.unambiguous_dna),
                                           0, 0)
            # exc_info.value returns the actual exception
            assert str(exc_info.value) == 'Invalid island size: 0'

        def test_negative_island_size(self, model):
            with pytest.raises(ValueError) as exc_info:
                model.annotate_cpg_islands(Seq('', IUPAC.unambiguous_dna),
                                           -1, 0)
            assert str(exc_info.value) == 'Invalid island size: -1'

        def test_island_size_less_than_sequence_size(self, model):
            """When the user submits an island size greater than the
            sequence size, they are shown an error.
            """
            with pytest.raises(ValueError) as exc_info:
                model.annotate_cpg_islands(Seq('ATATGCGC',
                                               IUPAC.unambiguous_dna),
                                           9, 0.5)
                assert ((str(exc_info.value)) ==
                        'Island size (9) must be less than or '
                        'equal to sequence length (8)')

        def test_negative_gc_ratio(self, model):
            """When the user submits a negative GC ratio, they are
            shown an error.
            """
            with pytest.raises(ValueError) as exc_info:
                model.annotate_cpg_islands(Seq('ATGC', IUPAC.unambiguous_dna),
                                           2, -1.5)
            assert (str(exc_info.value) ==
                    'Invalid GC ratio for ratio between zero and one: -1.5')

        def test_greater_than_one_gc_ratio(self, model):
            """When the user submits a GC ratio greater than one, they
            are shown an error.
            """
            with pytest.raises(ValueError) as exc_info:
                model.annotate_cpg_islands(Seq('ATGC', IUPAC.unambiguous_dna),
                                           2, 20)
            assert (str(exc_info.value) ==
                    'Invalid GC ratio for ratio between zero and one: 20')

        def test_single_base(self, model):
            locations_computed_callback = MagicMock()
            model.locations_computed.append(locations_computed_callback)
            seq = Seq('C', IUPAC.unambiguous_dna)
            model.annotate_cpg_islands(seq, 1, 1)
            assert locations_computed_callback.call_count == 1
            args = locations_computed_callback.call_args[0]
            assert len(args) == 1
            assert_features_equal(
                args[0],
                make_features([(0, 1)]),
                seq)

        def test_ratio_one(self, model):
            locations_computed_callback = MagicMock()
            model.locations_computed.append(locations_computed_callback)
            seq = Seq('ATGCCGATTTTA', IUPAC.unambiguous_dna)
            model.annotate_cpg_islands(seq, 4, 1)
            assert locations_computed_callback.call_count == 1
            args = locations_computed_callback.call_args[0]
            assert len(args) == 1
            assert_features_equal(
                args[0],
                make_features([(2, 6)]),
                seq)

        def test_ratio_half(self, model):
            locations_computed_callback = MagicMock()
            model.locations_computed.append(locations_computed_callback)
            seq = Seq('ATATGCTAAT', IUPAC.unambiguous_dna)
            model.annotate_cpg_islands(seq, 4, 0.5)
            assert locations_computed_callback.call_count == 1
            args = locations_computed_callback.call_args[0]
            assert len(args) == 1
            assert_features_equal(
                args[0],
                make_features([(2, 6), (3, 7), (4, 8)]),
                seq)

        def test_ratio_third(self, model):
            locations_computed_callback = MagicMock()
            model.locations_computed.append(locations_computed_callback)
            seq = Seq('TTATATGCTAATAT', IUPAC.unambiguous_dna)
            size = 6
            model.annotate_cpg_islands(seq, size, 1 / 3)
            assert locations_computed_callback.call_count == 1
            args = locations_computed_callback.call_args[0]
            assert len(args) == 1
            assert_features_equal(
                args[0],
                make_features([(i, i + size) for i in xrange(2, 7)]),
                seq)

        def test_seq_set(self, model):
            seq = Seq('TTATATGCTAATAT', IUPAC.unambiguous_dna)
            model.annotate_cpg_islands(seq, 6, 1 / 3)
            assert model._seq == seq

        def test_features_set(self, model):
            seq = Seq('ATATGCTAAT', IUPAC.unambiguous_dna)
            model.annotate_cpg_islands(seq, 4, 0.5)
            assert_features_equal(
                model._features,
                make_features([(2, 6), (3, 7), (4, 8)]),
                seq)

    class TestGetLocalSeq:
        def test_get_one(self, model):
            model._seq = Seq('GCTT', IUPAC.unambiguous_dna)
            model._features = [SeqFeature(FeatureLocation(0, 1))]
            computed = model.get_local_seq(0)
            assert computed == 'G'

    class TestGetGlobalSeq:
        def test_get_one(self, model):
            model._seq = Seq('GCTT', IUPAC.unambiguous_dna)
            computed = model.get_global_seq()
            assert computed == 'GCTT'
