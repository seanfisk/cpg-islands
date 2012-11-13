from __future__ import division

import pytest
from mock import MagicMock, sentinel, call
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC

from cpg_islands.models import ResultsModel


@pytest.fixture
def model():
    return ResultsModel()


class TestResultsModel:
    def test_set_results(self, model):
        callback = MagicMock()
        model.islands_computed.append(callback)
        model.set_results(sentinel.seq, sentinel.islands)
        assert callback.mock_calls == [call(sentinel.islands)]

    # TODO: These tests need to improve
    class TestGetLocalSeq:
        def test_get_one(self, model):
            model._seq = Seq('GCTT', IUPAC.unambiguous_dna)
            model._islands = [SeqFeature(FeatureLocation(0, 1))]
            assert model.get_local_seq(0) == 'G'

    # TODO: These tests need to improve
    class TestGetGlobalSeq:
        def test_get_one(self, model):
            model._seq = Seq('GCTT', IUPAC.unambiguous_dna)
            assert model.get_global_seq() == 'GCTT'
