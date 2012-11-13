from __future__ import division

import pytest
from mock import MagicMock, sentinel, call
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation

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

    class TestGetLocalSeq:
        def test_get_one(self, model):
            model._seq = Seq('GCTT', IUPAC.unambiguous_dna)
            model._islands = [SeqFeature(FeatureLocation(0, 1))]
            computed = model.get_local_seq(0)
            assert computed == 'G'

    class TestGetGlobalSeq:
        def test_get_one(self, model):
            model._seq = Seq('GCTT', IUPAC.unambiguous_dna)
            computed = model.get_global_seq()
            assert computed == 'GCTT'
