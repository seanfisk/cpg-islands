from __future__ import division

import pytest
from mock import MagicMock, sentinel, call

from cpg_islands.models import ResultsModel


@pytest.fixture
def model():
    return ResultsModel()


class TestResultsModel:
    def test_set_results(self, model):
        callback = MagicMock()
        model.islands_computed.append(callback)
        model.set_results(
            sentinel.seq_record, sentinel.algo_name, sentinel.exec_time)
        assert (callback.mock_calls ==
                [call(sentinel.seq_record,
                      sentinel.algo_name, sentinel.exec_time)])

    def test_get_seq_record(self, model):
        model.set_results(
            sentinel.seq_record, sentinel.algo_name, sentinel.exec_time)
        assert model.get_seq_record() == sentinel.seq_record
