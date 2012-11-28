from __future__ import division

import pytest
from mock import MagicMock, sentinel, call

from cpg_islands.models import ResultsModel, IslandInfo
from cpg_islands.algorithms import AlgoResults, IslandMetadata
from tests.helpers import make_seq_record


@pytest.fixture
def model():
    return ResultsModel()


class TestResultsModel:
    def test_set_results_islands_computed_called(self, model):
        seq_str = 'ATATCGCGCGCGCATATA'
        feature_tuples = [(0, 3), (5, 7), (8, 13)]
        seq_record = make_seq_record(seq_str, feature_tuples)
        results = AlgoResults(seq_record, [])

        callback = MagicMock()
        model.islands_computed.append(callback)
        model.set_results(results, sentinel.algo_name, sentinel.exec_time)
        assert (callback.mock_calls ==
                [call(seq_str,
                      feature_tuples,
                      sentinel.algo_name,
                      sentinel.exec_time)])

    def test_get_island_info(self, model):
        seq_record = make_seq_record(
            'ATATCGCGCGCGCATATA', [(0, 3), (5, 7), (8, 13)])
        island_metadata_list = [
            IslandMetadata(0.57, 0.89),
            IslandMetadata(0.65, 2.13),
            IslandMetadata(0.78, 1.3)]
        results = AlgoResults(seq_record, island_metadata_list)

        model.set_results(results, sentinel.algo_name, sentinel.exec_time)
        computed = model.get_island_info(1)
        expected = IslandInfo(5, 7, 2, 'GC', 0.65, 2.13)
        assert computed == expected
