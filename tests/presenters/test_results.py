from mock import create_autospec, call, sentinel
import pytest

from cpg_islands.models import MetaResultsModel, IslandInfo
from cpg_islands.views import BaseResultsView
from cpg_islands.presenters import ResultsPresenter


@pytest.fixture
def presenter():
    mock_model = create_autospec(MetaResultsModel, spec_set=True)
    mock_view = create_autospec(BaseResultsView, spec_set=True)
    presenter = ResultsPresenter(mock_model, mock_view)
    return presenter


class TestResultPresenter:
    def test_register_for_events(self, presenter):
        presenter.register_for_events()
        assert (presenter.model.mock_calls ==
                [call.islands_computed.append(
                    presenter._islands_computed)])
        assert (presenter.view.mock_calls ==
                [call.island_selected.append(presenter._island_selected)])

    def test_islands_computed(self, presenter):
        presenter._islands_computed(sentinel.global_seq,
                                    sentinel.feature_tuples,
                                    sentinel.algo_name,
                                    134.45)
        assert (presenter.view.mock_calls ==
                [call.clear_subseq(),
                 call.set_global_seq(sentinel.global_seq),
                 call.set_islands(sentinel.feature_tuples),
                 call.set_algo_name(sentinel.algo_name),
                 call.set_exec_time('134.45 seconds')
                 ])
        assert presenter.model.mock_calls == []

    def test_island_selected(self, presenter):
        presenter.model.get_island_info.return_value = \
            IslandInfo(45, 67, 22, 'ATGCTA', 0.567, 0.845)
        presenter._island_selected(1)
        assert presenter.model.mock_calls == [call.get_island_info(1)]
        assert (presenter.view.mock_calls ==
                [call.highlight_global_seq(45, 67),
                 call.set_start('45 (zero-indexed)'),
                 call.set_end('67 (zero-indexed)'),
                 call.set_length('22 bases'),
                 call.set_subseq('ATGCTA'),
                 call.set_gc_ratio('56.7%'),
                 call.set_obs_exp_cpg_ratio('0.845')])
