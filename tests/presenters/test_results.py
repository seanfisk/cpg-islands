from mock import create_autospec, call, sentinel
import pytest

from cpg_islands.models import MetaResultsModel
from cpg_islands.views import BaseResultsView
from cpg_islands.presenters import ResultsPresenter
from tests.helpers import make_features


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
        island_tuples = [(0, 5), (1, 6), (3, 8)]
        islands = make_features(island_tuples)
        presenter._islands_computed(islands)
        assert (presenter.view.mock_calls ==
                [call.set_islands(island_tuples)])
        assert presenter.model.mock_calls == []

    def test_island_selected(self, presenter):
        presenter.model.get_global_seq.return_value = sentinel.global_seq
        presenter.model.get_local_seq.return_value = sentinel.local_seq
        presenter._island_selected(sentinel.island_index)
        assert (presenter.model.mock_calls ==
                [call.get_global_seq(),
                 call.get_local_seq(sentinel.island_index)])
        assert (presenter.view.mock_calls ==
                [call.set_global_seq(sentinel.global_seq),
                 call.set_local_seq(sentinel.local_seq)])
