from mock import create_autospec, call
import pytest

from cpg_islands.models import MetaResultsModel
from cpg_islands.views import BaseResultsView
from cpg_islands.presenters import ResultsPresenter
from tests.helpers import make_seq_record


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
        seq_record = make_seq_record('', island_tuples)
        presenter._islands_computed(seq_record)
        assert (presenter.view.mock_calls ==
                [call.clear_local_seq(),
                 call.clear_global_seq(),
                 call.set_islands(island_tuples),
                 ])
        assert presenter.model.mock_calls == []

    def test_island_selected(self, presenter):
        seq_str = 'ATGCGCAT'
        presenter.model.get_results.return_value = \
            make_seq_record(seq_str, [(2, 4), (3, 5), (4, 6)])
        presenter._island_selected(1)
        assert presenter.model.mock_calls == [call.get_results()]
        assert (presenter.view.mock_calls ==
                [call.set_global_seq(seq_str, (3, 5)),
                 call.set_local_seq('CG')])
